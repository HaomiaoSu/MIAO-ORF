#!/usr/bin/env python3
"""Top-level runner for the MIAO analysis workflow.

This module deliberately remains a thin orchestrator.  The scientific logic
continues to live in the stage programs under ``src/miao_orf``.
"""

from __future__ import annotations

import argparse
import copy
import csv
import datetime as dt
import hashlib
import importlib.util
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


PROGRAM = "MIAO"
__version__ = "1.0.0"
COMPLETION_SCHEMA_VERSION = "1.0"
IDENTITY_HASH_LIMIT_BYTES = 64 * 1024 * 1024
STAGES = (
    "orfscan", "psite", "qc", "dm", "abundance", "context", "codon", "visualize",
)
MODE_SPECIFIC_STAGES = {"dm", "abundance", "context", "codon", "visualize"}
STAGE_TITLES = {
    "orfscan": "Scan transcriptome ORFs",
    "psite": "Build 1-nt P-site BAM",
    "qc": "Metagene QC and DM background",
    "dm": "Call intORFs with the DM model",
    "abundance": "Quantify model-allocated intORF abundance",
    "context": "Annotate same-gene annotated-CDS N-terminal reuse",
    "codon": "Export selected-ORF codon-level P-site data (optional internal tool)",
    "visualize": "Visualize final gene-level pure intORF results",
}


@dataclass(frozen=True)
class Paths:
    project: Path
    source: Path
    out_root: Path
    reference_name: str
    sample: str
    orf_prefix: Path
    per_chrom_dir: Path
    psite_prefix: Path
    qc_prefix: Path
    dm_prefix: Path
    abundance_prefix: Path
    context_prefix: Path
    codon_prefix: Path
    visualization_prefix: Path
    torf: Path
    psite_bam: Path
    dm_background: Path
    dm_results: Path
    dm_run_manifest: Path
    gene_context_results: Path


@dataclass
class StageSpec:
    name: str
    command: List[str]
    required_inputs: List[Path]
    expected_outputs: List[Path]
    log_path: Path


@dataclass(frozen=True)
class BatchReplicate:
    replicate_id: str
    bam: Path
    input_mode: str
    length_offsets: Optional[Dict[int, int]]
    ribotish_para: Optional[Path]
    ribotish_quality: Optional[Path]


@dataclass
class PipelineOutcome:
    sample: str
    input_mode: str
    status: str
    error: str
    failed_stage: str
    manifest: Optional[Path]
    mode_paths: Dict[str, Paths]


@dataclass
class IntegrationOutcome:
    dm_mode: str
    status: str
    error: str
    command: List[str]
    log: Path
    outputs: List[Path]
    visualization_status: str = ""
    visualization_error: str = ""
    visualization_command: Optional[List[str]] = None
    visualization_log: Optional[Path] = None
    visualization_outputs: Optional[List[Path]] = None
    context_status: str = ""
    context_error: str = ""
    context_command: Optional[List[str]] = None
    context_log: Optional[Path] = None
    context_outputs: Optional[List[Path]] = None


def csv_ints(value: str) -> List[int]:
    try:
        result = [int(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from exc
    if not result or any(item <= 0 for item in result):
        raise argparse.ArgumentTypeError("read lengths must be positive integers")
    return result


def csv_length_offsets(value: str) -> Dict[int, int]:
    result: Dict[int, int] = {}
    for item in value.split(","):
        token = item.strip()
        if not token:
            continue
        parts = token.split(":")
        if len(parts) != 2:
            raise argparse.ArgumentTypeError(
                "expected comma-separated LENGTH:OFFSET pairs, for example 28:12,29:12,30:12"
            )
        try:
            length, offset = (int(part.strip()) for part in parts)
        except ValueError as exc:
            raise argparse.ArgumentTypeError("lengths and offsets must be integers") from exc
        if length <= 0 or not 0 <= offset < length:
            raise argparse.ArgumentTypeError(
                f"invalid pair {token!r}; require length > 0 and 0 <= offset < length"
            )
        if length in result:
            raise argparse.ArgumentTypeError(f"duplicate read length: {length}")
        result[length] = offset
    if not result:
        raise argparse.ArgumentTypeError("at least one LENGTH:OFFSET pair is required")
    return result


def resolve_batch_path(value: str, batch_dir: Path) -> Path:
    expanded = Path(os.path.expandvars(os.path.expanduser(value.strip())))
    return expanded if expanded.is_absolute() else batch_dir / expanded


def read_batch_replicates(path: str | Path) -> List[BatchReplicate]:
    """Read and validate a replicate TSV for independent sample-level analysis."""
    batch_path = Path(path)
    required = {
        "replicate_id", "bam", "input_mode", "length_offsets",
        "ribotish_para", "ribotish_quality",
    }
    with batch_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = sorted(required - fields)
        if missing:
            raise ValueError(f"batch TSV is missing column(s): {', '.join(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError("batch TSV contains no replicate rows")

    result: List[BatchReplicate] = []
    seen: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        replicate_id = str(row.get("replicate_id", "")).strip()
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", replicate_id):
            raise ValueError(
                f"batch row {row_number}: replicate_id must use only letters, digits, '.', '_' or '-'"
            )
        if replicate_id.casefold() in seen:
            raise ValueError(f"batch row {row_number}: duplicate replicate_id {replicate_id!r}")
        seen.add(replicate_id.casefold())

        bam_value = str(row.get("bam", "")).strip()
        if not bam_value:
            raise ValueError(f"batch row {row_number}: bam is required")
        bam = resolve_batch_path(bam_value, batch_path.parent)
        input_mode = str(row.get("input_mode", "")).strip().casefold()
        length_text = str(row.get("length_offsets", "")).strip()
        para_text = str(row.get("ribotish_para", "")).strip()
        quality_text = str(row.get("ribotish_quality", "")).strip()

        if input_mode == "explicit":
            if not length_text:
                raise ValueError(
                    f"batch row {row_number}: explicit mode requires length_offsets"
                )
            if para_text or quality_text:
                raise ValueError(
                    f"batch row {row_number}: explicit mode cannot include Ribo-TISH files"
                )
            try:
                length_offsets = csv_length_offsets(length_text)
            except argparse.ArgumentTypeError as exc:
                raise ValueError(
                    f"batch row {row_number}: invalid length_offsets: {exc}"
                ) from exc
            para = quality = None
        elif input_mode == "ribotish":
            if length_text:
                raise ValueError(
                    f"batch row {row_number}: ribotish mode cannot include length_offsets"
                )
            if not para_text or not quality_text:
                raise ValueError(
                    f"batch row {row_number}: ribotish mode requires ribotish_para and ribotish_quality"
                )
            length_offsets = None
            para = resolve_batch_path(para_text, batch_path.parent)
            quality = resolve_batch_path(quality_text, batch_path.parent)
        else:
            raise ValueError(
                f"batch row {row_number}: input_mode must be 'explicit' or 'ribotish'"
            )

        result.append(
            BatchReplicate(
                replicate_id=replicate_id,
                bam=bam,
                input_mode=input_mode,
                length_offsets=length_offsets,
                ribotish_para=para,
                ribotish_quality=quality,
            )
        )
    return result


def stage_index(value: str) -> int:
    try:
        return STAGES.index(value)
    except ValueError as exc:
        raise ValueError(f"unknown stage: {value}") from exc


def path_arg(value: Optional[str]) -> Optional[Path]:
    if not value:
        return None
    return Path(os.path.expandvars(os.path.expanduser(value)))


def infer_reference_name(gtf: Optional[Path]) -> str:
    if gtf is None:
        return "reference"
    name = gtf.name
    for suffix in (".annotation.gtf.gz", ".annotation.gtf", ".gtf.gz", ".gtf"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return gtf.stem


def resolve_stage_source(project: Path) -> Path:
    source = project / "src" / "miao_orf"
    if source.is_dir():
        return source
    package_spec = importlib.util.find_spec("miao_orf")
    if package_spec is None or package_spec.origin is None:
        raise RuntimeError(
            "cannot locate installed miao_orf stage programs; reinstall the package"
        )
    return Path(package_spec.origin).resolve().parent


def add_boolean_option(
    parser: argparse.ArgumentParser,
    name: str,
    default: bool,
    help_text: str,
) -> None:
    parser.add_argument(
        name,
        action=argparse.BooleanOptionalAction,
        default=default,
        help=help_text,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="miao-orf",
        description=(
            "Run any contiguous portion of the MIAO workflow. "
            "Scientific stage programs remain "
            "independently runnable."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    identity = parser.add_mutually_exclusive_group(required=True)
    identity.add_argument("--sample", help="Single-sample label used in output paths")
    identity.add_argument(
        "--batch",
        help="Replicate TSV for sequential, failure-isolated batch processing",
    )
    parser.add_argument(
        "--out-root",
        required=True,
        help=(
            "Shared analysis root for reusable reference outputs and one or more "
            "samples; do not append the --sample name"
        ),
    )
    parser.add_argument("--from-stage", choices=STAGES, default="orfscan")
    parser.add_argument("--to-stage", choices=STAGES, default="visualize")
    parser.add_argument(
        "--only-stage",
        choices=STAGES,
        help="Run exactly one stage; overrides --from-stage and --to-stage",
    )
    parser.add_argument("--python", default=sys.executable, help="Python interpreter used for stage scripts")
    parser.add_argument("--workers", type=int, default=8, help="Default worker count for stages")
    parser.add_argument("--dry-run", action="store_true", help="Validate external inputs and print commands only")
    parser.add_argument("--force", action="store_true", help="Rerun selected stages even when outputs are complete")
    parser.add_argument("--status", action="store_true", help="Show expected primary outputs without running stages")

    inputs = parser.add_argument_group("input files and starting-point overrides")
    inputs.add_argument("--gtf", help="GENCODE/Ensembl GTF; required when orfscan is selected")
    inputs.add_argument("--fa", help="Indexed genome FASTA; required when orfscan is selected")
    inputs.add_argument("--bam", help="Sorted/indexed mapped Ribo-seq BAM; required when psite is selected")
    inputs.add_argument("--offsets", help="Automatic mode: Ribo-TISH mapped.para.py containing offdict")
    inputs.add_argument(
        "--ribotish-quality",
        help=(
            "Ribo-TISH mapped_qual.txt used to select 5'-matched read lengths with "
            "offset-corrected frame-0 proportion > 2/3"
        ),
    )
    inputs.add_argument(
        "--length-offsets", type=csv_length_offsets,
        help=(
            "Traditional mode: comma-separated read-length/P-site-offset pairs, "
            "for example 28:12,29:12,30:12"
        ),
    )
    inputs.add_argument("--torf", help="Existing tORF TSV when starting after orfscan")
    inputs.add_argument("--psite-bam", help="Existing 1-nt P-site BAM when starting after psite")
    inputs.add_argument("--dm-background", help="Existing validated DM background when starting at dm")
    inputs.add_argument(
        "--dm-results",
        help="Existing DM result TSV when starting at abundance or context",
    )
    inputs.add_argument(
        "--dm-run-manifest",
        help="Matching DM run manifest when starting directly at the codon stage",
    )
    inputs.add_argument(
        "--gene-context-results",
        help=(
            "Existing gene-context annotated result TSV when running only the "
            "codon or visualize stage"
        ),
    )
    inputs.add_argument("--reference-name", help="Reference label used in the default orfscan output path")

    outputs = parser.add_argument_group("optional output-prefix overrides")
    outputs.add_argument("--orfscan-out-prefix")
    outputs.add_argument("--psite-out-prefix")
    outputs.add_argument("--qc-out-prefix")
    outputs.add_argument("--dm-out-prefix")
    outputs.add_argument("--abundance-out-prefix")
    outputs.add_argument("--context-out-prefix")
    outputs.add_argument("--codon-out-prefix")
    outputs.add_argument("--visualization-out-prefix")

    scan = parser.add_argument_group("orfscan options")
    scan.add_argument("--min-aa", type=int, default=6)
    scan.add_argument("--start-codons", nargs="+", default=["ATG"])
    scan.add_argument("--near-cognate", action="store_true")
    add_boolean_option(scan, "--primary-only", True, "Restrict scanning to primary assembly contigs")
    add_boolean_option(scan, "--by-chrom", True, "Use chromosome-partitioned scanning")
    scan.add_argument("--orfscan-mp-chunksize", type=int, default=50)

    psite = parser.add_argument_group("P-site options")
    psite.add_argument("--min-mapq", type=int, default=20)
    psite.add_argument(
        "--min-frame0-proportion", type=float, default=2 / 3,
        help="Automatic mode: strict lower bound for corrected frame-0 proportion",
    )
    psite.add_argument(
        "--length-selection-policy",
        choices=("dominant_contiguous", "all_passing"),
        default="dominant_contiguous",
        help="Automatic mode: keep the read-richest contiguous block or all passing lengths",
    )
    psite.add_argument(
        "--keep-lengths", type=csv_ints,
        help=(
            "Optional comma-separated whitelist after Ribo-TISH quality screening; "
            "without --ribotish-quality this is an explicit legacy override"
        ),
    )
    add_boolean_option(psite, "--require-unique", True, "Require NH=1 when the tag is present")
    psite.add_argument("--bgzf-threads", type=int, default=2)
    psite.add_argument("--emit-bed", action="store_true")
    psite.add_argument("--emit-bedgraph", action="store_true")

    integration = parser.add_argument_group("replicate integration options")
    integration.add_argument(
        "--integration-min-replicates", type=int, default=2,
        help="Minimum primary-credible replicate count for reproducible consensus",
    )
    integration.add_argument(
        "--integration-min-fraction", type=float, default=0.5,
        help="Minimum primary-credible fraction among available replicate results",
    )

    model = parser.add_argument_group("QC and DM model options")
    model.add_argument("--pi-method", choices=("codon_equal", "trimmed", "read_weighted"), default="codon_equal")
    model.add_argument("--min-a0", type=float, default=1.0)
    model.add_argument("--min-intorf-aa", type=int, default=6)
    model.add_argument("--primary-min-intorf-aa", type=int, default=10)
    model.add_argument("--min-core-codons", type=int, default=5)
    model.add_argument("--min-active-core-codons", type=int, default=3)
    model.add_argument("--min-core-reads", type=int, default=15)
    model.add_argument("--min-credible-active-core-codons", type=int, default=5)
    model.add_argument("--min-credible-active-core-frac", type=float, default=0.15)
    model.add_argument("--min-credible-target-residual-frac", default="1/3")
    model.add_argument("--fdr-threshold", type=float, default=0.05)
    model.add_argument("--lambda-min", type=float, default=0.05)
    model.add_argument("--lambda-abs-diff-max", type=float, default=0.10)
    model.add_argument("--lambda-rel-diff-max", type=float, default=0.30)
    model.add_argument("--distance-to-segment-max", type=float, default=0.10)

    inference = parser.add_argument_group("formal DM inference options")
    inference.add_argument(
        "--dm-mode",
        choices=("accurate", "fast", "both"),
        default="accurate",
        help=(
            "Formal DM execution mode. 'both' runs shared upstream stages once, "
            "then writes independent accurate and fast DM/visualization outputs"
        ),
    )
    inference.add_argument("--bootstrap-reps", type=int, default=999)
    inference.add_argument("--importance-reps", type=int, default=9999)
    inference.add_argument("--seed", type=int, default=20260821)
    inference.add_argument("--dm-mp-chunksize", type=int, default=1)

    visual = parser.add_argument_group("visualization options")
    visual.add_argument("--formats", default="png,pdf")
    visual.add_argument("--dpi", type=int, default=220)
    visual.add_argument("--include-diagnostics", action="store_true")

    internal = parser.add_argument_group("optional internal ORF P-site export")
    internal.add_argument(
        "--export-orf-psites",
        action="store_true",
        help=(
            "Enable the normally disabled codon stage; without explicit IDs it exports "
            "the final credible ORFs"
        ),
    )
    internal.add_argument(
        "--orf-psite-selection",
        choices=("credible", "all"),
        default="credible",
        help="Fallback selection when no explicit ORF IDs are supplied",
    )
    internal.add_argument(
        "--orf-psite-id",
        action="append",
        default=[],
        help=(
            "Select candidate_key, gorf_id, or torf_id; repeat or provide comma-separated "
            "values. Supplying this option enables the codon stage"
        ),
    )
    internal.add_argument(
        "--orf-psite-list",
        help=(
            "Text/TSV ORF-ID list for batch export; supplying it enables the codon stage"
        ),
    )

    advanced = parser.add_argument_group("advanced pass-through arguments")
    for stage in STAGES:
        advanced.add_argument(
            f"--{stage}-extra",
            action="append",
            default=[],
            metavar="'ARGS'",
            help=f"Additional arguments appended to {stage}; repeat as needed",
        )
    return parser


def validate_cli(args: argparse.Namespace, parser: argparse.ArgumentParser) -> List[str]:
    if args.workers < 1:
        parser.error("--workers must be >= 1")
    if args.min_aa < 1 or args.min_intorf_aa < 1 or args.primary_min_intorf_aa < args.min_intorf_aa:
        parser.error("length thresholds must satisfy 1 <= min-intorf-aa <= primary-min-intorf-aa")
    if args.integration_min_replicates < 2:
        parser.error("--integration-min-replicates must be >= 2")
    if not 0 < args.integration_min_fraction <= 1:
        parser.error("--integration-min-fraction must be in (0, 1]")
    formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    if not formats or any(item not in {"png", "pdf"} for item in formats):
        parser.error("--formats must contain png and/or pdf")
    if args.dpi < 1:
        parser.error("--dpi must be >= 1")
    out_root = path_arg(args.out_root)
    if args.sample and out_root is not None and out_root.name.casefold() == args.sample.casefold():
        parser.error(
            "--out-root must be the shared analysis root, not the sample directory; "
            f"remove the trailing '/{args.sample}' because --sample adds sample-level paths"
        )
    if args.only_stage:
        args.from_stage = args.only_stage
        args.to_stage = args.only_stage
    if stage_index(args.from_stage) > stage_index(args.to_stage):
        parser.error("--from-stage must not come after --to-stage")
    selected = list(STAGES[stage_index(args.from_stage): stage_index(args.to_stage) + 1])
    codon_explicitly_requested = bool(
        args.export_orf_psites
        or args.orf_psite_id
        or args.orf_psite_list
        or args.orf_psite_selection != "credible"
        or args.only_stage == "codon"
        or args.from_stage == "codon"
        or args.to_stage == "codon"
        or args.codon_out_prefix
        or args.codon_extra
    )
    if "codon" in selected and not codon_explicitly_requested:
        selected.remove("codon")
    if args.batch:
        if "psite" not in selected:
            parser.error("--batch currently requires the psite stage to be selected")
        incompatible = {
            "--bam": args.bam,
            "--offsets": args.offsets,
            "--ribotish-quality": args.ribotish_quality,
            "--length-offsets": args.length_offsets,
            "--keep-lengths": args.keep_lengths,
            "--psite-bam": args.psite_bam,
            "--dm-background": args.dm_background,
            "--dm-results": args.dm_results,
            "--dm-run-manifest": args.dm_run_manifest,
            "--gene-context-results": args.gene_context_results,
            "--psite-out-prefix": args.psite_out_prefix,
            "--qc-out-prefix": args.qc_out_prefix,
            "--dm-out-prefix": args.dm_out_prefix,
            "--abundance-out-prefix": args.abundance_out_prefix,
            "--context-out-prefix": args.context_out_prefix,
            "--codon-out-prefix": args.codon_out_prefix,
            "--visualization-out-prefix": args.visualization_out_prefix,
        }
        supplied = [name for name, value in incompatible.items() if value is not None]
        if supplied:
            parser.error(
                "--batch obtains replicate-specific inputs and outputs from its TSV; "
                f"do not also provide: {', '.join(supplied)}"
            )
    elif "psite" in selected:
        if args.length_offsets is not None:
            incompatible = args.offsets or args.ribotish_quality or args.keep_lengths
            if incompatible:
                parser.error(
                    "--length-offsets is a complete traditional input mode and cannot be "
                    "combined with --offsets, --ribotish-quality or --keep-lengths"
                )
        elif not args.offsets:
            parser.error(
                "psite requires either traditional --length-offsets or Ribo-TISH --offsets"
            )
        elif not args.ribotish_quality and not args.keep_lengths:
            parser.error(
                "automatic Ribo-TISH mode requires --ribotish-quality; the legacy "
                "--offsets + --keep-lengths combination remains available for compatibility"
            )
    if args.dm_mode == "both" and (
        ("dm" in selected and args.dm_out_prefix)
        or ("abundance" in selected and args.abundance_out_prefix)
        or ("context" in selected and args.context_out_prefix)
        or ("codon" in selected and args.codon_out_prefix)
        or ("visualize" in selected and args.visualization_out_prefix)
    ):
        parser.error(
            "--dm-mode both cannot be combined with --dm-out-prefix, "
            "--abundance-out-prefix, --context-out-prefix, --codon-out-prefix or "
            "--visualization-out-prefix because the two "
            "modes require separate paths"
        )
    return selected


def build_paths(args: argparse.Namespace, selected: Sequence[str]) -> Paths:
    project = Path(__file__).resolve().parent
    source = resolve_stage_source(project)
    out_root = path_arg(args.out_root)
    assert out_root is not None
    gtf = path_arg(args.gtf)
    reference = args.reference_name or infer_reference_name(gtf)

    orf_prefix = path_arg(args.orfscan_out_prefix) or out_root / "01_orfscan" / reference
    psite_prefix = path_arg(args.psite_out_prefix) or out_root / "02_psite" / args.sample / args.sample
    qc_prefix = path_arg(args.qc_out_prefix) or out_root / "03_metagene_qc" / args.sample / args.sample
    dm_prefix = path_arg(args.dm_out_prefix) or out_root / "04_intorf_dm" / args.sample / args.dm_mode / args.sample
    abundance_prefix = path_arg(args.abundance_out_prefix) or dm_prefix
    context_prefix = path_arg(args.context_out_prefix) or dm_prefix
    codon_prefix = path_arg(args.codon_out_prefix) or dm_prefix
    visualization_prefix = (
        path_arg(args.visualization_out_prefix)
        or out_root / "05_visualization" / args.sample / args.dm_mode / args.sample
    )

    torf_generated = Path(f"{orf_prefix}.torf.tsv")
    psite_generated = Path(f"{psite_prefix}.psite.bam")
    background_generated = Path(f"{qc_prefix}.dm_background.tsv")
    dm_results_generated = Path(f"{dm_prefix}.intorf_dm_results.tsv")
    dm_run_manifest_generated = Path(f"{dm_prefix}.run_manifest.json")
    gene_context_generated = Path(f"{context_prefix}.gene_cds_context.tsv")

    torf = torf_generated if "orfscan" in selected else (path_arg(args.torf) or torf_generated)
    psite_bam = psite_generated if "psite" in selected else (path_arg(args.psite_bam) or psite_generated)
    dm_background = (
        background_generated if "qc" in selected else (path_arg(args.dm_background) or background_generated)
    )
    dm_results = dm_results_generated if "dm" in selected else (path_arg(args.dm_results) or dm_results_generated)
    dm_run_manifest = (
        dm_run_manifest_generated
        if "dm" in selected
        else (path_arg(args.dm_run_manifest) or dm_run_manifest_generated)
    )
    gene_context_results = (
        gene_context_generated
        if "context" in selected
        else (path_arg(args.gene_context_results) or gene_context_generated)
    )

    return Paths(
        project=project,
        source=source,
        out_root=out_root,
        reference_name=reference,
        sample=args.sample,
        orf_prefix=orf_prefix,
        per_chrom_dir=orf_prefix.parent / "per_chrom",
        psite_prefix=psite_prefix,
        qc_prefix=qc_prefix,
        dm_prefix=dm_prefix,
        abundance_prefix=abundance_prefix,
        context_prefix=context_prefix,
        codon_prefix=codon_prefix,
        visualization_prefix=visualization_prefix,
        torf=torf,
        psite_bam=psite_bam,
        dm_background=dm_background,
        dm_results=dm_results,
        dm_run_manifest=dm_run_manifest,
        gene_context_results=gene_context_results,
    )


def append_bool(command: List[str], flag: str, enabled: bool, negative: Optional[str] = None) -> None:
    if enabled:
        command.append(flag)
    elif negative:
        command.append(negative)


def append_extra(command: List[str], values: Iterable[str]) -> None:
    for value in values:
        command.extend(shlex.split(value))


def bam_index_candidates(path: Path) -> List[Path]:
    candidates = [Path(f"{path}.bai")]
    if path.suffix == ".bam":
        candidates.append(path.with_suffix(".bai"))
    return candidates


def existing_bam_index(path: Path) -> Path:
    for candidate in bam_index_candidates(path):
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    return bam_index_candidates(path)[0]


def build_stage_specs(args: argparse.Namespace, paths: Paths) -> Dict[str, StageSpec]:
    python = args.python
    gtf = path_arg(args.gtf)
    fa = path_arg(args.fa)
    bam = path_arg(args.bam)
    offsets = path_arg(args.offsets)
    ribotish_quality = path_arg(args.ribotish_quality)

    orf_cmd = [
        python, str(paths.source / "orf_scan_transcriptome.py"),
        "--gtf", str(gtf or ""), "--fa", str(fa or ""),
        "--out-prefix", str(paths.orf_prefix), "--min-aa", str(args.min_aa),
        "--start", *args.start_codons, "--workers", str(args.workers),
        "--mp-chunksize", str(args.orfscan_mp_chunksize),
    ]
    append_bool(orf_cmd, "--primary-only", args.primary_only)
    if args.near_cognate:
        orf_cmd.append("--near-cognate")
    if args.by_chrom:
        orf_cmd.extend(["--by-chrom", "--perchrom-outdir", str(paths.per_chrom_dir)])
    append_extra(orf_cmd, args.orfscan_extra)

    psite_cmd = [
        python, str(paths.source / "psite-caller.py"),
        "--bam", str(bam or ""),
        "--out-prefix", str(paths.psite_prefix), "--min-mapq", str(args.min_mapq),
        "--workers", str(args.workers), "--bgzf-threads", str(args.bgzf_threads), "--merge",
    ]
    if args.length_offsets is not None:
        psite_cmd.extend(
            ["--length-offsets", *[f"{length}:{offset}" for length, offset in sorted(args.length_offsets.items())]]
        )
    else:
        psite_cmd.extend(
            [
                "--offsets", str(offsets or ""),
                "--min-frame0-proportion", str(args.min_frame0_proportion),
                "--length-selection-policy", args.length_selection_policy,
            ]
        )
    if ribotish_quality is not None:
        psite_cmd.extend(["--ribotish-quality", str(ribotish_quality)])
    if args.keep_lengths is not None:
        psite_cmd.extend(["--keep-lengths", *[str(value) for value in args.keep_lengths]])
    append_bool(psite_cmd, "--require-unique", args.require_unique, "--no-require-unique")
    if not args.emit_bed:
        psite_cmd.append("--no-bed")
    if not args.emit_bedgraph:
        psite_cmd.append("--no-bedgraph")
    append_extra(psite_cmd, args.psite_extra)

    qc_cmd = [
        python, str(paths.source / "ribo_metagene_qc.py"),
        "--psite-bam", str(paths.psite_bam), "--torf", str(paths.torf),
        "--out-prefix", str(paths.qc_prefix), "--pi-method", args.pi_method,
        "--min-a0", str(args.min_a0), "--workers", str(args.workers),
    ]
    append_extra(qc_cmd, args.qc_extra)

    dm_cmd = [
        python, str(paths.source / "ribo_intorf_dm_caller.py"),
        "--psite-bam", str(paths.psite_bam), "--torf", str(paths.torf),
        "--dm-background", str(paths.dm_background), "--out-prefix", str(paths.dm_prefix),
        "--pi-method", args.pi_method, "--min-a0", str(args.min_a0),
        "--min-intorf-aa", str(args.min_intorf_aa),
        "--primary-min-intorf-aa", str(args.primary_min_intorf_aa),
        "--min-core-codons", str(args.min_core_codons),
        "--min-active-core-codons", str(args.min_active_core_codons),
        "--min-core-reads", str(args.min_core_reads),
        "--min-credible-active-core-codons", str(args.min_credible_active_core_codons),
        "--min-credible-active-core-frac", str(args.min_credible_active_core_frac),
        "--min-credible-target-residual-frac", str(args.min_credible_target_residual_frac),
        "--bootstrap-gate-p", "0.20", "--bootstrap-reps", str(args.bootstrap_reps),
        "--bootstrap-engine", "adaptive_importance", "--importance-mode", args.dm_mode,
        "--importance-reps", str(args.importance_reps),
        "--importance-iid-exceedance-threshold", "10",
        "--importance-etas", "0,0.05,0.1,0.2,0.4,0.7,1.0",
        "--importance-min-tail-ess", "30", "--importance-max-relative-se", "0.25",
        "--disable-block-bootstrap", "--seed", str(args.seed),
        "--fdr-threshold", str(args.fdr_threshold), "--lambda-min", str(args.lambda_min),
        "--template-separation-min", "0.05", "--lambda-grid-size", "101",
        "--lambda-abs-diff-max", str(args.lambda_abs_diff_max),
        "--lambda-rel-diff-max", str(args.lambda_rel_diff_max),
        "--lambda-rel-eps", "0.01",
        "--distance-to-segment-max", str(args.distance_to_segment_max),
        "--candidate-dedup", "gorf", "--require-gorf-id", "--member-preview-n", "20",
        "--workers", str(args.workers), "--mp-chunksize", str(args.dm_mp_chunksize),
    ]
    append_extra(dm_cmd, args.dm_extra)

    abundance_cmd = [
        python, str(paths.source / "quantify_intorf_abundance.py"),
        "--dm-results", str(paths.dm_results),
        "--psite-bam", str(paths.psite_bam),
        "--out-prefix", str(paths.abundance_prefix),
        "--sample", paths.sample,
    ]
    append_extra(abundance_cmd, args.abundance_extra)

    context_cmd = [
        python, str(paths.source / "annotate_gene_cds_context.py"),
        "--torf", str(paths.torf),
        "--input-tsv", str(paths.dm_results),
        "--out-prefix", str(paths.context_prefix),
    ]
    append_extra(context_cmd, args.context_extra)

    codon_cmd = [
        python, str(paths.source / "export_orf_psites.py"),
        "--input-results", str(paths.gene_context_results),
        "--torf", str(paths.torf),
        "--psite-bam", str(paths.psite_bam),
        "--dm-run-manifest", str(paths.dm_run_manifest),
        "--out-prefix", str(paths.codon_prefix),
        "--sample", paths.sample,
        "--selection", args.orf_psite_selection,
    ]
    for value in args.orf_psite_id:
        codon_cmd.extend(["--orf-id", value])
    if args.orf_psite_list:
        codon_cmd.extend(["--orf-list", args.orf_psite_list])
    if (
        args.orf_psite_selection == "credible"
        and not args.orf_psite_id
        and not args.orf_psite_list
    ):
        codon_cmd.append("--gene-level-pure-intorf-only")
    append_extra(codon_cmd, args.codon_extra)

    viz_cmd = [
        python, str(paths.source / "visualize_intorf_dm_results.py"),
        "--dm-results", str(paths.gene_context_results),
        "--gene-level-pure-intorf-only",
        "--out-prefix", str(paths.visualization_prefix),
        "--fdr-threshold", str(args.fdr_threshold), "--lambda-min", str(args.lambda_min),
        "--lambda-abs-diff-max", str(args.lambda_abs_diff_max),
        "--lambda-rel-diff-max", str(args.lambda_rel_diff_max),
        "--distance-to-segment-max", str(args.distance_to_segment_max),
        "--min-active-core-codons", str(args.min_credible_active_core_codons),
        "--min-active-core-frac", str(args.min_credible_active_core_frac),
        "--min-target-residual-frac", str(args.min_credible_target_residual_frac),
        "--formats", args.formats, "--dpi", str(args.dpi), "--seed", str(args.seed),
    ]
    if args.include_diagnostics:
        viz_cmd.append("--include-diagnostics")
    append_extra(viz_cmd, args.visualize_extra)

    visualization_stems = [
        "model_geometry_ternary", "model_geometry_ternary_by_frame",
        "gate_waterfall", "breadth_gate_plane", "effect_significance",
        "gate_combinations", "length_stratification", "depth_effects",
    ]
    if args.include_diagnostics:
        visualization_stems.extend(
            ["branch_aligned_geometry", "lambda_gain_vs_drop", "gate_margin_heatmap"]
        )
    visualization_formats = [
        item.strip().lower() for item in args.formats.split(",") if item.strip()
    ]
    visualization_outputs = [
        Path(f"{paths.visualization_prefix}.{stem}.{extension}")
        for stem in visualization_stems
        for extension in visualization_formats
    ] + [
        Path(f"{paths.visualization_prefix}.plot_data.tsv"),
        Path(f"{paths.visualization_prefix}.gate_waterfall.tsv"),
        Path(f"{paths.visualization_prefix}.gate_combinations.tsv"),
        Path(f"{paths.visualization_prefix}.length_stratification.tsv"),
        Path(f"{paths.visualization_prefix}.gate_margin_candidates.tsv"),
        Path(f"{paths.visualization_prefix}.manifest.tsv"),
    ]

    return {
        "orfscan": StageSpec(
            "orfscan", orf_cmd,
            [item for item in (gtf, fa, Path(f"{fa}.fai") if fa else None) if item is not None],
            [paths.torf, Path(f"{paths.orf_prefix}.gorf.tsv")],
            paths.orf_prefix.parent / "run.log",
        ),
        "psite": StageSpec(
            "psite", psite_cmd,
            [
                item for item in (
                    bam,
                    existing_bam_index(bam) if bam else None,
                    offsets,
                    ribotish_quality,
                ) if item is not None
            ],
            [
                paths.psite_bam,
                Path(f"{paths.psite_bam}.bai"),
                Path(f"{paths.psite_prefix}.offset_selection.tsv"),
            ],
            paths.psite_prefix.parent / "run.log",
        ),
        "qc": StageSpec(
            "qc", qc_cmd,
            [paths.psite_bam, existing_bam_index(paths.psite_bam), paths.torf],
            [
                paths.dm_background,
                Path(f"{paths.qc_prefix}.qc_summary.tsv"),
                Path(f"{paths.qc_prefix}.codon_frame_ternary.tsv"),
                Path(f"{paths.qc_prefix}.codon_frame_ternary.png"),
                Path(f"{paths.qc_prefix}.codon_frame_ternary.pdf"),
            ],
            paths.qc_prefix.parent / "run.log",
        ),
        "dm": StageSpec(
            "dm", dm_cmd,
            [paths.psite_bam, existing_bam_index(paths.psite_bam), paths.torf, paths.dm_background],
            [
                paths.dm_results,
                Path(f"{paths.dm_prefix}.summary.tsv"),
                paths.dm_run_manifest,
            ],
            paths.dm_prefix.parent / "run.log",
        ),
        "abundance": StageSpec(
            "abundance", abundance_cmd,
            [paths.dm_results, paths.psite_bam, existing_bam_index(paths.psite_bam)],
            [
                Path(f"{paths.abundance_prefix}.intorf_abundance.tsv"),
                Path(f"{paths.abundance_prefix}.intorf_abundance_summary.tsv"),
                Path(f"{paths.abundance_prefix}.intorf_abundance_manifest.json"),
            ],
            paths.abundance_prefix.parent / "abundance.log",
        ),
        "context": StageSpec(
            "context", context_cmd,
            [paths.dm_results, paths.torf],
            [
                paths.gene_context_results,
                Path(f"{paths.context_prefix}.gene_cds_context_only.tsv"),
                Path(f"{paths.context_prefix}.gene_cds_context_summary.tsv"),
                Path(f"{paths.context_prefix}.gene_cds_context_manifest.json"),
            ],
            paths.context_prefix.parent / "gene_context.log",
        ),
        "codon": StageSpec(
            "codon", codon_cmd,
            [
                paths.gene_context_results,
                paths.torf,
                paths.psite_bam,
                existing_bam_index(paths.psite_bam),
                paths.dm_run_manifest,
            ],
            [
                Path(f"{paths.codon_prefix}.orf_psite_codons.tsv.gz"),
                Path(f"{paths.codon_prefix}.orf_psite_ternary.tsv"),
                Path(f"{paths.codon_prefix}.orf_psite_summary.tsv"),
                Path(f"{paths.codon_prefix}.orf_psite_manifest.json"),
            ],
            paths.codon_prefix.parent / "codon_ternary.log",
        ),
        "visualize": StageSpec(
            "visualize", viz_cmd,
            [paths.gene_context_results],
            visualization_outputs,
            paths.visualization_prefix.parent / "run.log",
        ),
    }


def nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def completion_file_identity(path: Path, force_hash: bool = False) -> Dict[str, object]:
    absolute = path.absolute()
    identity: Dict[str, object] = {"path": str(absolute), "exists": absolute.is_file()}
    if not identity["exists"]:
        return identity
    stat = absolute.stat()
    identity.update({"size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    hashable_suffixes = {".csv", ".json", ".py", ".tsv", ".txt", ".yml", ".yaml"}
    should_hash = force_hash or (
        stat.st_size <= IDENTITY_HASH_LIMIT_BYTES
        and absolute.suffix.casefold() in hashable_suffixes
    )
    identity["sha256"] = sha256_file(absolute) if should_hash else None
    identity["sha256_status"] = "computed" if should_hash else "not_computed"
    return identity


def validate_output_file(path: Path) -> List[str]:
    """Return content-validation errors for one expected stage output."""
    if not nonempty(path):
        return [f"missing/empty: {path}"]
    suffix = path.suffix.casefold()
    try:
        if path.name.endswith(".orf_psite_codons.tsv.gz"):
            import gzip
            with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
                header = next(csv.reader(handle, delimiter="\t"), [])
            required_columns = {
                "candidate_key", "gorf_id", "torf_id", "core_codon_index0",
                "phase0_reads", "phase1_reads", "phase2_reads", "codon_reads",
                "phase0_fraction", "phase1_fraction", "phase2_fraction",
                "phase0_percent", "phase1_percent", "phase2_percent",
            }
            missing_columns = sorted(required_columns - set(header))
            if missing_columns:
                return [
                    f"missing required column(s) {', '.join(missing_columns)}: {path}"
                ]
        elif suffix in {".tsv", ".csv", ".txt"}:
            delimiter = "\t" if suffix in {".tsv", ".txt"} else ","
            with path.open("r", encoding="utf-8-sig", newline="") as handle:
                header = next(csv.reader(handle, delimiter=delimiter), [])
            if not header or not any(str(value).strip() for value in header):
                return [f"missing tabular header: {path}"]
            required_columns: set[str] = set()
            name = path.name
            if name.endswith(".torf.tsv"):
                required_columns = {"torf_id", "gorf_id"}
            elif name.endswith(".gorf.tsv"):
                required_columns = {"gorf_id"}
            elif name.endswith(".offset_selection.tsv"):
                required_columns = {"length", "offset", "selected", "reason"}
            elif name.endswith(".codon_frame_ternary.tsv"):
                required_columns = {
                    "grid_P0_percent", "grid_P1_percent", "grid_P2_percent",
                    "codon_count", "total_codons", "sample_percentage",
                }
            elif name.endswith(".intorf_dm_results.tsv"):
                required_columns = {
                    "gorf_id", "overlap_type", "classification", "p_final", "q_BH",
                    "lambda_hat", "core_reads",
                }
            elif name.endswith(".intorf_abundance.tsv"):
                required_columns = {
                    "gorf_id", "effective_core_nt", "usable_library_psites",
                    "model_expected_intorf_core_reads", "intorf_pFPKM",
                    "abundance_status",
                }
            elif name.endswith(".gene_cds_context.tsv"):
                required_columns = {
                    "gorf_id", "gene_id", "gene_level_orf_class",
                    "gene_level_pure_intorf_eligible", "gene_cds_nterm_match",
                    "gene_cds_nterm_coordinate_prefix_complete_codons",
                }
            elif name.endswith(".gene_cds_context_only.tsv"):
                required_columns = {
                    "gorf_id", "gene_id", "gene_level_orf_class",
                    "gene_level_pure_intorf_eligible", "gene_cds_nterm_match",
                }
            elif name.endswith(".orf_psite_ternary.tsv"):
                required_columns = {
                    "candidate_key", "gorf_id", "torf_id", "phase_ratio_P0",
                    "phase_ratio_P1", "phase_ratio_P2", "exact_P0_percent",
                    "exact_P1_percent", "exact_P2_percent", "codon_count",
                    "total_plotted_codons", "candidate_percentage",
                }
            elif name.endswith(".orf_psite_summary.tsv"):
                required_columns = {
                    "candidate_key", "gorf_id", "torf_id", "n_core_codons",
                    "n_active_core_codons", "n_plotted_codons", "core_reads",
                    "codon_equal_mean_P0_percent", "codon_equal_mean_P1_percent",
                    "codon_equal_mean_P2_percent", "count_reconciliation_status",
                }
            elif name.endswith(".pfpkm_matrix.tsv"):
                required_columns = {
                    "candidate_key", "gorf_id", "overlap_type", "gene_name",
                    "n_primary_credible",
                }
            elif name.endswith(".lambda_matrix.tsv"):
                required_columns = {
                    "candidate_key", "gorf_id", "overlap_type", "gene_name",
                    "n_primary_credible", "max_core_reads",
                }
            elif name.endswith(".pfpkm_correlations.tsv"):
                required_columns = {
                    "replicate_left", "replicate_right", "shared_candidates",
                    "pearson_log2p1_pFPKM", "spearman_log2p1_pFPKM",
                }
            elif name.endswith(
                (
                    ".summary.tsv", ".qc_summary.tsv", ".intorf_abundance_summary.tsv",
                    ".gene_cds_context_summary.tsv",
                )
            ):
                required_columns = {"metric", "value"}
            elif name.endswith(".manifest.tsv"):
                required_columns = {"key", "value"}
            missing_columns = sorted(required_columns - set(header))
            if missing_columns:
                return [
                    f"missing required column(s) {', '.join(missing_columns)}: {path}"
                ]
        elif suffix == ".json":
            with path.open("r", encoding="utf-8") as handle:
                value = json.load(handle)
            if not isinstance(value, (dict, list)):
                return [f"JSON root is not an object/array: {path}"]
        elif suffix == ".bam":
            samtools = shutil.which("samtools")
            if samtools is None:
                return [f"cannot validate BAM without samtools on PATH: {path}"]
            checked = subprocess.run(
                [samtools, "quickcheck", "-v", str(path)],
                capture_output=True,
                text=True,
                timeout=120,
                check=False,
            )
            if checked.returncode != 0:
                detail = (checked.stderr or checked.stdout).strip()
                return [f"samtools quickcheck failed: {path}: {detail}"]
        elif suffix == ".bai":
            with path.open("rb") as handle:
                if handle.read(4) != b"BAI\x01":
                    return [f"invalid BAI magic: {path}"]
        elif suffix == ".png":
            with path.open("rb") as handle:
                if handle.read(8) != b"\x89PNG\r\n\x1a\n":
                    return [f"invalid PNG signature: {path}"]
        elif suffix == ".pdf":
            with path.open("rb") as handle:
                if handle.read(5) != b"%PDF-":
                    return [f"invalid PDF signature: {path}"]
    except (OSError, UnicodeError, csv.Error, json.JSONDecodeError, subprocess.SubprocessError) as exc:
        return [f"unreadable/invalid output: {path}: {exc}"]
    return []


def validate_stage_outputs(spec: StageSpec) -> List[str]:
    errors = [error for path in spec.expected_outputs for error in validate_output_file(path)]
    if spec.name == "psite" and len(spec.expected_outputs) >= 2 and not errors:
        samtools = shutil.which("samtools")
        if samtools is None:
            errors.append("cannot validate BAM index without samtools on PATH")
        else:
            try:
                indexed = subprocess.run(
                    [samtools, "idxstats", str(spec.expected_outputs[0])],
                    capture_output=True,
                    text=True,
                    timeout=120,
                    check=False,
                )
            except subprocess.SubprocessError as exc:
                errors.append(f"samtools idxstats failed: {exc}")
            else:
                if indexed.returncode != 0 or not indexed.stdout.strip():
                    detail = (indexed.stderr or indexed.stdout).strip()
                    errors.append(f"BAM index validation failed: {detail}")
    if spec.name == "dm" and len(spec.expected_outputs) >= 3 and not errors:
        manifest_path = spec.expected_outputs[2]
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"invalid DM run manifest: {manifest_path}: {exc}")
        else:
            required = {
                "program", "version", "output_schema_version", "statistical_engine_id",
                "parameters", "inputs", "outputs", "result_rows",
            }
            missing = sorted(required - set(manifest))
            if missing:
                errors.append(
                    f"DM run manifest missing field(s) {', '.join(missing)}: {manifest_path}"
                )
    if spec.name == "abundance" and len(spec.expected_outputs) >= 3 and not errors:
        manifest_path = spec.expected_outputs[2]
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"invalid abundance manifest: {manifest_path}: {exc}")
        else:
            required = {
                "program", "version", "output_schema_version", "formula",
                "semantics", "inputs", "outputs", "usable_library_psites",
                "result_rows",
            }
            missing = sorted(required - set(manifest))
            if missing:
                errors.append(
                    "abundance manifest missing field(s) "
                    f"{', '.join(missing)}: {manifest_path}"
                )
    if spec.name == "codon" and len(spec.expected_outputs) >= 4 and not errors:
        manifest_path = spec.expected_outputs[3]
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"invalid credible-intORF codon manifest: {manifest_path}: {exc}")
        else:
            required = {
                "program", "version", "schema_version", "selection", "parameters",
                "semantics", "inputs", "outputs", "counts",
            }
            missing = sorted(required - set(manifest))
            if missing:
                errors.append(
                    "credible-intORF codon manifest missing field(s) "
                    f"{', '.join(missing)}: {manifest_path}"
                )
    return errors


def stage_completion_path(spec: StageSpec) -> Path:
    return Path(f"{spec.expected_outputs[0]}.complete.json")


def stage_contract(spec: StageSpec) -> Dict[str, object]:
    program_path = Path(spec.command[1]) if len(spec.command) > 1 else Path()
    return {
        "completion_schema_version": COMPLETION_SCHEMA_VERSION,
        "program": PROGRAM,
        "program_version": __version__,
        "stage": spec.name,
        "command": [str(value) for value in spec.command],
        "stage_program": completion_file_identity(program_path, force_hash=True),
        "required_inputs": [completion_file_identity(path) for path in spec.required_inputs],
        "expected_outputs": [completion_file_identity(path) for path in spec.expected_outputs],
    }


def write_stage_completion(spec: StageSpec) -> Path:
    errors = validate_stage_outputs(spec)
    if errors:
        raise RuntimeError("\n".join(errors))
    destination = stage_completion_path(spec)
    destination.parent.mkdir(parents=True, exist_ok=True)
    record = {
        **stage_contract(spec),
        "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "content_validation": "passed",
    }
    temporary = Path(f"{destination}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")
    temporary.replace(destination)
    return destination


def stage_completion_state(spec: StageSpec) -> tuple[bool, str]:
    completion = stage_completion_path(spec)
    if not completion.is_file():
        if any(nonempty(path) for path in spec.expected_outputs):
            return False, "outputs exist but completion record is missing"
        return False, "outputs are missing"
    try:
        recorded = json.loads(completion.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return False, f"completion record is unreadable: {exc}"
    errors = validate_stage_outputs(spec)
    if errors:
        return False, errors[0]
    current = stage_contract(spec)
    for field in (
        "completion_schema_version", "program", "program_version", "stage", "command",
        "stage_program", "required_inputs", "expected_outputs",
    ):
        if recorded.get(field) != current.get(field):
            return False, f"completion contract changed: {field}"
    if recorded.get("content_validation") != "passed":
        return False, "completion record lacks successful content validation"
    return True, "validated completion record"


def command_text(command: Sequence[str]) -> str:
    return shlex.join(str(value) for value in command)


def print_status(selected: Sequence[str], specs: Dict[str, StageSpec]) -> None:
    print(f"{PROGRAM} status")
    for stage in selected:
        spec = specs[stage]
        complete, reason = stage_completion_state(spec)
        state = "COMPLETE" if complete else ("STALE" if any(nonempty(path) for path in spec.expected_outputs) else "MISSING")
        print(f"  {stage:10s} {state} ({reason})")
        print(f"    completion: {stage_completion_path(spec)}")
        for path in spec.expected_outputs:
            state = f"{path.stat().st_size:,} bytes" if nonempty(path) else "missing/empty"
            print(f"    {path} [{state}]")


def execution_plan(
    args: argparse.Namespace,
    selected: Sequence[str],
) -> tuple[
    Dict[str, argparse.Namespace],
    Dict[str, Paths],
    Dict[str, Dict[str, StageSpec]],
    List[tuple[str, str, StageSpec]],
]:
    """Build one shared upstream path followed by one or two formal DM paths."""
    modes = ("accurate", "fast") if args.dm_mode == "both" else (args.dm_mode,)
    mode_args: Dict[str, argparse.Namespace] = {}
    mode_paths: Dict[str, Paths] = {}
    mode_specs: Dict[str, Dict[str, StageSpec]] = {}
    for mode in modes:
        current = copy.copy(args)
        current.dm_mode = mode
        paths = build_paths(current, selected)
        mode_args[mode] = current
        mode_paths[mode] = paths
        mode_specs[mode] = build_stage_specs(current, paths)

    primary_mode = modes[0]
    plan: List[tuple[str, str, StageSpec]] = []
    for stage in selected:
        stage_modes = modes if stage in MODE_SPECIFIC_STAGES else (primary_mode,)
        for mode in stage_modes:
            plan.append((stage, mode, mode_specs[mode][stage]))
    return mode_args, mode_paths, mode_specs, plan


def print_execution_status(plan: Sequence[tuple[str, str, StageSpec]], two_modes: bool) -> None:
    print(f"{PROGRAM} status")
    for stage, mode, spec in plan:
        label = f"{stage}:{mode}" if two_modes and stage in MODE_SPECIFIC_STAGES else stage
        complete, reason = stage_completion_state(spec)
        state = "COMPLETE" if complete else ("STALE" if any(nonempty(path) for path in spec.expected_outputs) else "MISSING")
        print(f"  {label:17s} {state} ({reason})")
        print(f"    completion: {stage_completion_path(spec)}")
        for path in spec.expected_outputs:
            state = f"{path.stat().st_size:,} bytes" if nonempty(path) else "missing/empty"
            print(f"    {path} [{state}]")


def require_inputs(stage: str, spec: StageSpec, args: argparse.Namespace) -> None:
    missing = [path for path in spec.required_inputs if not nonempty(path)]
    if missing:
        joined = "\n".join(f"  - {path}" for path in missing)
        raise RuntimeError(f"{stage}: missing or empty required input(s):\n{joined}")
    if stage == "orfscan" and (not args.gtf or not args.fa):
        raise RuntimeError("orfscan: --gtf and --fa are required")
    if stage == "psite" and not args.bam:
        raise RuntimeError("psite: --bam is required")
    if stage == "psite" and shutil.which("samtools") is None:
        raise RuntimeError(
            "psite: 'samtools' is not available on PATH; activate the complete "
            "runtime environment before starting so this is detected before "
            "per-chromosome P-site generation"
        )


def stream_command(command: Sequence[str], log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8", newline="") as log:
        log.write(f"COMMAND: {command_text(command)}\n\n")
        log.flush()
        process = subprocess.Popen(
            list(command),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
        return process.wait()


def write_manifest(paths: Paths, record: Dict[str, object]) -> Path:
    run_dir = paths.out_root / "pipeline_runs" / paths.sample
    run_dir.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    destination = run_dir / f"{stamp}.pipeline_run.json"
    temporary = Path(f"{destination}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    temporary.replace(destination)
    return destination


def read_metric_tsv(path: Path) -> Dict[str, str]:
    if not nonempty(path):
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"metric", "value"}.issubset(reader.fieldnames):
            return {}
        return {
            str(row.get("metric", "")): str(row.get("value", ""))
            for row in reader
            if row.get("metric")
        }


def selected_offset_text(path: Path) -> str:
    if not nonempty(path):
        return ""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        selected = []
        for row in reader:
            if str(row.get("selected", "")).strip().casefold() not in {"1", "true", "yes"}:
                continue
            try:
                selected.append((int(row["length"]), int(row["offset"])))
            except (KeyError, TypeError, ValueError):
                continue
    return ",".join(f"{length}:{offset}" for length, offset in sorted(selected))


def bam_indexed_alignment_count(path: Path) -> object:
    if not nonempty(path) or shutil.which("samtools") is None:
        return ""
    try:
        result = subprocess.run(
            ["samtools", "idxstats", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
        return sum(
            int(fields[2])
            for line in result.stdout.splitlines()
            if len(fields := line.split("\t")) >= 3
        )
    except (OSError, subprocess.CalledProcessError, ValueError):
        return ""


def prepare_batch_args(
    base: argparse.Namespace, replicate: BatchReplicate, index: int
) -> argparse.Namespace:
    current = copy.copy(base)
    current.sample = replicate.replicate_id
    current.bam = str(replicate.bam)
    current.length_offsets = (
        dict(replicate.length_offsets) if replicate.length_offsets is not None else None
    )
    current.offsets = str(replicate.ribotish_para) if replicate.ribotish_para else None
    current.ribotish_quality = (
        str(replicate.ribotish_quality) if replicate.ribotish_quality else None
    )
    current.keep_lengths = None
    current.batch_replicate_index = index
    return current


def build_batch_summary_rows(outcomes: Sequence[PipelineOutcome]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    alignment_counts: Dict[Path, object] = {}
    for outcome in outcomes:
        if not outcome.mode_paths:
            rows.append(
                {
                    "replicate_id": outcome.sample,
                    "input_mode": outcome.input_mode,
                    "dm_mode": "",
                    "status": outcome.status,
                    "failed_stage": outcome.failed_stage,
                    "error": outcome.error,
                }
            )
            continue
        for mode, paths in outcome.mode_paths.items():
            offset_audit = Path(f"{paths.psite_prefix}.offset_selection.tsv")
            qc_summary = Path(f"{paths.qc_prefix}.qc_summary.tsv")
            dm_summary = Path(f"{paths.dm_prefix}.summary.tsv")
            qc = read_metric_tsv(qc_summary)
            dm = read_metric_tsv(dm_summary)
            if paths.psite_bam not in alignment_counts:
                alignment_counts[paths.psite_bam] = bam_indexed_alignment_count(paths.psite_bam)
            rows.append(
                {
                    "replicate_id": outcome.sample,
                    "input_mode": outcome.input_mode,
                    "dm_mode": mode,
                    "status": outcome.status,
                    "failed_stage": outcome.failed_stage,
                    "error": outcome.error,
                    "selected_length_offsets": selected_offset_text(offset_audit),
                    "psite_bam": str(paths.psite_bam),
                    "psite_bam_bytes": paths.psite_bam.stat().st_size if nonempty(paths.psite_bam) else "",
                    "psite_alignments": alignment_counts[paths.psite_bam],
                    "qc_template_status": qc.get("template_status", ""),
                    "qc_frame0_prop": qc.get("frame0_prop", ""),
                    "qc_A0": qc.get("A0", ""),
                    "candidate_intorf_altframe_total": dm.get("candidate_intorf_altframe_total", ""),
                    "primary_credible_calls": dm.get(
                        "classification::credible_extra_ORF_like_signal",
                        "0" if dm else "",
                    ),
                    "dm_results": str(paths.dm_results),
                    "gene_context_results": str(paths.gene_context_results),
                    "pipeline_manifest": str(outcome.manifest or ""),
                }
            )
    return rows


def write_batch_records(
    args: argparse.Namespace,
    outcomes: Sequence[PipelineOutcome],
    integrations: Sequence[IntegrationOutcome],
    started_at: str,
) -> tuple[Path, Path]:
    out_root = path_arg(args.out_root)
    assert out_root is not None
    destination = out_root / "batch_runs"
    destination.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    stem = Path(args.batch).stem
    summary_path = destination / f"{stem}.{stamp}.batch_summary.tsv"
    manifest_path = destination / f"{stem}.{stamp}.batch_run.json"
    rows = build_batch_summary_rows(outcomes)
    fields = [
        "replicate_id", "input_mode", "dm_mode", "status", "failed_stage", "error",
        "selected_length_offsets", "psite_bam", "psite_bam_bytes",
        "psite_alignments",
        "qc_template_status", "qc_frame0_prop", "qc_A0",
        "candidate_intorf_altframe_total", "primary_credible_calls",
        "dm_results", "gene_context_results", "pipeline_manifest",
    ]
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    failed = [outcome for outcome in outcomes if outcome.status == "failed"]
    record = {
        "program": PROGRAM,
        "version": __version__,
        "run_type": "replicate_batch",
        "batch_tsv": str(Path(args.batch).resolve()),
        "started_at": started_at,
        "finished_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "partial_failed" if failed else "completed",
        "replicate_count": len(outcomes),
        "failed_replicate_count": len(failed),
        "command": [sys.executable, *sys.argv],
        "summary_tsv": str(summary_path),
        "replicates": rows,
        "integrations": [
            {
                "dm_mode": item.dm_mode,
                "status": item.status,
                "error": item.error,
                "command": item.command,
                "command_shell": command_text(item.command),
                "log": str(item.log),
                "outputs": [str(path) for path in item.outputs],
                "context_status": item.context_status,
                "context_error": item.context_error,
                "context_command": item.context_command or [],
                "context_command_shell": (
                    command_text(item.context_command) if item.context_command else ""
                ),
                "context_log": str(item.context_log or ""),
                "context_outputs": [str(path) for path in (item.context_outputs or [])],
                "visualization_status": item.visualization_status,
                "visualization_error": item.visualization_error,
                "visualization_command": item.visualization_command or [],
                "visualization_command_shell": (
                    command_text(item.visualization_command)
                    if item.visualization_command else ""
                ),
                "visualization_log": str(item.visualization_log or ""),
                "visualization_outputs": [
                    str(path) for path in (item.visualization_outputs or [])
                ],
            }
            for item in integrations
        ],
    }
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return summary_path, manifest_path


def run_replicate_integrations(
    args: argparse.Namespace,
    outcomes: Sequence[PipelineOutcome],
) -> List[IntegrationOutcome]:
    out_root = path_arg(args.out_root)
    assert out_root is not None
    source = resolve_stage_source(Path(__file__).resolve().parent)
    script = source / "integrate_replicates.py"
    context_script = source / "annotate_gene_cds_context.py"
    visualization_script = source / "visualize_replicate_integration.py"
    if not nonempty(script):
        return [
            IntegrationOutcome(
                dm_mode="all",
                status="failed",
                error=f"missing integration program: {script}",
                command=[],
                log=out_root / "06_replicate_integration" / "integration.log",
                outputs=[],
            )
        ]
    expected = [outcome.sample for outcome in outcomes]
    modes = ["accurate", "fast"] if args.dm_mode == "both" else [args.dm_mode]
    records: List[IntegrationOutcome] = []
    for mode in modes:
        prefix = out_root / "06_replicate_integration" / mode / Path(args.batch).stem
        sample_metadata_path = Path(f"{prefix}.sample_metadata.tsv")
        outputs = [
            Path(f"{prefix}.replicate_long.tsv"),
            Path(f"{prefix}.consensus.tsv"),
            Path(f"{prefix}.summary.tsv"),
            Path(f"{prefix}.manifest.json"),
        ]
        log = prefix.parent / "run.log"
        command = [
            args.python,
            str(script),
            "--out-prefix", str(prefix),
            "--min-replicates", str(args.integration_min_replicates),
            "--min-fraction", str(args.integration_min_fraction),
            "--fdr-threshold", str(args.fdr_threshold),
        ]
        for replicate in expected:
            command.extend(["--expected-replicate", replicate])
        for outcome in outcomes:
            paths = outcome.mode_paths.get(mode)
            dm_is_current = (
                outcome.status == "completed"
                or outcome.failed_stage in {"abundance", "context", "codon", "visualize"}
            )
            if paths is not None and dm_is_current and nonempty(paths.dm_results):
                command.extend(["--result", f"{outcome.sample}={paths.dm_results}"])
                run_manifest = Path(f"{paths.dm_prefix}.run_manifest.json")
                if nonempty(run_manifest):
                    command.extend(["--run-manifest", f"{outcome.sample}={run_manifest}"])
        if "--result" not in command:
            records.append(
                IntegrationOutcome(mode, "skipped_no_results", "", command, log, outputs)
            )
            continue
        print(f"\n===== replicate integration:{mode} =====")
        print(command_text(command))
        code = stream_command(command, log)
        output_errors = [
            error for path in outputs for error in validate_output_file(path)
        ]
        if code != 0:
            records.append(
                IntegrationOutcome(
                    mode, "failed", f"integration exited with status {code}", command, log, outputs
                )
            )
        elif output_errors:
            records.append(
                IntegrationOutcome(
                    mode,
                    "failed_output_validation",
                    "invalid output(s): " + "; ".join(output_errors),
                    command,
                    log,
                    outputs,
                )
            )
        else:
            context_prefix = Path(f"{prefix}.consensus")
            context_outputs = [
                Path(f"{context_prefix}.gene_cds_context.tsv"),
                Path(f"{context_prefix}.gene_cds_context_only.tsv"),
                Path(f"{context_prefix}.gene_cds_context_summary.tsv"),
                Path(f"{context_prefix}.gene_cds_context_manifest.json"),
            ]
            context_log = prefix.parent / "gene_context.log"
            torf_path = next(
                (
                    paths.torf
                    for outcome in outcomes
                    for paths in [outcome.mode_paths.get(mode)]
                    if paths is not None and nonempty(paths.torf)
                ),
                None,
            )
            context_command = [
                args.python,
                str(context_script),
                "--torf", str(torf_path or ""),
                "--input-tsv", str(outputs[1]),
                "--out-prefix", str(context_prefix),
            ]
            if not nonempty(context_script):
                context_status = "failed"
                context_error = f"missing gene-context program: {context_script}"
            elif torf_path is None:
                context_status = "failed"
                context_error = "no complete tORF table available for replicate consensus"
            else:
                print(f"\n===== replicate consensus gene context:{mode} =====")
                print(command_text(context_command))
                context_code = stream_command(context_command, context_log)
                context_validation_errors = [
                    error for path in context_outputs for error in validate_output_file(path)
                ]
                if context_code != 0:
                    context_status = "failed"
                    context_error = f"gene-context annotation exited with status {context_code}"
                elif context_validation_errors:
                    context_status = "failed_output_validation"
                    context_error = "invalid output(s): " + "; ".join(context_validation_errors)
                else:
                    context_status = "completed"
                    context_error = ""
            metadata_rows = []
            abundance_inputs: Dict[str, Path] = {}
            outcomes_by_sample = {outcome.sample: outcome for outcome in outcomes}
            for replicate in expected:
                outcome = outcomes_by_sample.get(replicate)
                paths = outcome.mode_paths.get(mode) if outcome is not None else None
                qc = (
                    read_metric_tsv(Path(f"{paths.qc_prefix}.qc_summary.tsv"))
                    if paths is not None else {}
                )
                metadata_rows.append({
                    "replicate_id": replicate,
                    "selected_length_offsets": (
                        selected_offset_text(
                            Path(f"{paths.psite_prefix}.offset_selection.tsv")
                        )
                        if paths is not None else ""
                    ),
                    "psite_alignments": (
                        bam_indexed_alignment_count(paths.psite_bam)
                        if paths is not None and nonempty(paths.psite_bam) else ""
                    ),
                    "frame0_prop": qc.get("frame0_prop", ""),
                    "A0": qc.get("A0", ""),
                })
                abundance_is_current = bool(
                    outcome is not None
                    and (
                        outcome.status == "completed"
                        or outcome.failed_stage in {"context", "codon", "visualize"}
                    )
                )
                abundance_path = (
                    Path(f"{paths.abundance_prefix}.intorf_abundance.tsv")
                    if paths is not None else None
                )
                if (
                    abundance_is_current
                    and abundance_path is not None
                    and nonempty(abundance_path)
                ):
                    abundance_inputs[replicate] = abundance_path
            sample_metadata_path.parent.mkdir(parents=True, exist_ok=True)
            with sample_metadata_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "replicate_id", "selected_length_offsets", "psite_alignments",
                        "frame0_prop", "A0",
                    ],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(metadata_rows)
            formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
            visualization_outputs = [
                Path(f"{prefix}.{figure}.{extension}")
                for figure in (
                    "call_counts", "support_combinations", "lambda_concordance",
                    "lambda_heatmap", "primary_reproducibility",
                )
                for extension in formats
            ] + [
                sample_metadata_path,
                Path(f"{prefix}.plot_summary.tsv"),
                Path(f"{prefix}.plot_manifest.json"),
                Path(f"{prefix}.lambda_matrix.tsv"),
            ]
            if abundance_inputs:
                visualization_outputs.extend([
                    *[
                        Path(f"{prefix}.{figure}.{extension}")
                        for figure in (
                            "pfpkm_correlation", "pfpkm_correlation_matrix",
                            "pfpkm_heatmap",
                        )
                        for extension in formats
                    ],
                    Path(f"{prefix}.pfpkm_matrix.tsv"),
                    Path(f"{prefix}.pfpkm_correlations.tsv"),
                ])
            visualization_log = prefix.parent / "visualization.log"
            visualization_command = [
                args.python,
                str(visualization_script),
                "--replicate-long", str(outputs[0]),
                "--consensus", str(context_outputs[0]),
                "--gene-level-pure-intorf-only",
                "--sample-metadata", str(sample_metadata_path),
                "--out-prefix", str(prefix),
                "--fdr-threshold", str(args.fdr_threshold),
                "--formats", args.formats,
                "--dpi", str(args.dpi),
            ]
            for replicate, abundance_path in abundance_inputs.items():
                visualization_command.extend([
                    "--abundance", f"{replicate}={abundance_path}"
                ])
            if not nonempty(visualization_script):
                record = IntegrationOutcome(
                    mode, "completed", "", command, log, outputs,
                    "failed", f"missing visualization program: {visualization_script}",
                    visualization_command, visualization_log, visualization_outputs,
                )
                record.context_status = context_status
                record.context_error = context_error
                record.context_command = context_command
                record.context_log = context_log
                record.context_outputs = context_outputs
                records.append(record)
                continue
            print(f"\n===== replicate integration visualization:{mode} =====")
            print(command_text(visualization_command))
            visualization_code = stream_command(visualization_command, visualization_log)
            visualization_errors = [
                error
                for path in visualization_outputs
                for error in validate_output_file(path)
            ]
            if visualization_code != 0:
                visualization_status = "failed"
                visualization_error = f"visualization exited with status {visualization_code}"
            elif visualization_errors:
                visualization_status = "failed_output_validation"
                visualization_error = "invalid output(s): " + "; ".join(visualization_errors)
            else:
                visualization_status = "completed"
                visualization_error = ""
            record = IntegrationOutcome(
                mode, "completed", "", command, log, outputs,
                visualization_status, visualization_error,
                visualization_command, visualization_log, visualization_outputs,
            )
            record.context_status = context_status
            record.context_error = context_error
            record.context_command = context_command
            record.context_log = context_log
            record.context_outputs = context_outputs
            records.append(record)
    return records


def run_single(
    args: argparse.Namespace,
    selected: Sequence[str],
    parser: argparse.ArgumentParser,
) -> PipelineOutcome:
    mode_args, mode_paths, mode_specs, plan = execution_plan(args, selected)
    primary_mode = "accurate" if args.dm_mode == "both" else args.dm_mode
    paths = mode_paths[primary_mode]
    if args.length_offsets is not None:
        input_mode = "explicit"
    elif args.ribotish_quality:
        input_mode = "ribotish"
    elif args.offsets:
        input_mode = "legacy_para_keep_lengths"
    else:
        input_mode = "not_applicable"

    for script in (
        "orf_scan_transcriptome.py", "psite-caller.py", "ribo_metagene_qc.py",
        "ribo_intorf_dm_caller.py", "quantify_intorf_abundance.py",
        "annotate_gene_cds_context.py",
        "export_orf_psites.py", "export_intorf_codon_ternary.py",
        "visualize_intorf_dm_results.py",
    ):
        if not nonempty(paths.source / script):
            parser.error(f"missing stage program: {paths.source / script}")

    if args.status:
        print_execution_status(plan, args.dm_mode == "both")
        return PipelineOutcome(paths.sample, input_mode, "status", "", "", None, mode_paths)

    record: Dict[str, object] = {
        "program": PROGRAM,
        "version": __version__,
        "started_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "sample": paths.sample,
        "reference_name": paths.reference_name,
        "from_stage": selected[0],
        "to_stage": selected[-1],
        "dm_mode": args.dm_mode,
        "dm_modes_executed": list(mode_paths),
        "dry_run": bool(args.dry_run),
        "force": bool(args.force),
        "command": [sys.executable, *sys.argv],
        "parameters": {
            key: value
            for key, value in sorted(vars(args).items())
        },
        "stages": [],
    }
    stage_records: List[Dict[str, object]] = record["stages"]  # type: ignore[assignment]
    failed_stage = ""

    try:
        for number, (stage, mode, spec) in enumerate(plan, start=1):
            failed_stage = stage
            complete, completion_reason = stage_completion_state(spec)
            label = (
                f"{stage}:{mode}"
                if args.dm_mode == "both" and stage in MODE_SPECIFIC_STAGES
                else stage
            )
            item: Dict[str, object] = {
                "stage": stage,
                "mode": mode if stage in MODE_SPECIFIC_STAGES else None,
                "title": STAGE_TITLES[stage],
                "command": spec.command,
                "command_shell": command_text(spec.command),
                "log": str(spec.log_path),
                "expected_outputs": [str(path) for path in spec.expected_outputs],
                "completion_record": str(stage_completion_path(spec)),
                "completion_state": completion_reason,
            }
            stage_records.append(item)
            print(f"\n===== [{number}/{len(plan)}] {label}: {STAGE_TITLES[stage]} =====")

            force_stage = bool(args.force) and not (
                stage == "orfscan" and getattr(args, "batch_replicate_index", 0) > 0
            )
            if complete and not force_stage:
                item["status"] = "skipped_complete"
                print("SKIP: validated completion record matches command, code, inputs and outputs")
                continue

            if args.dry_run:
                # External inputs are checked; outputs supplied by an earlier selected
                # stage are intentionally allowed to be absent in a dry run.
                earlier_outputs = {
                    output
                    for _, _, earlier_spec in plan[:number - 1]
                    for output in earlier_spec.expected_outputs
                }
                missing = [
                    path for path in spec.required_inputs
                    if not nonempty(path) and path not in earlier_outputs
                ]
                if missing:
                    raise RuntimeError(
                        f"{stage}: missing external input(s) for dry run:\n"
                        + "\n".join(f"  - {path}" for path in missing)
                    )
                item["status"] = "dry_run"
                print(command_text(spec.command))
                continue

            completion_path = stage_completion_path(spec)
            if completion_path.exists():
                completion_path.unlink()
            require_inputs(stage, spec, mode_args[mode])
            print(command_text(spec.command))
            code = stream_command(spec.command, spec.log_path)
            item["exit_code"] = code
            if code != 0:
                item["status"] = "failed"
                raise RuntimeError(f"{stage}: command exited with status {code}; log: {spec.log_path}")
            output_errors = validate_stage_outputs(spec)
            if output_errors:
                item["status"] = "failed_output_validation"
                raise RuntimeError(
                    f"{stage}: command finished but output validation failed:\n"
                    + "\n".join(f"  - {error}" for error in output_errors)
                )
            completion_path = write_stage_completion(spec)
            item["completion_record"] = str(completion_path)
            item["completion_state"] = "validated completion record"
            item["status"] = "completed"

        record["status"] = "dry_run" if args.dry_run else "completed"
        record["finished_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
    except Exception as exc:
        record["status"] = "failed"
        record["error"] = str(exc)
        record["finished_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
        manifest = None
        if not args.dry_run:
            manifest = write_manifest(paths, record)
            sys.stderr.write(f"Pipeline record: {manifest}\n")
        sys.stderr.write(f"ERROR: {exc}\n")
        return PipelineOutcome(
            paths.sample, input_mode, "failed", str(exc), failed_stage, manifest, mode_paths
        )

    manifest = None
    if not args.dry_run:
        manifest = write_manifest(paths, record)
        print(f"\nPipeline record: {manifest}")
    print("DONE")
    return PipelineOutcome(
        paths.sample,
        input_mode,
        "dry_run" if args.dry_run else "completed",
        "",
        "",
        manifest,
        mode_paths,
    )


def run_batch(
    args: argparse.Namespace,
    selected: Sequence[str],
    parser: argparse.ArgumentParser,
) -> int:
    try:
        replicates = read_batch_replicates(args.batch)
    except (OSError, ValueError, argparse.ArgumentTypeError) as exc:
        parser.error(f"invalid --batch TSV: {exc}")

    modes = sorted({replicate.input_mode for replicate in replicates})
    if len(modes) > 1:
        print(
            "WARN: batch mixes explicit and Ribo-TISH input modes; "
            "the actual mode is retained in every audit and summary row"
        )
    started_at = dt.datetime.now(dt.timezone.utc).isoformat()
    outcomes: List[PipelineOutcome] = []
    for index, replicate in enumerate(replicates):
        print(
            f"\n######## replicate [{index + 1}/{len(replicates)}] "
            f"{replicate.replicate_id} ({replicate.input_mode}) ########"
        )
        current = prepare_batch_args(args, replicate, index)
        try:
            outcome = run_single(current, selected, parser)
        except Exception as exc:
            # Keep a malformed replicate from preventing diagnostics for later rows.
            outcome = PipelineOutcome(
                replicate.replicate_id,
                replicate.input_mode,
                "failed",
                str(exc),
                "setup",
                None,
                {},
            )
            sys.stderr.write(f"ERROR [{replicate.replicate_id}]: {exc}\n")
        outcomes.append(outcome)

    print("\nBatch replicate status:")
    for outcome in outcomes:
        suffix = f" ({outcome.failed_stage}: {outcome.error})" if outcome.error else ""
        print(f"  {outcome.sample}\t{outcome.input_mode}\t{outcome.status}{suffix}")

    integrations: List[IntegrationOutcome] = []
    if not args.dry_run and not args.status and "dm" in selected:
        integrations = run_replicate_integrations(args, outcomes)
        print("\nReplicate integration status:")
        for item in integrations:
            suffix = f" ({item.error})" if item.error else ""
            context = f"\tcontext={item.context_status}" if item.context_status else ""
            if item.context_error:
                context += f" ({item.context_error})"
            visualization = (
                f"\tvisualization={item.visualization_status}"
                if item.visualization_status else ""
            )
            if item.visualization_error:
                visualization += f" ({item.visualization_error})"
            print(f"  {item.dm_mode}\t{item.status}{suffix}{context}{visualization}")

    if not args.dry_run and not args.status:
        summary_path, manifest_path = write_batch_records(
            args, outcomes, integrations, started_at
        )
        print(f"Batch summary: {summary_path}")
        print(f"Batch record: {manifest_path}")
    failed = any(outcome.status == "failed" for outcome in outcomes) or any(
        item.status.startswith("failed")
        or item.context_status.startswith("failed")
        or item.visualization_status.startswith("failed")
        for item in integrations
    )
    return 1 if failed else 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    selected = validate_cli(args, parser)
    if args.batch:
        return run_batch(args, selected, parser)
    outcome = run_single(args, selected, parser)
    return 1 if outcome.status == "failed" else 0


if __name__ == "__main__":
    raise SystemExit(main())
