#!/usr/bin/env python3
"""Export auditable codon-level P-site data for selected ORF calls.

The tool is intentionally optional and is designed for internal inspection of
one ORF, a supplied ORF list, or the default set of primary credible calls. The
output contains observed P-site composition in each candidate-relative codon.
Because an intORF overlaps its host CDS, these counts describe the observed
mixture and must not be interpreted as hard read-level assignment to the
intORF component.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import importlib.util
import json
import math
import os
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

try:
    from . import ribo_intorf_dm_caller as dm
except ImportError:  # Direct script execution.
    dm_path = Path(__file__).with_name("ribo_intorf_dm_caller.py")
    dm_spec = importlib.util.spec_from_file_location("miao_orf_codon_ternary_dm", dm_path)
    if dm_spec is None or dm_spec.loader is None:
        raise ImportError(f"cannot load DM caller from {dm_path}")
    dm = importlib.util.module_from_spec(dm_spec)
    import sys
    sys.modules[dm_spec.name] = dm
    dm_spec.loader.exec_module(dm)


PROGRAM = "miao-orf-orf-psite-export"
VERSION = "1.0.0"
SCHEMA_VERSION = "1.2"
SUPPORTED_ID_COLUMNS = ("orf_id", "candidate_key", "gorf_id", "torf_id")

LONG_COLUMNS = [
    "sample", "candidate_key", "gorf_id", "torf_id", "gene_id", "gene_name",
    "transcript_id", "overlap_type", "classification", "q_BH", "lambda_hat",
    "chrom", "strand", "orf_codon_index0", "core_codon_index0",
    "phase0_genomic_pos1", "phase1_genomic_pos1", "phase2_genomic_pos1",
    "phase0_reads", "phase1_reads", "phase2_reads", "codon_reads",
    "phase0_fraction", "phase1_fraction", "phase2_fraction",
    "phase0_percent", "phase1_percent", "phase2_percent",
]

TERNARY_COLUMNS = [
    "sample", "candidate_key", "gorf_id", "torf_id", "gene_id", "gene_name",
    "transcript_id", "overlap_type", "classification", "q_BH", "lambda_hat",
    "phase_ratio_P0", "phase_ratio_P1", "phase_ratio_P2",
    "exact_P0_fraction", "exact_P1_fraction", "exact_P2_fraction",
    "exact_P0_percent", "exact_P1_percent", "exact_P2_percent",
    "codon_count", "total_plotted_codons", "candidate_percentage",
]

SUMMARY_COLUMNS = [
    "sample", "candidate_key", "gorf_id", "torf_id", "gene_id", "gene_name",
    "transcript_id", "overlap_type", "classification", "q_BH", "lambda_hat",
    "gene_level_pure_intorf_eligible", "n_core_codons", "n_active_core_codons",
    "n_plotted_codons", "core_reads", "plotted_codon_reads",
    "codon_equal_mean_P0_percent", "codon_equal_mean_P1_percent",
    "codon_equal_mean_P2_percent", "P0_vertex_percent", "P1_vertex_percent",
    "P2_vertex_percent", "count_reconciliation_status",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export exact per-codon P-site data for selected ORFs."
    )
    parser.add_argument("--input-results", required=True, help="DM or gene-context result TSV")
    parser.add_argument("--torf", required=True, help="tORF TSV, directory, glob, or prefix")
    parser.add_argument("--psite-bam", required=True, help="Sample-matched indexed 1-nt P-site BAM")
    parser.add_argument(
        "--dm-run-manifest",
        help=(
            "DM run manifest containing the exact core exclusions and optional read-length "
            "filter; inferred beside *.intorf_dm_results.tsv when omitted"
        ),
    )
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--sample", default="sample")
    parser.add_argument("--min-codon-reads", type=int, default=1)
    parser.add_argument(
        "--selection",
        choices=("credible", "all"),
        default="credible",
        help="Selection used when no explicit ORF IDs are supplied",
    )
    parser.add_argument(
        "--orf-id",
        action="append",
        default=[],
        help=(
            "Select candidate_key, gorf_id, or torf_id values; repeat this option "
            "or provide comma-separated values"
        ),
    )
    parser.add_argument(
        "--orf-list",
        help=(
            "Text/TSV file containing ORF IDs. Accepts one ID per line or columns named "
            "orf_id, candidate_key, gorf_id, and/or torf_id"
        ),
    )
    parser.add_argument(
        "--gene-level-pure-intorf-only",
        action="store_true",
        help="Also require gene_level_pure_intorf_eligible=1",
    )
    parser.add_argument(
        "--exclude-start-codons",
        type=int,
        help="Override the DM-manifest start-codon exclusion",
    )
    parser.add_argument(
        "--exclude-stop-codons",
        type=int,
        help="Override the DM-manifest stop-codon exclusion",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return parser


def truthy(value: object) -> bool:
    return str(value).strip().casefold() in {"1", "true", "yes", "y"}


def finite_int(value: object) -> Optional[int]:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return int(number) if math.isfinite(number) and number.is_integer() else None


def finite_float(value: object) -> Optional[float]:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def format_value(value: object) -> object:
    if isinstance(value, float):
        return "" if not math.isfinite(value) else f"{value:.10g}"
    return value


def candidate_key(row: Mapping[str, str]) -> str:
    gorf_id = str(row.get("gorf_id", "")).strip()
    overlap = str(row.get("overlap_type", "")).strip()
    return f"{gorf_id}|{overlap}" if gorf_id else str(row.get("torf_id", "")).strip()


def infer_manifest_path(input_results: str) -> Path:
    path = Path(input_results)
    suffix = ".intorf_dm_results.tsv"
    if path.name.endswith(suffix):
        return path.with_name(path.name[:-len(suffix)] + ".run_manifest.json")
    name = path.name
    for suffix in (".gene_cds_context.tsv", ".gene_cds_context_only.tsv"):
        if name.endswith(suffix):
            return path.with_name(name[:-len(suffix)] + ".run_manifest.json")
    raise SystemExit(
        "ERROR: cannot infer --dm-run-manifest from this result filename; supply it explicitly"
    )


def load_run_contract(args: argparse.Namespace) -> Tuple[Path, int, int, Optional[int], str]:
    manifest_path = Path(args.dm_run_manifest) if args.dm_run_manifest else infer_manifest_path(args.input_results)
    if not manifest_path.is_file():
        raise SystemExit(f"ERROR: DM run manifest not found: {manifest_path}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise SystemExit(f"ERROR: invalid DM run manifest {manifest_path}: {exc}") from exc
    parameters = manifest.get("parameters")
    if not isinstance(parameters, dict):
        raise SystemExit(f"ERROR: DM run manifest has no parameters object: {manifest_path}")

    start = args.exclude_start_codons
    stop = args.exclude_stop_codons
    if start is None:
        start = finite_int(parameters.get("exclude_start_codons"))
    if stop is None:
        stop = finite_int(parameters.get("exclude_stop_codons"))
    if start is None or stop is None or start < 0 or stop < 0:
        raise SystemExit(
            "ERROR: valid start/stop codon exclusions must come from the DM manifest or explicit overrides"
        )
    rl_raw = parameters.get("rl")
    rl = finite_int(rl_raw) if rl_raw not in (None, "", "pooled") else None
    rl_source = str(parameters.get("rl_source", "auto") or "auto")
    return manifest_path, start, stop, rl, rl_source


def reduced_phase_ratio(counts: Sequence[int]) -> Tuple[int, int, int]:
    values = tuple(int(value) for value in counts)
    if len(values) != 3 or any(value < 0 for value in values) or sum(values) <= 0:
        raise ValueError("exact ternary ratio requires three non-negative counts with a positive sum")
    divisor = math.gcd(math.gcd(values[0], values[1]), values[2])
    return values[0] // divisor, values[1] // divisor, values[2] // divisor


def split_requested_ids(values: Iterable[str]) -> set[str]:
    requested: set[str] = set()
    for value in values:
        requested.update(token.strip() for token in str(value).split(",") if token.strip())
    return requested


def load_orf_list(path: str) -> set[str]:
    list_path = Path(path)
    if not list_path.is_file():
        raise SystemExit(f"ERROR: ORF list not found: {list_path}")
    with list_path.open("r", encoding="utf-8-sig", newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise SystemExit(f"ERROR: ORF list contains no IDs: {list_path}")

    first_fields = [value.strip() for value in lines[0].rstrip("\r\n").split("\t")]
    supported_headers = [value for value in first_fields if value in SUPPORTED_ID_COLUMNS]
    if supported_headers:
        reader = csv.DictReader(lines, delimiter="\t")
        requested = {
            str(row.get(column, "")).strip()
            for row in reader
            for column in supported_headers
            if str(row.get(column, "")).strip()
        }
    else:
        requested = {
            line.rstrip("\r\n").split("\t", 1)[0].strip()
            for line in lines
            if line.rstrip("\r\n").split("\t", 1)[0].strip()
        }
    if not requested:
        raise SystemExit(f"ERROR: ORF list contains no usable IDs: {list_path}")
    return requested


def row_identifiers(row: Mapping[str, str]) -> set[str]:
    return {
        value
        for value in (
            candidate_key(row),
            str(row.get("gorf_id", "")).strip(),
            str(row.get("torf_id", "")).strip(),
        )
        if value
    }


def load_selected_rows(
    path: str,
    pure_only: bool,
    selection: str = "credible",
    requested_ids: Optional[set[str]] = None,
) -> Tuple[List[Dict[str, str]], List[str]]:
    requested = set(requested_ids or set())
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        required = {"torf_id", "gorf_id", "overlap_type"}
        if not requested and selection == "credible":
            required.add("primary_credible_call")
        missing = sorted(required - set(fields))
        if missing:
            raise SystemExit(f"ERROR: input results missing columns: {', '.join(missing)}")
        if pure_only and not requested and "gene_level_pure_intorf_eligible" not in fields:
            raise SystemExit(
                "ERROR: --gene-level-pure-intorf-only requires gene_level_pure_intorf_eligible"
            )
        all_rows = [dict(row) for row in reader]

    if requested:
        rows = [row for row in all_rows if row_identifiers(row) & requested]
        matched = {value for row in rows for value in row_identifiers(row) if value in requested}
        missing_requested = sorted(requested - matched)
        if missing_requested:
            preview = ", ".join(missing_requested[:10])
            raise SystemExit(
                f"ERROR: {len(missing_requested)} requested ORF ID(s) were not found: {preview}"
            )
    elif selection == "credible":
        rows = [
            row for row in all_rows
            if truthy(row.get("primary_credible_call", ""))
            and (not pure_only or truthy(row.get("gene_level_pure_intorf_eligible", "")))
        ]
    else:
        rows = [
            row for row in all_rows
            if not pure_only or truthy(row.get("gene_level_pure_intorf_eligible", ""))
        ]

    if not rows:
        raise SystemExit("ERROR: ORF selection produced no rows")
    seen: set[str] = set()
    unique: List[Dict[str, str]] = []
    for row in rows:
        key = candidate_key(row)
        if key in seen:
            raise SystemExit(f"ERROR: duplicate selected candidate key in input results: {key}")
        seen.add(key)
        unique.append(row)
    return unique, fields


def load_torf_lookup(spec: str, wanted_ids: set[str]) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    paths = dm.list_input_torfs(spec)
    lookup: Dict[str, Dict[str, str]] = {}
    required = {"torf_id", "chrom", "strand", "blockSizes", "genomic_block_starts1"}
    for path in paths:
        with open(path, "r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = sorted(required - set(reader.fieldnames or []))
            if missing:
                raise SystemExit(f"ERROR: tORF table {path} missing columns: {', '.join(missing)}")
            for row in reader:
                torf_id = str(row.get("torf_id", ""))
                if torf_id not in wanted_ids:
                    continue
                current = lookup.get(torf_id)
                if current is not None and current != row:
                    raise SystemExit(f"ERROR: conflicting tORF rows for {torf_id}")
                lookup[torf_id] = dict(row)
    missing_ids = sorted(wanted_ids - set(lookup))
    if missing_ids:
        preview = ", ".join(missing_ids[:5])
        raise SystemExit(f"ERROR: {len(missing_ids)} selected torf_id values were not found: {preview}")
    return lookup, paths


def build_orf(raw: Dict[str, str]) -> dm.ORFRow:
    blocks = dm.parse_blocks_01(raw.get("blockSizes", ""), raw.get("genomic_block_starts1", ""))
    if not blocks:
        raise SystemExit(f"ERROR: invalid genomic blocks for torf_id={raw.get('torf_id', '')}")
    return dm.ORFRow(
        raw=raw,
        torf_id=raw.get("torf_id", ""),
        gene_id=raw.get("gene_id", ""),
        gene_name=raw.get("gene_name", ""),
        transcript_id=raw.get("transcript_id", ""),
        chrom=raw.get("chrom", ""),
        strand=raw.get("strand", "+") or "+",
        blocks=blocks,
        t_start=finite_int(raw.get("t_start")) or 0,
        t_end=finite_int(raw.get("t_end")) or 0,
        aa_len=finite_int(raw.get("aa_len")) or 0,
        peptide_len=finite_int(raw.get("peptide_len", raw.get("aa_len"))) or 0,
        cds_first_t=finite_int(raw.get("cds_first_t")),
        orf_biotype=raw.get("orf_biotype", "intORF_altframe"),
        is_annotated_cds=False,
        has_mane=truthy(raw.get("has_mane", "")),
        gorf_id=raw.get("gorf_id", ""),
    )


def genomic_position_for_tindex(row: dm.ORFRow, tindex: int) -> int:
    remaining = int(tindex)
    blocks = row.blocks if row.strand == "+" else list(reversed(row.blocks))
    for start, end in blocks:
        length = end - start
        if remaining < length:
            return start + remaining if row.strand == "+" else end - 1 - remaining
        remaining -= length
    raise IndexError(f"tindex outside ORF blocks: {tindex}")


def reconcile_counts(result: Mapping[str, str], matrix: object) -> List[str]:
    errors: List[str] = []
    phase_totals = matrix.sum(axis=0)
    total_reads = int(phase_totals.sum())
    active = int((matrix.sum(axis=1) > 0).sum())
    checks = {
        "core_reads": total_reads,
        "analyzed_core_reads": total_reads,
        "n_active_core_codons": active,
        "n_core_codons": int(len(matrix)),
        "observed_phase0_core_reads": int(phase_totals[0]),
        "observed_phase1_core_reads": int(phase_totals[1]),
        "observed_phase2_core_reads": int(phase_totals[2]),
    }
    for column, observed in checks.items():
        expected = finite_int(result.get(column, ""))
        if expected is not None and expected != observed:
            errors.append(f"{column}: result={expected}, recomputed={observed}")
    return errors


def metadata_values(result: Mapping[str, str]) -> Dict[str, object]:
    return {
        "sample": "",
        "candidate_key": candidate_key(result),
        "gorf_id": result.get("gorf_id", ""),
        "torf_id": result.get("torf_id", ""),
        "gene_id": result.get("gene_id", ""),
        "gene_name": result.get("gene_name", ""),
        "transcript_id": result.get("transcript_id", ""),
        "overlap_type": result.get("overlap_type", ""),
        "classification": result.get("classification", ""),
        "q_BH": result.get("q_BH", ""),
        "lambda_hat": result.get("lambda_hat", ""),
    }


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    temporary = Path(f"{path}.tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: format_value(row.get(column, "")) for column in columns})
    os.replace(temporary, path)


def write_gzip_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    temporary = Path(f"{path}.tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(temporary, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: format_value(row.get(column, "")) for column in columns})
    os.replace(temporary, path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(path: Path) -> Dict[str, object]:
    absolute = path.resolve()
    stat = absolute.stat()
    return {
        "path": str(absolute),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256_file(absolute) if stat.st_size <= 64 * 1024 * 1024 else None,
    }


def existing_bam_index(path: Path) -> Path:
    candidates = [Path(f"{path}.bai")]
    if path.suffix.casefold() == ".bam":
        candidates.append(path.with_suffix(".bai"))
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    raise SystemExit(f"ERROR: indexed P-site BAM is required; checked: {', '.join(map(str, candidates))}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.min_codon_reads < 1:
        parser.error("--min-codon-reads must be >= 1")
    requested_ids = split_requested_ids(args.orf_id)
    if args.orf_list:
        requested_ids.update(load_orf_list(args.orf_list))
    dm.require_dependencies()

    input_results = Path(args.input_results)
    psite_bam = Path(args.psite_bam)
    if not input_results.is_file():
        raise SystemExit(f"ERROR: input results not found: {input_results}")
    if not psite_bam.is_file():
        raise SystemExit(f"ERROR: P-site BAM not found: {psite_bam}")
    psite_bam_index = existing_bam_index(psite_bam)
    manifest_path, exclude_start, exclude_stop, rl, rl_source = load_run_contract(args)
    selected_rows, _ = load_selected_rows(
        str(input_results),
        args.gene_level_pure_intorf_only,
        args.selection,
        requested_ids,
    )
    torf_lookup, torf_paths = load_torf_lookup(
        args.torf, {row["torf_id"] for row in selected_rows}
    )

    long_rows: List[Dict[str, object]] = []
    ternary_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []
    with dm.pysam.AlignmentFile(str(psite_bam), "rb") as bam:
        for result in selected_rows:
            orf = build_orf(torf_lookup[result["torf_id"]])
            tcounts = dm.count_psites_by_tindex(bam, orf, rl, rl_source)
            matrix_all = dm.codon_count_matrix(tcounts, orf.orf_len_nt)
            core_start = min(exclude_start, len(matrix_all))
            core_end = max(core_start, len(matrix_all) - exclude_stop)
            matrix_core = matrix_all[core_start:core_end]
            mismatches = reconcile_counts(result, matrix_core)
            if mismatches:
                raise SystemExit(
                    f"ERROR: count reconciliation failed for {candidate_key(result)}: "
                    + "; ".join(mismatches)
                )

            base = metadata_values(result)
            base["sample"] = args.sample
            ternary_counter: Counter[Tuple[int, int, int]] = Counter()
            mean_sum = [0.0, 0.0, 0.0]
            plotted_reads = 0
            plotted_codons = 0
            for core_index, counts_array in enumerate(matrix_core):
                counts = [int(value) for value in counts_array]
                depth = sum(counts)
                if depth < args.min_codon_reads:
                    continue
                ratio_key = reduced_phase_ratio(counts)
                ternary_counter[ratio_key] += 1
                fractions = [value / depth for value in counts]
                for phase in range(3):
                    mean_sum[phase] += fractions[phase]
                plotted_reads += depth
                plotted_codons += 1
                orf_index = core_start + core_index
                positions = [
                    genomic_position_for_tindex(orf, 3 * orf_index + phase) + 1
                    for phase in range(3)
                ]
                row: Dict[str, object] = dict(base)
                row.update({
                    "chrom": orf.chrom,
                    "strand": orf.strand,
                    "orf_codon_index0": orf_index,
                    "core_codon_index0": core_index,
                    "phase0_genomic_pos1": positions[0],
                    "phase1_genomic_pos1": positions[1],
                    "phase2_genomic_pos1": positions[2],
                    "phase0_reads": counts[0],
                    "phase1_reads": counts[1],
                    "phase2_reads": counts[2],
                    "codon_reads": depth,
                    "phase0_fraction": fractions[0],
                    "phase1_fraction": fractions[1],
                    "phase2_fraction": fractions[2],
                    "phase0_percent": 100.0 * fractions[0],
                    "phase1_percent": 100.0 * fractions[1],
                    "phase2_percent": 100.0 * fractions[2],
                })
                long_rows.append(row)

            for key, count in sorted(ternary_counter.items(), key=lambda item: (-item[1], item[0])):
                ratio_total = sum(key)
                row = dict(base)
                row.update({
                    "phase_ratio_P0": key[0],
                    "phase_ratio_P1": key[1],
                    "phase_ratio_P2": key[2],
                    "exact_P0_fraction": key[0] / ratio_total,
                    "exact_P1_fraction": key[1] / ratio_total,
                    "exact_P2_fraction": key[2] / ratio_total,
                    "exact_P0_percent": 100.0 * key[0] / ratio_total,
                    "exact_P1_percent": 100.0 * key[1] / ratio_total,
                    "exact_P2_percent": 100.0 * key[2] / ratio_total,
                    "codon_count": count,
                    "total_plotted_codons": plotted_codons,
                    "candidate_percentage": 100.0 * count / plotted_codons if plotted_codons else 0.0,
                })
                ternary_rows.append(row)

            vertices = [
                ternary_counter.get((1, 0, 0), 0),
                ternary_counter.get((0, 1, 0), 0),
                ternary_counter.get((0, 0, 1), 0),
            ]
            summary = dict(base)
            summary.update({
                "gene_level_pure_intorf_eligible": result.get("gene_level_pure_intorf_eligible", ""),
                "n_core_codons": len(matrix_core),
                "n_active_core_codons": int((matrix_core.sum(axis=1) > 0).sum()),
                "n_plotted_codons": plotted_codons,
                "core_reads": int(matrix_core.sum()),
                "plotted_codon_reads": plotted_reads,
                "codon_equal_mean_P0_percent": 100.0 * mean_sum[0] / plotted_codons if plotted_codons else 0.0,
                "codon_equal_mean_P1_percent": 100.0 * mean_sum[1] / plotted_codons if plotted_codons else 0.0,
                "codon_equal_mean_P2_percent": 100.0 * mean_sum[2] / plotted_codons if plotted_codons else 0.0,
                "P0_vertex_percent": 100.0 * vertices[0] / plotted_codons if plotted_codons else 0.0,
                "P1_vertex_percent": 100.0 * vertices[1] / plotted_codons if plotted_codons else 0.0,
                "P2_vertex_percent": 100.0 * vertices[2] / plotted_codons if plotted_codons else 0.0,
                "count_reconciliation_status": "exact_match",
            })
            summary_rows.append(summary)

    prefix = Path(args.out_prefix)
    long_path = Path(f"{prefix}.orf_psite_codons.tsv.gz")
    ternary_path = Path(f"{prefix}.orf_psite_ternary.tsv")
    summary_path = Path(f"{prefix}.orf_psite_summary.tsv")
    manifest_output = Path(f"{prefix}.orf_psite_manifest.json")
    write_gzip_tsv(long_path, LONG_COLUMNS, long_rows)
    write_tsv(ternary_path, TERNARY_COLUMNS, ternary_rows)
    write_tsv(summary_path, SUMMARY_COLUMNS, summary_rows)
    manifest = {
        "program": PROGRAM,
        "version": VERSION,
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "sample": args.sample,
        "selection": {
            "mode": "explicit_ids" if requested_ids else args.selection,
            "requested_ids": sorted(requested_ids),
            "primary_credible_call": not requested_ids and args.selection == "credible",
            "gene_level_pure_intorf_only_requested": bool(
                args.gene_level_pure_intorf_only
            ),
            "gene_level_pure_intorf_only_effective": bool(
                args.gene_level_pure_intorf_only and not requested_ids
            ),
            "explicit_ids_override_default_filters": bool(requested_ids),
        },
        "parameters": {
            "min_codon_reads": args.min_codon_reads,
            "exclude_start_codons": exclude_start,
            "exclude_stop_codons": exclude_stop,
            "rl": rl,
            "rl_source": rl_source,
        },
        "semantics": {
            "phase_coordinates": "candidate-intORF-relative P0/P1/P2",
            "counts": "observed P-sites in the DM-analyzed candidate core",
            "source_assignment": (
                "host-CDS and intORF-like components remain unresolved at read level; "
                "rows describe their observed overlap mixture"
            ),
            "ternary_denominator": "codons with codon_reads >= min_codon_reads",
            "ternary_coordinates": "exact normalized per-codon P0/P1/P2 counts",
            "exact_point_aggregation": (
                "codons are grouped only when their reduced integer P0:P1:P2 ratios match; "
                "no binning or coordinate rounding"
            ),
        },
        "counts": {
            "selected_candidates": len(selected_rows),
            "candidates_exported": len(summary_rows),
            "codon_rows": len(long_rows),
            "ternary_rows": len(ternary_rows),
            "exactly_reconciled_candidates": len(summary_rows),
        },
        "inputs": {
            "input_results": file_identity(input_results),
            "psite_bam": file_identity(psite_bam),
            "psite_bam_index": file_identity(psite_bam_index),
            "dm_run_manifest": file_identity(manifest_path),
            "torf_files": [file_identity(Path(path)) for path in torf_paths],
        },
        "outputs": {
            "codon_phase": str(long_path.resolve()),
            "ternary_exact": str(ternary_path.resolve()),
            "candidate_summary": str(summary_path.resolve()),
        },
    }
    temporary = Path(f"{manifest_output}.tmp")
    manifest_output.parent.mkdir(parents=True, exist_ok=True)
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, manifest_output)
    print(f"Wrote {long_path}")
    print(f"Wrote {ternary_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {manifest_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
