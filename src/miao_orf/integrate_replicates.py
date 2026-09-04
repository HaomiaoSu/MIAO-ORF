#!/usr/bin/env python3
"""Integrate independent MIAO DM results across replicates.

This module does not pool reads, refit the DM model, or combine p-values.
It summarizes replicate-level evidence for the formal candidate identity
``gorf_id + overlap_type`` and emits auditable long, matrix and consensus tables.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


PROGRAM = "miao-orf-integrate-replicates"
VERSION = "1.0.0"
PRIMARY_CREDIBLE_CLASS = "credible_extra_ORF_like_signal"

COMPATIBILITY_FIELDS = (
    "program", "version", "output_schema_version", "statistical_engine_id",
    "benchmark_certification_id",
)
SEMANTIC_PARAMETER_KEYS = (
    "pi_method", "min_a0", "rl", "rl_source", "template_min_cds_nt",
    "template_min_density", "template_trim_start_nt", "template_trim_stop_nt",
    "template_min_codon_reads", "trim_fraction", "trim_reservoir_size",
    "min_intorf_aa", "primary_min_intorf_aa",
    "min_core_codons", "min_active_core_codons", "min_core_reads",
    "min_credible_active_core_codons", "min_credible_active_core_frac",
    "min_credible_target_residual_frac", "exclude_start_codons",
    "exclude_stop_codons", "include_annotated_stop_confounded",
    "bootstrap_gate_p", "bootstrap_reps",
    "bootstrap_engine", "importance_mode", "importance_reps",
    "importance_iid_exceedance_threshold", "importance_etas",
    "importance_pilot_reps", "importance_confirm_min_reps",
    "importance_allocation_safety_factor", "importance_use_for_p_final",
    "importance_min_tail_ess", "importance_max_relative_se",
    "seed", "disable_block_bootstrap", "fdr_threshold", "lambda_min",
    "template_separation_min", "lambda_grid_size", "lambda_abs_diff_max",
    "lambda_rel_diff_max", "lambda_rel_eps", "distance_to_segment_max",
    "lag1_threshold", "min_active_for_lag1", "block_sizes", "review_gate_window",
    "candidate_dedup", "require_gorf_id", "candidate_mane_only",
    "no_preload_psites", "no_coverage_prefilter", "member_preview_n",
)

IDENTITY_COLUMNS = (
    "gorf_id", "overlap_type", "gene_id", "gene_name", "chrom", "strand",
    "orf_biotype", "aa_len", "peptide_len",
)
EVIDENCE_COLUMNS = (
    "classification", "primary_credible_call", "q_BH", "q_BY", "p_final",
    "lambda_hat", "core_reads", "n_active_core_codons", "active_core_codon_frac",
    "mixture_geometry_consistent", "distance_to_mixture_segment",
)


def parse_assignment(value: str, label: str) -> Tuple[str, str]:
    if "=" not in value:
        raise ValueError(f"{label} must use REPLICATE_ID=VALUE")
    replicate, assigned = value.split("=", 1)
    replicate = replicate.strip()
    assigned = assigned.strip()
    if not replicate or not assigned:
        raise ValueError(f"{label} must use non-empty REPLICATE_ID=VALUE")
    return replicate, assigned


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(path: Path, hash_limit_bytes: int = 64 * 1024 * 1024) -> Dict[str, object]:
    absolute = path.absolute()
    stat = absolute.stat()
    compute_hash = stat.st_size <= hash_limit_bytes
    return {
        "path": str(absolute),
        "size_bytes": stat.st_size,
        "sha256": sha256_file(absolute) if compute_hash else None,
        "sha256_status": "computed" if compute_hash else f"skipped_larger_than_{hash_limit_bytes}_bytes",
    }


def default_run_manifest_path(result_path: Path) -> Path:
    suffix = ".intorf_dm_results.tsv"
    rendered = str(result_path)
    if not rendered.endswith(suffix):
        raise ValueError(
            f"cannot infer run manifest from nonstandard DM result name: {result_path}; "
            "supply --run-manifest REPLICATE_ID=PATH"
        )
    return Path(rendered[:-len(suffix)] + ".run_manifest.json")


def load_run_manifest(path: Path, replicate_id: str) -> Dict[str, object]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"run manifest is missing or empty for {replicate_id}: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"run manifest is unreadable for {replicate_id}: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"run manifest root must be an object for {replicate_id}: {path}")
    return value


def normalized_torf_identities(manifest: Mapping[str, object]) -> List[Tuple[object, ...]]:
    inputs = manifest.get("inputs")
    torf_files = inputs.get("torf_files") if isinstance(inputs, dict) else None
    if not isinstance(torf_files, list) or not torf_files:
        raise ValueError("run manifest does not contain inputs.torf_files")
    normalized = []
    for item in torf_files:
        if not isinstance(item, dict) or not item.get("exists"):
            raise ValueError("run manifest contains an invalid tORF input identity")
        digest = item.get("sha256")
        if digest:
            normalized.append(("sha256", digest, item.get("size_bytes")))
        else:
            normalized.append(
                ("path_stat", item.get("path"), item.get("size_bytes"), item.get("mtime_utc"))
            )
    return sorted(normalized, key=lambda item: tuple(str(value) for value in item))


def compatibility_signature(manifest: Mapping[str, object]) -> Dict[str, object]:
    parameters = manifest.get("parameters")
    source = manifest.get("source")
    caller = source.get("caller") if isinstance(source, dict) else None
    if not isinstance(parameters, dict):
        raise ValueError("run manifest does not contain a parameters object")
    if not isinstance(caller, dict) or not caller.get("sha256"):
        raise ValueError("run manifest does not contain source.caller.sha256")
    missing = [field for field in COMPATIBILITY_FIELDS if field not in manifest]
    if missing:
        raise ValueError(f"run manifest missing compatibility field(s): {', '.join(missing)}")
    missing_parameters = [key for key in SEMANTIC_PARAMETER_KEYS if key not in parameters]
    if missing_parameters:
        raise ValueError(
            "run manifest missing semantic parameter(s): " + ", ".join(missing_parameters)
        )
    schema = manifest.get("output_schema")
    if not isinstance(schema, list) or not schema:
        raise ValueError("run manifest does not contain a non-empty output_schema")
    return {
        **{field: manifest[field] for field in COMPATIBILITY_FIELDS},
        "caller_sha256": caller["sha256"],
        "semantic_parameters": {key: parameters[key] for key in SEMANTIC_PARAMETER_KEYS},
        "torf_inputs": normalized_torf_identities(manifest),
        "output_schema": schema,
    }


def validate_result_against_manifest(
    replicate_id: str, result_path: Path, manifest: Mapping[str, object]
) -> None:
    with result_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        actual_schema = list(reader.fieldnames or [])
        actual_rows = sum(1 for _row in reader)
    if manifest.get("output_schema") != actual_schema:
        raise ValueError(f"DM result schema no longer matches its run manifest for {replicate_id}")
    if safe_int(manifest.get("result_rows")) != actual_rows:
        raise ValueError(f"DM result row count no longer matches its run manifest for {replicate_id}")
    outputs = manifest.get("outputs")
    if not isinstance(outputs, list):
        raise ValueError(f"run manifest outputs are invalid for {replicate_id}")
    result_name = str(result_path.absolute())
    candidates = [
        item for item in outputs
        if isinstance(item, dict)
        and str(item.get("path", "")) == result_name
    ]
    if not candidates:
        candidates = [
            item for item in outputs
            if isinstance(item, dict)
            and str(item.get("path", "")).endswith(".intorf_dm_results.tsv")
        ]
    if len(candidates) != 1:
        raise ValueError(
            f"run manifest must identify exactly one DM result for {replicate_id}: {result_path}"
        )
    recorded = candidates[0]
    current_size = result_path.stat().st_size
    if recorded.get("size_bytes") != current_size:
        raise ValueError(f"DM result size no longer matches its run manifest for {replicate_id}")
    recorded_sha256 = recorded.get("sha256")
    if recorded_sha256 and recorded_sha256 != sha256_file(result_path):
        raise ValueError(f"DM result SHA-256 no longer matches its run manifest for {replicate_id}")


def validate_manifest_compatibility(
    expected_replicates: Sequence[str],
    result_paths: Mapping[str, Path],
    manifest_paths: Mapping[str, Path],
) -> Tuple[Dict[str, Dict[str, object]], Dict[str, object]]:
    manifests: Dict[str, Dict[str, object]] = {}
    signatures: Dict[str, Dict[str, object]] = {}
    for replicate in expected_replicates:
        if replicate not in result_paths:
            continue
        manifest_path = manifest_paths.get(replicate)
        if manifest_path is None:
            raise ValueError(f"no run manifest was supplied or inferred for {replicate}")
        manifest = load_run_manifest(manifest_path, replicate)
        validate_result_against_manifest(replicate, result_paths[replicate], manifest)
        try:
            signature = compatibility_signature(manifest)
        except ValueError as exc:
            raise ValueError(f"incompatible run manifest for {replicate}: {exc}") from exc
        manifests[replicate] = manifest
        signatures[replicate] = signature
    baseline_replicate = next(iter(signatures))
    baseline = signatures[baseline_replicate]
    for replicate, signature in signatures.items():
        differing = [key for key in baseline if signature.get(key) != baseline.get(key)]
        if differing:
            raise ValueError(
                f"replicate {replicate} is incompatible with {baseline_replicate}; "
                f"differing contract field(s): {', '.join(differing)}"
            )
    return manifests, {
        "status": "passed",
        "baseline_replicate": baseline_replicate,
        "validated_replicates": list(signatures),
        "signature": baseline,
    }


def safe_float(value: object) -> float:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return float("nan")
    return result if math.isfinite(result) else float("nan")


def safe_int(value: object) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return 0


def fmt_number(value: float) -> str:
    return "" if not math.isfinite(value) else f"{value:.10g}"


def median_finite(values: Iterable[object]) -> float:
    finite = [value for item in values if math.isfinite(value := safe_float(item))]
    return float(statistics.median(finite)) if finite else float("nan")


def min_finite(values: Iterable[object]) -> float:
    finite = [value for item in values if math.isfinite(value := safe_float(item))]
    return min(finite) if finite else float("nan")


def load_result(path: Path, replicate_id: str) -> Dict[Tuple[str, str], Dict[str, str]]:
    required = {"gorf_id", "overlap_type", "classification", "q_BH", "lambda_hat", "core_reads"}
    rows: Dict[Tuple[str, str], Dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"{path} is missing column(s): {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            gorf_id = str(row.get("gorf_id", "")).strip()
            overlap_type = str(row.get("overlap_type", "")).strip()
            if not gorf_id or not overlap_type:
                raise ValueError(
                    f"{path}:{line_number} has empty gorf_id or overlap_type; "
                    "formal replicate integration requires the DM deduplication identity"
                )
            key = (gorf_id, overlap_type)
            if key in rows:
                raise ValueError(f"{path}:{line_number} duplicates candidate {gorf_id!r}/{overlap_type!r}")
            normalized = {column: str(row.get(column, "")) for column in (*IDENTITY_COLUMNS, *EVIDENCE_COLUMNS)}
            normalized["replicate_id"] = replicate_id
            normalized["candidate_key"] = f"{gorf_id}|{overlap_type}"
            if "primary_credible_call" not in (reader.fieldnames or []):
                normalized["primary_credible_call"] = str(
                    int(normalized["classification"] == PRIMARY_CREDIBLE_CLASS)
                )
            rows[key] = normalized
    return rows


def metadata_for_candidate(
    observations: Sequence[Mapping[str, str]],
) -> Tuple[Dict[str, str], str]:
    first = observations[0]
    metadata = {column: str(first.get(column, "")) for column in IDENTITY_COLUMNS}
    conflicts = []
    for column in IDENTITY_COLUMNS:
        values = {str(row.get(column, "")) for row in observations}
        if len(values) > 1:
            conflicts.append(column)
    return metadata, ",".join(conflicts)


def integrate(
    expected_replicates: Sequence[str],
    result_paths: Mapping[str, Path],
    min_replicates: int,
    min_fraction: float,
    fdr_threshold: float,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], Dict[str, object]]:
    expected = list(expected_replicates)
    if len(set(expected)) != len(expected):
        raise ValueError("--expected-replicate contains duplicates")
    if not expected:
        raise ValueError("at least one --expected-replicate is required")
    unknown = sorted(set(result_paths) - set(expected))
    if unknown:
        raise ValueError(f"result supplied for unexpected replicate(s): {', '.join(unknown)}")
    if not result_paths:
        raise ValueError("no readable replicate result was supplied")

    by_replicate = {
        replicate: load_result(path, replicate)
        for replicate, path in result_paths.items()
    }
    available = [replicate for replicate in expected if replicate in by_replicate]
    keys = sorted({key for rows in by_replicate.values() for key in rows})
    long_rows: List[Dict[str, object]] = []
    consensus_rows: List[Dict[str, object]] = []
    label_counts: Counter[str] = Counter()

    for key in keys:
        observations = [by_replicate[replicate][key] for replicate in available if key in by_replicate[replicate]]
        metadata, conflicts = metadata_for_candidate(observations)
        observed_by = {str(row["replicate_id"]): row for row in observations}
        primary_count = sum(safe_int(row.get("primary_credible_call", 0)) == 1 for row in observations)
        significant_count = sum(
            math.isfinite(q := safe_float(row.get("q_BH"))) and q <= fdr_threshold
            for row in observations
        )
        support_fraction = primary_count / len(available)
        if len(available) < min_replicates:
            label = "single_replicate_primary_credible" if primary_count else "not_primary_credible"
        elif primary_count >= min_replicates and support_fraction >= min_fraction:
            label = "reproducible_primary_credible"
        elif primary_count:
            label = "primary_credible_not_reproducible"
        else:
            label = "not_primary_credible"
        label_counts[label] += 1

        for replicate in available:
            if replicate not in observed_by:
                continue
            row = observed_by[replicate]
            long_rows.append(
                {
                    "candidate_key": row["candidate_key"],
                    "replicate_id": replicate,
                    **metadata,
                    **{column: row.get(column, "") for column in EVIDENCE_COLUMNS},
                }
            )

        consensus: Dict[str, object] = {
            "candidate_key": f"{key[0]}|{key[1]}",
            **metadata,
            "metadata_conflict_columns": conflicts,
            "n_replicates_expected": len(expected),
            "n_replicates_available": len(available),
            "n_replicates_observed": len(observations),
            "n_primary_credible": primary_count,
            "primary_credible_fraction_available": fmt_number(support_fraction),
            "n_bh_significant": significant_count,
            "min_q_BH": fmt_number(min_finite(row.get("q_BH") for row in observations)),
            "median_q_BH": fmt_number(median_finite(row.get("q_BH") for row in observations)),
            "median_lambda_hat": fmt_number(median_finite(row.get("lambda_hat") for row in observations)),
            "max_core_reads": max((safe_int(row.get("core_reads")) for row in observations), default=0),
            "consensus_classification": label,
            "primary_credible_replicates": ",".join(
                replicate for replicate in available
                if replicate in observed_by and safe_int(observed_by[replicate].get("primary_credible_call")) == 1
            ),
        }
        for replicate in expected:
            row = observed_by.get(replicate)
            prefix = f"replicate::{replicate}::"
            consensus[prefix + "present"] = int(row is not None)
            consensus[prefix + "classification"] = row.get("classification", "") if row else ""
            consensus[prefix + "primary_credible_call"] = row.get("primary_credible_call", "") if row else ""
            consensus[prefix + "q_BH"] = row.get("q_BH", "") if row else ""
            consensus[prefix + "lambda_hat"] = row.get("lambda_hat", "") if row else ""
            consensus[prefix + "core_reads"] = row.get("core_reads", "") if row else ""
        consensus_rows.append(consensus)

    summary: Dict[str, object] = {
        "program": PROGRAM,
        "version": VERSION,
        "integration_policy": "replicate_evidence_only_no_read_pooling_no_pvalue_combination",
        "candidate_identity": "gorf_id+overlap_type",
        "expected_replicates": ",".join(expected),
        "available_replicates": ",".join(available),
        "unavailable_replicates": ",".join(item for item in expected if item not in by_replicate),
        "n_replicates_expected": len(expected),
        "n_replicates_available": len(available),
        "min_replicates": min_replicates,
        "min_fraction": min_fraction,
        "fdr_threshold": fdr_threshold,
        "candidate_union_total": len(consensus_rows),
        "observed_candidate_replicate_rows": len(long_rows),
    }
    for label in (
        "reproducible_primary_credible", "single_replicate_primary_credible",
        "primary_credible_not_reproducible", "not_primary_credible",
    ):
        summary[f"consensus::{label}"] = label_counts[label]
    return long_rows, consensus_rows, summary


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{path}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog=PROGRAM, description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--expected-replicate", action="append", required=True)
    parser.add_argument("--result", action="append", default=[], metavar="REPLICATE_ID=DM_RESULTS.tsv")
    parser.add_argument(
        "--run-manifest", action="append", default=[], metavar="REPLICATE_ID=RUN_MANIFEST.json",
        help="DM run manifest; inferred from a standard result filename when omitted",
    )
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--min-replicates", type=int, default=2)
    parser.add_argument("--min-fraction", type=float, default=0.5)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.min_replicates < 2:
        raise SystemExit("ERROR: --min-replicates must be >= 2")
    if not 0 < args.min_fraction <= 1:
        raise SystemExit("ERROR: --min-fraction must be in (0, 1]")
    if not 0 < args.fdr_threshold < 1:
        raise SystemExit("ERROR: --fdr-threshold must be in (0, 1)")
    try:
        expected = [str(item).strip() for item in args.expected_replicate]
        assignments = [parse_assignment(value, "--result") for value in args.result]
        manifest_assignments = [
            parse_assignment(value, "--run-manifest") for value in args.run_manifest
        ]
        result_paths: Dict[str, Path] = {}
        for replicate, value in assignments:
            if replicate in result_paths:
                raise ValueError(f"duplicate --result for replicate {replicate!r}")
            path = Path(value)
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"result is missing or empty for {replicate}: {path}")
            result_paths[replicate] = path
        manifest_paths: Dict[str, Path] = {}
        for replicate, value in manifest_assignments:
            if replicate in manifest_paths:
                raise ValueError(f"duplicate --run-manifest for replicate {replicate!r}")
            manifest_paths[replicate] = Path(value)
        unknown_manifests = sorted(set(manifest_paths) - set(result_paths))
        if unknown_manifests:
            raise ValueError(
                "run manifest supplied without a result for replicate(s): "
                + ", ".join(unknown_manifests)
            )
        for replicate, result_path in result_paths.items():
            if replicate not in manifest_paths:
                manifest_paths[replicate] = default_run_manifest_path(result_path)
        _run_manifests, compatibility = validate_manifest_compatibility(
            expected, result_paths, manifest_paths
        )
        long_rows, consensus_rows, summary = integrate(
            expected, result_paths, args.min_replicates, args.min_fraction, args.fdr_threshold
        )
        summary["compatibility_status"] = compatibility["status"]
        summary["compatibility_baseline_replicate"] = compatibility["baseline_replicate"]
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

    prefix = Path(args.out_prefix)
    long_path = Path(f"{prefix}.replicate_long.tsv")
    consensus_path = Path(f"{prefix}.consensus.tsv")
    summary_path = Path(f"{prefix}.summary.tsv")
    manifest_path = Path(f"{prefix}.manifest.json")
    long_fields = ["candidate_key", "replicate_id", *IDENTITY_COLUMNS, *EVIDENCE_COLUMNS]
    consensus_fields = list(consensus_rows[0]) if consensus_rows else [
        "candidate_key", *IDENTITY_COLUMNS, "metadata_conflict_columns",
        "n_replicates_expected", "n_replicates_available", "n_replicates_observed",
        "n_primary_credible", "primary_credible_fraction_available", "n_bh_significant",
        "min_q_BH", "median_q_BH", "median_lambda_hat", "max_core_reads",
        "consensus_classification", "primary_credible_replicates",
    ]
    write_tsv(long_path, long_rows, long_fields)
    write_tsv(consensus_path, consensus_rows, consensus_fields)
    write_tsv(
        summary_path,
        [{"metric": key, "value": value} for key, value in summary.items()],
        ["metric", "value"],
    )
    manifest = {
        **summary,
        "compatibility": compatibility,
        "result_inputs": {
            key: {
                "result": file_identity(result_paths[key]),
                "run_manifest": file_identity(manifest_paths[key]),
            }
            for key in result_paths
        },
        "outputs": {
            "replicate_long": str(long_path),
            "consensus": str(consensus_path),
            "summary": str(summary_path),
        },
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_temporary = Path(f"{manifest_path}.tmp")
    with manifest_temporary.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    manifest_temporary.replace(manifest_path)
    print(f"[ok] replicate_long={long_path}")
    print(f"[ok] consensus={consensus_path}")
    print(f"[ok] summary={summary_path}")
    print(f"[ok] manifest={manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
