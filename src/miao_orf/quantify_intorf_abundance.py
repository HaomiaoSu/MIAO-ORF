#!/usr/bin/env python3
"""Quantify model-allocated intORF translation abundance from DM results.

This is a read-only post-processing stage.  It does not refit the DM model or
change any call, p-value, FDR result, threshold, or classification.  For each
DM row with a model-expected intORF component it reports

    intorf_pFPKM = 1e9 * model_expected_intorf_core_reads
                    / (effective_core_nt * usable_library_psites)

where ``effective_core_nt`` is ``3 * n_core_codons`` and the library size is
the number of mapped 1-nt P-site alignments in the sample P-site BAM.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import statistics
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


PROGRAM = "miao-orf-quantify-intorf-abundance"
VERSION = "1.0.0"
OUTPUT_SCHEMA_VERSION = "1.0"
FORMULA = (
    "intorf_pFPKM = 1e9 * model_expected_intorf_core_reads / "
    "(3 * n_core_codons * mapped_1nt_psite_alignments)"
)

IDENTITY_COLUMNS = [
    "torf_id", "gorf_id", "gene_id", "gene_name", "transcript_id",
    "chrom", "strand", "overlap_type", "classification",
    "primary_credible_call", "qc_status", "filter_reason",
]
OUTPUT_COLUMNS = IDENTITY_COLUMNS + [
    "n_core_codons", "effective_core_nt", "usable_library_psites",
    "analyzed_core_reads", "lambda_hat", "lambda_profile_ci95_low",
    "lambda_profile_ci95_high", "model_expected_intorf_core_reads",
    "model_expected_host_cds_core_reads", "observed_core_pFPKM",
    "intorf_psite_RPM", "intorf_pFPKM", "intorf_pFPKM_ci95_low",
    "intorf_pFPKM_ci95_high", "host_component_pFPKM",
    "intorf_to_host_ratio", "abundance_status",
]
REQUIRED_DM_COLUMNS = {
    "gorf_id", "n_core_codons", "analyzed_core_reads", "lambda_hat",
    "lambda_profile_ci95_low", "lambda_profile_ci95_high",
    "model_expected_intorf_core_reads",
    "model_expected_host_cds_core_reads",
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog=PROGRAM, description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--dm-results", required=True, help="DM result TSV")
    parser.add_argument(
        "--psite-bam",
        required=True,
        help="Sample 1-nt P-site BAM used by the DM caller",
    )
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--sample", default="", help="Optional sample label for provenance")
    return parser.parse_args(argv)


def finite_float(value: object) -> float:
    try:
        result = float(str(value).strip())
    except (TypeError, ValueError):
        return float("nan")
    return result if math.isfinite(result) else float("nan")


def finite_int(value: object) -> int | None:
    number = finite_float(value)
    if not math.isfinite(number) or number != math.floor(number):
        return None
    return int(number)


def format_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.12g}" if math.isfinite(value) else ""
    return str(value)


def abundance_value(expected_reads: float, effective_nt: int, library_psites: int) -> float:
    return 1.0e9 * expected_reads / (effective_nt * library_psites)


def quantify_row(row: Mapping[str, str], library_psites: int) -> Dict[str, object]:
    output: Dict[str, object] = {column: row.get(column, "") for column in IDENTITY_COLUMNS}
    output.update({column: "" for column in OUTPUT_COLUMNS if column not in output})
    output["usable_library_psites"] = library_psites

    n_core_codons = finite_int(row.get("n_core_codons", ""))
    output["n_core_codons"] = "" if n_core_codons is None else n_core_codons
    if n_core_codons is None or n_core_codons <= 0:
        output["abundance_status"] = "invalid_core_length"
        return output

    effective_nt = 3 * n_core_codons
    output["effective_core_nt"] = effective_nt
    analyzed_reads = finite_float(row.get("analyzed_core_reads", ""))
    lambda_hat = finite_float(row.get("lambda_hat", ""))
    lambda_low = finite_float(row.get("lambda_profile_ci95_low", ""))
    lambda_high = finite_float(row.get("lambda_profile_ci95_high", ""))
    intorf_reads = finite_float(row.get("model_expected_intorf_core_reads", ""))
    host_reads = finite_float(row.get("model_expected_host_cds_core_reads", ""))

    for name, value in (
        ("analyzed_core_reads", analyzed_reads),
        ("lambda_hat", lambda_hat),
        ("lambda_profile_ci95_low", lambda_low),
        ("lambda_profile_ci95_high", lambda_high),
        ("model_expected_intorf_core_reads", intorf_reads),
        ("model_expected_host_cds_core_reads", host_reads),
    ):
        output[name] = value

    if math.isfinite(analyzed_reads) and analyzed_reads >= 0:
        output["observed_core_pFPKM"] = abundance_value(
            analyzed_reads, effective_nt, library_psites
        )

    if not math.isfinite(intorf_reads):
        output["abundance_status"] = "not_model_quantifiable"
        return output
    if intorf_reads < 0 or (math.isfinite(host_reads) and host_reads < 0):
        output["abundance_status"] = "invalid_model_expected_reads"
        return output

    output["intorf_psite_RPM"] = 1.0e6 * intorf_reads / library_psites
    output["intorf_pFPKM"] = abundance_value(intorf_reads, effective_nt, library_psites)
    if math.isfinite(host_reads):
        output["host_component_pFPKM"] = abundance_value(
            host_reads, effective_nt, library_psites
        )
        if host_reads > 0:
            output["intorf_to_host_ratio"] = intorf_reads / host_reads

    if (
        math.isfinite(analyzed_reads)
        and analyzed_reads >= 0
        and math.isfinite(lambda_low)
        and math.isfinite(lambda_high)
        and 0 <= lambda_low <= lambda_high <= 1
    ):
        output["intorf_pFPKM_ci95_low"] = abundance_value(
            analyzed_reads * lambda_low, effective_nt, library_psites
        )
        output["intorf_pFPKM_ci95_high"] = abundance_value(
            analyzed_reads * lambda_high, effective_nt, library_psites
        )

    output["abundance_status"] = "quantified"
    return output


def quantify_rows(
    rows: Iterable[Mapping[str, str]], library_psites: int
) -> List[Dict[str, object]]:
    if library_psites <= 0:
        raise ValueError("usable P-site library size must be positive")
    return [quantify_row(row, library_psites) for row in rows]


def read_dm_results(path: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        missing = sorted(REQUIRED_DM_COLUMNS - set(fieldnames))
        if missing:
            raise ValueError(
                "DM result table is missing required column(s): " + ", ".join(missing)
            )
        rows = list(reader)
    return rows, fieldnames


def count_mapped_psites(path: Path) -> Tuple[int, str]:
    """Count mapped alignments in the already filtered 1-nt P-site BAM."""
    try:
        import pysam  # type: ignore
    except ImportError as exc:
        raise RuntimeError("pysam is required to inspect the P-site BAM") from exc

    with pysam.AlignmentFile(str(path), "rb") as bam:
        try:
            mapped = int(bam.mapped)
            method = "BAM_index_mapped_alignments"
        except (ValueError, OSError):
            mapped = sum(1 for read in bam.fetch(until_eof=True) if not read.is_unmapped)
            method = "sequential_mapped_alignment_scan"
    if mapped <= 0:
        raise ValueError("P-site BAM contains no mapped alignments")
    return mapped, method


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(path: Path, include_hash: bool) -> Dict[str, object]:
    stat = path.stat()
    result: Dict[str, object] = {
        "path": str(path.resolve()),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    result["sha256"] = sha256_file(path) if include_hash else None
    result["sha256_status"] = "computed" if include_hash else "not_computed_large_binary"
    return result


def atomic_write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{path}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: format_value(row.get(column, "")) for column in columns})
    os.replace(temporary, path)


def atomic_write_json(path: Path, value: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{path}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def summary_rows(
    sample: str,
    quantified: Sequence[Mapping[str, object]],
    library_psites: int,
    count_method: str,
) -> List[Dict[str, object]]:
    statuses = Counter(str(row.get("abundance_status", "")) for row in quantified)
    values = [
        finite_float(row.get("intorf_pFPKM", ""))
        for row in quantified
        if str(row.get("abundance_status", "")) == "quantified"
    ]
    values = [value for value in values if math.isfinite(value)]
    metrics: List[Tuple[str, object]] = [
        ("program", PROGRAM),
        ("version", VERSION),
        ("output_schema_version", OUTPUT_SCHEMA_VERSION),
        ("sample", sample),
        ("input_rows", len(quantified)),
        ("usable_library_psites", library_psites),
        ("library_count_method", count_method),
        ("quantified_rows", statuses.get("quantified", 0)),
        ("not_model_quantifiable_rows", statuses.get("not_model_quantifiable", 0)),
        ("invalid_core_length_rows", statuses.get("invalid_core_length", 0)),
        ("invalid_model_expected_reads_rows", statuses.get("invalid_model_expected_reads", 0)),
        ("median_intorf_pFPKM", statistics.median(values) if values else ""),
        ("max_intorf_pFPKM", max(values) if values else ""),
        ("formula", FORMULA),
        ("numerator_semantics", "DM model-expected intORF-like P-sites in analyzed core"),
        ("effective_length_semantics", "3 * DM n_core_codons"),
        ("library_semantics", "mapped alignments in the sample 1-nt P-site BAM"),
    ]
    return [{"metric": key, "value": value} for key, value in metrics]


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    dm_results = Path(args.dm_results)
    psite_bam = Path(args.psite_bam)
    if not dm_results.is_file() or dm_results.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty DM results: {dm_results}")
    if not psite_bam.is_file() or psite_bam.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty P-site BAM: {psite_bam}")

    rows, input_columns = read_dm_results(dm_results)
    library_psites, count_method = count_mapped_psites(psite_bam)
    quantified = quantify_rows(rows, library_psites)

    prefix = Path(args.out_prefix)
    table_path = Path(f"{prefix}.intorf_abundance.tsv")
    summary_path = Path(f"{prefix}.intorf_abundance_summary.tsv")
    manifest_path = Path(f"{prefix}.intorf_abundance_manifest.json")
    atomic_write_tsv(table_path, OUTPUT_COLUMNS, quantified)
    summaries = summary_rows(args.sample, quantified, library_psites, count_method)
    atomic_write_tsv(summary_path, ["metric", "value"], summaries)

    statuses = Counter(str(row.get("abundance_status", "")) for row in quantified)
    manifest: Dict[str, object] = {
        "program": PROGRAM,
        "version": VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "sample": args.sample,
        "formula": FORMULA,
        "semantics": {
            "post_processing_only": True,
            "dm_model_refit": False,
            "dm_calls_or_thresholds_changed": False,
            "numerator": "model_expected_intorf_core_reads",
            "effective_length_nt": "3 * n_core_codons",
            "library_size": "mapped alignments in the sample 1-nt P-site BAM",
            "length_subset": "inherited from the supplied P-site BAM",
        },
        "inputs": {
            "dm_results": file_identity(dm_results, include_hash=True),
            "psite_bam": file_identity(psite_bam, include_hash=False),
        },
        "input_dm_columns": input_columns,
        "output_columns": OUTPUT_COLUMNS,
        "usable_library_psites": library_psites,
        "library_count_method": count_method,
        "result_rows": len(quantified),
        "status_counts": dict(sorted(statuses.items())),
        "outputs": {
            "abundance_table": file_identity(table_path, include_hash=True),
            "summary": file_identity(summary_path, include_hash=True),
        },
    }
    atomic_write_json(manifest_path, manifest)
    print(f"Wrote {table_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
