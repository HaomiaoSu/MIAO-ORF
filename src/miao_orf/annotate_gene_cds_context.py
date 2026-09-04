#!/usr/bin/env python3
"""Annotate called ORFs with same-gene annotated-CDS N-terminal reuse.

This is deliberately a postprocessing layer.  It preserves every input result
column and never changes transcript-level ORF types, DM evidence, p/q values,
lambda estimates, or primary-call decisions.  A match requires the candidate
and an annotated CDS from the same versioned gene and strand to share the
translation-oriented genomic coordinate path from the candidate start at an
annotated-CDS codon boundary.  Interval overlap alone is not considered a
match.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


PROGRAM = "miao-orf-gene-cds-context"
VERSION = "1.0.0"

CONTEXT_FIELDS = (
    "gene_level_orf_class",
    "gene_level_pure_intorf_eligible",
    "gene_cds_nterm_match",
    "gene_cds_nterm_context",
    "gene_cds_nterm_candidate_torf_id",
    "gene_cds_nterm_best_transcript_id",
    "gene_cds_nterm_best_transcript_name",
    "gene_cds_nterm_best_transcript_tags",
    "gene_cds_nterm_tied_transcript_ids",
    "gene_cds_nterm_coordinate_prefix_nt",
    "gene_cds_nterm_coordinate_prefix_complete_codons",
    "gene_cds_nterm_coordinate_prefix_fraction",
    "gene_cds_nterm_peptide_prefix_aa",
    "gene_cds_nterm_peptide_prefix_fraction",
    "gene_cds_nterm_break_reason",
    "gene_cds_nterm_min_prefix_codons",
)


@dataclass(frozen=True)
class TorfRecord:
    torf_id: str
    gorf_id: str
    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str
    transcript_tags: str
    chrom: str
    strand: str
    peptide: str
    peptide_len: int
    is_annotated_cds: bool
    coding_path: Tuple[int, ...]


@dataclass(frozen=True)
class Match:
    candidate: TorfRecord
    cds: TorfRecord
    cds_start_index: int
    coordinate_prefix_nt: int
    peptide_prefix_aa: int
    break_reason: str

    @property
    def complete_codons(self) -> int:
        return self.coordinate_prefix_nt // 3


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog=PROGRAM,
        description=(
            "Add an auditable same-gene annotated-CDS N-terminal reuse layer "
            "to an existing candidate/result TSV."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--torf", required=True, help="Full tORF TSV from MIAO orfscan")
    parser.add_argument(
        "--input-tsv",
        required=True,
        help="DM result or replicate-consensus TSV containing gorf_id and gene_id",
    )
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument(
        "--min-prefix-codons",
        type=int,
        default=5,
        help="Minimum coordinate-identical complete N-terminal codons for a reported reuse match",
    )
    return parser.parse_args(argv)


def split_csv_ints(value: str, label: str) -> List[int]:
    try:
        result = [int(token.strip()) for token in str(value).split(",") if token.strip()]
    except ValueError as exc:
        raise ValueError(f"invalid integer list in {label}: {value!r}") from exc
    if not result:
        raise ValueError(f"empty integer list in {label}")
    return result


def truthy(value: object) -> bool:
    return str(value).strip().casefold() in {"1", "true", "yes", "y"}


def translation_path(row: Mapping[str, str], line_number: int) -> Tuple[int, ...]:
    starts = split_csv_ints(row.get("genomic_block_starts1", ""), "genomic_block_starts1")
    ends = split_csv_ints(row.get("genomic_block_ends1", ""), "genomic_block_ends1")
    if len(starts) != len(ends):
        raise ValueError(f"tORF line {line_number}: block start/end count mismatch")
    blocks = []
    for start, end in zip(starts, ends):
        if start > end:
            raise ValueError(f"tORF line {line_number}: block start exceeds end")
        blocks.append((start, end))
    strand = str(row.get("strand", "")).strip()
    if strand == "+":
        ordered = sorted(blocks)
        path = [coordinate for start, end in ordered for coordinate in range(start, end + 1)]
    elif strand == "-":
        ordered = sorted(blocks, reverse=True)
        path = [coordinate for start, end in ordered for coordinate in range(end, start - 1, -1)]
    else:
        raise ValueError(f"tORF line {line_number}: invalid strand {strand!r}")
    try:
        peptide_len = int(str(row.get("peptide_len", "")).strip())
    except ValueError as exc:
        raise ValueError(f"tORF line {line_number}: invalid peptide_len") from exc
    coding_nt = peptide_len * 3
    if peptide_len < 1 or len(path) < coding_nt:
        raise ValueError(
            f"tORF line {line_number}: genomic path has {len(path)} nt but peptide needs {coding_nt} nt"
        )
    return tuple(path[:coding_nt])


def read_input(path: Path) -> Tuple[List[str], List[Dict[str, str]], set[str], set[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        required = {"gorf_id", "gene_id"}
        missing = sorted(required - set(fields))
        if missing:
            raise ValueError(f"input TSV is missing column(s): {', '.join(missing)}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError("input TSV contains no candidate rows")
    target_gorfs = {str(row.get("gorf_id", "")).strip() for row in rows}
    target_genes = {str(row.get("gene_id", "")).strip() for row in rows}
    if "" in target_gorfs or "" in target_genes:
        raise ValueError("input TSV contains blank gorf_id or gene_id")
    return fields, rows, target_gorfs, target_genes


def to_torf_record(row: Mapping[str, str], line_number: int) -> TorfRecord:
    peptide = str(row.get("peptide", "")).strip()
    peptide_len = int(str(row.get("peptide_len", "")).strip())
    if len(peptide) != peptide_len:
        raise ValueError(
            f"tORF line {line_number}: peptide length {len(peptide)} != peptide_len {peptide_len}"
        )
    return TorfRecord(
        torf_id=str(row.get("torf_id", "")).strip(),
        gorf_id=str(row.get("gorf_id", "")).strip(),
        gene_id=str(row.get("gene_id", "")).strip(),
        gene_name=str(row.get("gene_name", "")).strip(),
        transcript_id=str(row.get("transcript_id", "")).strip(),
        transcript_name=str(row.get("transcript_name", "")).strip(),
        transcript_tags=str(row.get("tx_tags", "")).strip(),
        chrom=str(row.get("chrom", "")).strip(),
        strand=str(row.get("strand", "")).strip(),
        peptide=peptide,
        peptide_len=peptide_len,
        is_annotated_cds=truthy(row.get("is_annotated_cds", "")),
        coding_path=translation_path(row, line_number),
    )


def read_relevant_torfs(
    path: Path,
    target_gorfs: set[str],
    target_genes: set[str],
) -> Tuple[Dict[str, List[TorfRecord]], Dict[str, List[TorfRecord]], int]:
    candidate_by_gorf: Dict[str, List[TorfRecord]] = {}
    cds_by_gene: Dict[str, List[TorfRecord]] = {}
    total = 0
    required = {
        "torf_id", "gorf_id", "gene_id", "transcript_id", "chrom", "strand",
        "is_annotated_cds", "genomic_block_starts1", "genomic_block_ends1",
        "peptide_len", "peptide",
    }
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"tORF TSV is missing column(s): {', '.join(missing)}")
        for line_number, row in enumerate(reader, start=2):
            total += 1
            gorf_id = str(row.get("gorf_id", "")).strip()
            gene_id = str(row.get("gene_id", "")).strip()
            is_cds = truthy(row.get("is_annotated_cds", ""))
            if not ((gorf_id in target_gorfs and not is_cds) or (gene_id in target_genes and is_cds)):
                continue
            record = to_torf_record(row, line_number)
            if record.is_annotated_cds:
                cds_by_gene.setdefault(record.gene_id, []).append(record)
            else:
                candidate_by_gorf.setdefault(record.gorf_id, []).append(record)
    return candidate_by_gorf, cds_by_gene, total


def annotation_priority(cds: TorfRecord) -> Tuple[int, str]:
    tags = cds.transcript_tags.casefold()
    if "mane_select" in tags:
        rank = 0
    elif "ccds" in tags:
        rank = 1
    elif "ensembl_canonical" in tags or "gencode_primary" in tags:
        rank = 2
    elif "appris" in tags:
        rank = 3
    else:
        rank = 4
    return rank, cds.transcript_id


def common_prefix(left: Sequence[object], right: Sequence[object]) -> int:
    count = 0
    for left_value, right_value in zip(left, right):
        if left_value != right_value:
            break
        count += 1
    return count


def compare(candidate: TorfRecord, cds: TorfRecord) -> Optional[Match]:
    if (
        candidate.gene_id != cds.gene_id
        or candidate.chrom != cds.chrom
        or candidate.strand != cds.strand
        or not candidate.coding_path
    ):
        return None
    try:
        cds_start_index = cds.coding_path.index(candidate.coding_path[0])
    except ValueError:
        return None
    if cds_start_index % 3 != 0:
        return None
    cds_suffix = cds.coding_path[cds_start_index:]
    coordinate_prefix_nt = common_prefix(candidate.coding_path, cds_suffix)
    if coordinate_prefix_nt == 0:
        return None
    cds_peptide_suffix = cds.peptide[cds_start_index // 3:]
    peptide_prefix_aa = common_prefix(candidate.peptide, cds_peptide_suffix)
    if coordinate_prefix_nt >= len(candidate.coding_path):
        reason = "candidate_fully_contained_in_annotated_cds_path"
    elif coordinate_prefix_nt >= len(cds_suffix):
        reason = "annotated_cds_path_ended"
    else:
        reason = "splice_or_coordinate_path_divergence"
    return Match(
        candidate=candidate,
        cds=cds,
        cds_start_index=cds_start_index,
        coordinate_prefix_nt=coordinate_prefix_nt,
        peptide_prefix_aa=peptide_prefix_aa,
        break_reason=reason,
    )


def empty_context(min_prefix_codons: int, candidate_torf_id: str = "") -> Dict[str, object]:
    return {
        "gene_level_orf_class": "pure_intorf_no_annotated_cds_nterm_reuse",
        "gene_level_pure_intorf_eligible": 1,
        "gene_cds_nterm_match": 0,
        "gene_cds_nterm_context": "none",
        "gene_cds_nterm_candidate_torf_id": candidate_torf_id,
        "gene_cds_nterm_best_transcript_id": "",
        "gene_cds_nterm_best_transcript_name": "",
        "gene_cds_nterm_best_transcript_tags": "",
        "gene_cds_nterm_tied_transcript_ids": "",
        "gene_cds_nterm_coordinate_prefix_nt": 0,
        "gene_cds_nterm_coordinate_prefix_complete_codons": 0,
        "gene_cds_nterm_coordinate_prefix_fraction": "0",
        "gene_cds_nterm_peptide_prefix_aa": 0,
        "gene_cds_nterm_peptide_prefix_fraction": "0",
        "gene_cds_nterm_break_reason": "no_same_gene_same_strand_codon_aligned_coordinate_prefix",
        "gene_cds_nterm_min_prefix_codons": min_prefix_codons,
    }


def context_for_candidate(
    candidates: Sequence[TorfRecord],
    annotated_cdss: Sequence[TorfRecord],
    min_prefix_codons: int,
) -> Dict[str, object]:
    candidate_torf_id = candidates[0].torf_id if candidates else ""
    context = empty_context(min_prefix_codons, candidate_torf_id)
    matches = [
        match
        for candidate in candidates
        for cds in annotated_cdss
        for match in [compare(candidate, cds)]
        if match is not None
    ]
    if not matches:
        if not candidates:
            context["gene_cds_nterm_break_reason"] = "candidate_gorf_absent_from_torf"
        elif not annotated_cdss:
            context["gene_cds_nterm_break_reason"] = "no_annotated_cds_for_same_versioned_gene"
        return context
    matches.sort(
        key=lambda item: (
            -item.coordinate_prefix_nt,
            -item.peptide_prefix_aa,
            annotation_priority(item.cds),
            item.candidate.torf_id,
        )
    )
    best = matches[0]
    tied = sorted({
        item.cds.transcript_id
        for item in matches
        if item.coordinate_prefix_nt == best.coordinate_prefix_nt
        and item.peptide_prefix_aa == best.peptide_prefix_aa
    })
    coordinate_fraction = best.coordinate_prefix_nt / len(best.candidate.coding_path)
    peptide_fraction = best.peptide_prefix_aa / best.candidate.peptide_len
    significant = best.complete_codons >= min_prefix_codons
    full_coordinate = best.coordinate_prefix_nt >= len(best.candidate.coding_path)
    full_peptide = best.peptide_prefix_aa >= best.candidate.peptide_len
    if significant and full_coordinate:
        level_class = "annotated_cds_derived_full_coordinate_path"
        level_context = "full_coordinate_path"
    elif significant and full_peptide:
        level_class = "annotated_cds_full_peptide_alternative_coordinate_path"
        level_context = "full_peptide_alternative_path"
    elif significant:
        level_class = "annotated_cds_nterm_reuse_with_splice_derived_cterm"
        level_context = "partial_nterm_reuse"
    else:
        level_class = "pure_intorf_no_annotated_cds_nterm_reuse"
        level_context = "below_minimum_prefix"
    context.update({
        "gene_level_orf_class": level_class,
        "gene_level_pure_intorf_eligible": 0 if significant else 1,
        "gene_cds_nterm_match": 1 if significant else 0,
        "gene_cds_nterm_context": level_context,
        "gene_cds_nterm_candidate_torf_id": best.candidate.torf_id,
        "gene_cds_nterm_best_transcript_id": best.cds.transcript_id,
        "gene_cds_nterm_best_transcript_name": best.cds.transcript_name,
        "gene_cds_nterm_best_transcript_tags": best.cds.transcript_tags,
        "gene_cds_nterm_tied_transcript_ids": ",".join(tied),
        "gene_cds_nterm_coordinate_prefix_nt": best.coordinate_prefix_nt,
        "gene_cds_nterm_coordinate_prefix_complete_codons": best.complete_codons,
        "gene_cds_nterm_coordinate_prefix_fraction": f"{coordinate_fraction:.12g}",
        "gene_cds_nterm_peptide_prefix_aa": best.peptide_prefix_aa,
        "gene_cds_nterm_peptide_prefix_fraction": f"{peptide_fraction:.12g}",
        "gene_cds_nterm_break_reason": best.break_reason,
    })
    return context


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(path: Path) -> Dict[str, object]:
    stat = path.stat()
    return {
        "path": str(path.absolute()),
        "size_bytes": stat.st_size,
        "sha256": sha256_file(path),
    }


def write_tsv(path: Path, rows: Iterable[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{path}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def run(args: argparse.Namespace) -> Dict[str, Path]:
    if args.min_prefix_codons < 1:
        raise ValueError("--min-prefix-codons must be >= 1")
    torf_path = Path(args.torf)
    input_path = Path(args.input_tsv)
    for path in (torf_path, input_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing/empty input: {path}")
    input_fields, input_rows, target_gorfs, target_genes = read_input(input_path)
    candidate_by_gorf, cds_by_gene, torf_rows_scanned = read_relevant_torfs(
        torf_path, target_gorfs, target_genes
    )
    candidate_keys = {
        (str(row["gorf_id"]).strip(), str(row["gene_id"]).strip())
        for row in input_rows
    }
    missing_candidate_keys = sorted(
        (gorf_id, gene_id)
        for gorf_id, gene_id in candidate_keys
        if not any(
            candidate.gene_id == gene_id
            for candidate in candidate_by_gorf.get(gorf_id, [])
        )
    )
    if missing_candidate_keys:
        examples = ", ".join(
            f"{gorf_id}/{gene_id}" for gorf_id, gene_id in missing_candidate_keys[:10]
        )
        raise ValueError(
            f"{len(missing_candidate_keys)} input gORF/gene key(s) are absent from the supplied tORF; "
            f"the result/reference pair is incompatible (examples: {examples})"
        )
    context_by_key: Dict[Tuple[str, str], Dict[str, object]] = {}
    for row in input_rows:
        gorf_id = str(row["gorf_id"]).strip()
        gene_id = str(row["gene_id"]).strip()
        key = (gorf_id, gene_id)
        if key in context_by_key:
            continue
        context_by_key[key] = context_for_candidate(
            [
                candidate
                for candidate in candidate_by_gorf.get(gorf_id, [])
                if candidate.gene_id == gene_id
            ],
            cds_by_gene.get(gene_id, []),
            args.min_prefix_codons,
        )
    annotated_rows = [
        {
            **row,
            **context_by_key[(str(row["gorf_id"]).strip(), str(row["gene_id"]).strip())],
        }
        for row in input_rows
    ]
    compact_fields = [
        field for field in ("gorf_id", "overlap_type", "gene_id", "gene_name", "chrom", "strand")
        if field in input_fields
    ] + list(CONTEXT_FIELDS)
    compact_rows = [
        {
            **{field: row.get(field, "") for field in compact_fields},
            **context_by_key[(str(row["gorf_id"]).strip(), str(row["gene_id"]).strip())],
        }
        for row in input_rows
    ]
    class_counts: Dict[str, int] = {}
    primary_class_counts: Dict[str, int] = {}
    for row in annotated_rows:
        level_class = str(row["gene_level_orf_class"])
        class_counts[level_class] = class_counts.get(level_class, 0) + 1
        if truthy(row.get("primary_credible_call", "")):
            primary_class_counts[level_class] = primary_class_counts.get(level_class, 0) + 1
    summary: List[Dict[str, object]] = [
        {"metric": "program", "value": PROGRAM},
        {"metric": "version", "value": VERSION},
        {"metric": "candidate_rows", "value": len(input_rows)},
        {"metric": "unique_gorfs", "value": len(target_gorfs)},
        {"metric": "unique_gorf_gene_keys", "value": len(context_by_key)},
        {"metric": "torf_rows_scanned", "value": torf_rows_scanned},
        {"metric": "min_prefix_codons", "value": args.min_prefix_codons},
    ]
    for level_class in sorted(class_counts):
        summary.append({"metric": f"gene_level_orf_class::{level_class}", "value": class_counts[level_class]})
    for level_class in sorted(primary_class_counts):
        summary.append({"metric": f"primary_gene_level_orf_class::{level_class}", "value": primary_class_counts[level_class]})

    prefix = Path(args.out_prefix)
    outputs = {
        "annotated": Path(f"{prefix}.gene_cds_context.tsv"),
        "context_only": Path(f"{prefix}.gene_cds_context_only.tsv"),
        "summary": Path(f"{prefix}.gene_cds_context_summary.tsv"),
        "manifest": Path(f"{prefix}.gene_cds_context_manifest.json"),
    }
    annotated_fields = [field for field in input_fields if field not in CONTEXT_FIELDS]
    write_tsv(outputs["annotated"], annotated_rows, [*annotated_fields, *CONTEXT_FIELDS])
    write_tsv(outputs["context_only"], compact_rows, compact_fields)
    write_tsv(outputs["summary"], summary, ("metric", "value"))
    manifest = {
        "program": PROGRAM,
        "version": VERSION,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "algorithm": {
            "gene_identity": "exact versioned gene_id",
            "strand_required": True,
            "candidate_start_requires_annotated_cds_codon_boundary": True,
            "comparison": "translation-oriented genomic coordinate path",
            "interval_overlap_alone_is_match": False,
            "min_prefix_codons": args.min_prefix_codons,
            "transcript_level_classification_modified": False,
            "dm_statistics_modified": False,
        },
        "inputs": {"torf": file_identity(torf_path), "candidate_table": file_identity(input_path)},
        "outputs": {key: str(path.absolute()) for key, path in outputs.items()},
        "counts": {
            "candidate_rows": len(input_rows),
            "unique_gorfs": len(target_gorfs),
            "unique_gorf_gene_keys": len(context_by_key),
            "class_counts": class_counts,
            "primary_class_counts": primary_class_counts,
        },
    }
    outputs["manifest"].parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{outputs['manifest']}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    os.replace(temporary, outputs["manifest"])
    return outputs


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        outputs = run(args)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    for label, path in outputs.items():
        print(f"[ok] {label}={path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
