#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Sample-level Ribo-seq metagene quality control over MANE CDS annotations.

The script reads a 1-nt P-site BAM and one or more tORF TSV files produced by
``orf_scan_transcriptome.py``. It selects annotated MANE CDS records, evaluates
sample-level 3-nt periodicity, produces start/stop metagene profiles, writes
per-CDS QC metrics, and estimates a canonical translation template
``pi_trans`` together with a Dirichlet-multinomial concentration parameter
``A0``.

This tool is intended for sample-level QC. It does not call intORFs.

Author: Haomiao Su
Contact: suhaomiao@csu.edu.cn
License: MIT
"""

from __future__ import annotations

import argparse
import collections
import csv
import glob
import math
import os
import sys
from dataclasses import dataclass
from multiprocessing import Pool, cpu_count
from typing import Any, Counter, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, TYPE_CHECKING

PROGRAM = "miao-orf-metagene-qc"
VERSION = "1.0.0"
AUTHOR = "Haomiao Su"
CONTACT = "suhaomiao@csu.edu.cn"
CODON_TERNARY_GRID_PERCENT = 5.0
CODON_TERNARY_RESOLUTION = int(round(100.0 / CODON_TERNARY_GRID_PERCENT))
CODON_TERNARY_AREA_PER_PERCENT = 18.0

if TYPE_CHECKING:
    import numpy as np
    import pysam
else:
    np = Any
    pysam = Any

# Populated after dependency validation.
_np = None
_pysam = None
_plt = None
_gammaln = None
_minimize_scalar = None


@dataclass(frozen=True)
class CDSRecord:
    torf_id: str
    gene_id: str
    gene_name: str
    transcript_id: str
    chrom: str
    strand: str
    blocks: Tuple[Tuple[int, int], ...]
    flank5: Tuple[Tuple[int, int], ...]
    flank3: Tuple[Tuple[int, int], ...]
    start_anchor: int
    stop_anchor: int

    @property
    def orf_len_nt(self) -> int:
        return sum(e - s for s, e in self.blocks)


@dataclass
class TemplateStats:
    sum_counts: List[int]
    sum_props: List[float]
    patterns: collections.Counter
    n_codons: int
    n_cds_kept: int

    @classmethod
    def empty(cls) -> "TemplateStats":
        return cls([0, 0, 0], [0.0, 0.0, 0.0], collections.Counter(), 0, 0)

    def add_pattern(self, c0: int, c1: int, c2: int) -> None:
        n = c0 + c1 + c2
        if n <= 0:
            return
        self.sum_counts[0] += c0
        self.sum_counts[1] += c1
        self.sum_counts[2] += c2
        self.sum_props[0] += c0 / n
        self.sum_props[1] += c1 / n
        self.sum_props[2] += c2 / n
        self.patterns[(c0, c1, c2)] += 1
        self.n_codons += 1

    def merge(self, other: "TemplateStats") -> None:
        for i in range(3):
            self.sum_counts[i] += other.sum_counts[i]
            self.sum_props[i] += other.sum_props[i]
        self.patterns.update(other.patterns)
        self.n_codons += other.n_codons
        self.n_cds_kept += other.n_cds_kept


def require_worker_dependencies() -> None:
    """Import the lightweight dependencies required inside worker processes.

    This initializer is required for multiprocessing start methods such as
    ``spawn`` (used by Windows and recent macOS Python builds), where workers
    re-import the module instead of inheriting initialized globals.
    """
    global _np, _pysam
    if _np is None:
        try:
            import numpy as numpy_module
        except ImportError as exc:
            raise SystemExit("ERROR: missing dependency 'numpy'. Install with: conda install numpy") from exc
        _np = numpy_module
    if _pysam is None:
        try:
            import pysam as pysam_module
        except ImportError as exc:
            raise SystemExit("ERROR: missing dependency 'pysam'. Install with: conda install -c bioconda pysam") from exc
        _pysam = pysam_module


def require_dependencies() -> None:
    """Import runtime dependencies after argparse handles --help/--version."""
    global _plt, _gammaln, _minimize_scalar
    require_worker_dependencies()
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as pyplot_module
    except ImportError as exc:
        raise SystemExit("ERROR: missing dependency 'matplotlib'. Install with: conda install matplotlib") from exc
    try:
        from scipy.optimize import minimize_scalar as minimize_scalar_func
        from scipy.special import gammaln as gammaln_func
    except ImportError as exc:
        raise SystemExit("ERROR: missing dependency 'scipy'. Install with: conda install scipy") from exc
    _plt = pyplot_module
    _gammaln = gammaln_func
    _minimize_scalar = minimize_scalar_func


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        prog=PROGRAM,
        description=(
            "MIAO sample-level Ribo-seq metagene QC on annotated MANE CDS records. "
            "Accepts one final tORF TSV, a directory, or a glob."
        )
    )
    ap.add_argument("--version", action="version", version=f"{PROGRAM} {VERSION}")
    ap.add_argument("--psite-bam", required=True, help="Sorted and indexed 1-nt P-site BAM")
    ap.add_argument(
        "--torf",
        required=True,
        help="Final tORF TSV, directory containing *.torf.tsv, or glob such as PREFIX.chr*.torf.tsv",
    )
    ap.add_argument("--out-prefix", required=True, help="Output prefix")
    ap.add_argument("--window", type=int, default=30, help="Start/stop metagene window in nt [default: 30]")
    ap.add_argument("--bins", type=int, default=100, help="Scaled CDS profile bins [default: 100]")
    ap.add_argument("--min-cds-nt", type=int, default=90, help="Minimum MANE CDS length [default: 90]")
    ap.add_argument("--min-reads", type=int, default=1, help="Minimum CDS reads for frame/profile metrics [default: 1]")
    ap.add_argument(
        "--min-read-density",
        type=float,
        default=0.05,
        help="Minimum P-sites per CDS nt for QC inclusion [default: 0.05]",
    )
    ap.add_argument(
        "--uniformity-thr",
        type=float,
        default=1.0 / 3.0,
        help="Per-codon frame-0 fraction threshold for Uniformity [default: 1/3]",
    )
    ap.add_argument("--rl", type=int, default=None, help="Optional RPF read-length filter")
    ap.add_argument(
        "--rl-source",
        choices=["auto", "tag", "query_length"],
        default="auto",
        help="Source for --rl filtering. auto prefers BAM tag RL, then query_length [default: auto]",
    )
    ap.add_argument("--workers", type=int, default=max(1, cpu_count() // 2), help="Worker processes [default: half CPUs]")
    ap.add_argument("--mp-chunksize", type=int, default=1, help="Multiprocessing chunksize [default: 1]")
    ap.add_argument("--template-trim-start-nt", type=int, default=45, help="Trim MANE CDS start for template estimation [default: 45]")
    ap.add_argument("--template-trim-stop-nt", type=int, default=30, help="Trim MANE CDS stop for template estimation [default: 30]")
    ap.add_argument("--template-min-codon-reads", type=int, default=3, help="Minimum reads per CDS codon for template estimation [default: 3]")
    ap.add_argument(
        "--min-a0",
        type=float,
        default=1.0,
        help="Minimum usable DM concentration; lower values fail the background QC gate [default: 1.0]",
    )
    ap.add_argument(
        "--pi-method",
        choices=["codon_equal", "read_weighted"],
        default="codon_equal",
        help="Canonical translation template estimator [default: codon_equal]",
    )
    return ap.parse_args(argv)


def list_input_torfs(spec: str) -> List[str]:
    if any(ch in spec for ch in "*?"):
        files = sorted(glob.glob(spec))
    elif os.path.isdir(spec):
        files = sorted(glob.glob(os.path.join(spec, "*.torf.tsv")))
    elif os.path.isfile(spec):
        files = [spec]
    else:
        files = sorted(glob.glob(f"{spec}.*.torf.tsv"))
    if not files:
        raise SystemExit(f"ERROR: no tORF TSV files found for --torf {spec!r}")
    return files


def parse_int_list_1based_closed(starts_s: str, ends_s: str) -> Tuple[Tuple[int, int], ...]:
    if not starts_s or not ends_s:
        return tuple()
    starts = [int(x) for x in str(starts_s).split(",") if x != ""]
    ends = [int(x) for x in str(ends_s).split(",") if x != ""]
    if len(starts) != len(ends):
        return tuple()
    return tuple((s - 1, e) for s, e in zip(starts, ends) if e >= s)


def parse_blocks_01(block_sizes: str, starts1: str) -> Tuple[Tuple[int, int], ...]:
    sizes = [int(x) for x in str(block_sizes).split(",") if x != ""]
    starts = [int(x) for x in str(starts1).split(",") if x != ""]
    if not sizes or len(sizes) != len(starts):
        return tuple()
    return tuple((s1 - 1, s1 - 1 + size) for s1, size in zip(starts, sizes) if size > 0)


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def load_mane_cds(paths: Sequence[str], min_cds_nt: int) -> Tuple[List[CDSRecord], Dict[str, int]]:
    required = {
        "orf_biotype", "is_annotated_cds", "has_MANE_Select_tag", "chrom", "strand",
        "t_start", "t_end", "chromStart1", "chromEnd1", "blockSizes", "genomic_block_starts1",
        "flank5_genomic_starts1", "flank5_genomic_ends1", "flank3_genomic_starts1", "flank3_genomic_ends1",
    }
    records: List[CDSRecord] = []
    seen: set[str] = set()
    stats = collections.Counter()
    fallback_idx = 0
    for path in paths:
        with open(path, "r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = sorted(required - set(reader.fieldnames or []))
            if missing:
                raise SystemExit(f"ERROR: columns missing in {path}: {', '.join(missing)}")
            for row in reader:
                stats["torf_rows_scanned"] += 1
                if not (row.get("orf_biotype") == "CDS" and is_true(row.get("is_annotated_cds")) and is_true(row.get("has_MANE_Select_tag"))):
                    continue
                stats["mane_cds_rows"] += 1
                try:
                    t_start = int(row["t_start"])
                    t_end = int(row["t_end"])
                    chrom_start1 = int(row["chromStart1"])
                    chrom_end1 = int(row["chromEnd1"])
                except (TypeError, ValueError):
                    stats["skip_invalid_coordinates"] += 1
                    continue
                if t_end - t_start < min_cds_nt:
                    stats["skip_short_cds"] += 1
                    continue
                blocks = parse_blocks_01(row.get("blockSizes", ""), row.get("genomic_block_starts1", ""))
                if not blocks:
                    stats["skip_invalid_blocks"] += 1
                    continue
                strand = str(row.get("strand", "+") or "+").strip()
                if strand not in {"+", "-"}:
                    stats["skip_invalid_strand"] += 1
                    continue
                torf_id = str(row.get("torf_id", "")).strip()
                if not torf_id:
                    fallback_idx += 1
                    torf_id = f"CDS_fallback_{fallback_idx}"
                if torf_id in seen:
                    stats["skip_duplicate_torf_id"] += 1
                    continue
                seen.add(torf_id)
                orf_len_nt = sum(end - start for start, end in blocks)
                if orf_len_nt < 3:
                    stats["skip_invalid_orf_length"] += 1
                    continue
                start_anchor = gpos_at_orf_tindex(blocks, strand, 0)
                stop_anchor = gpos_at_orf_tindex(blocks, strand, orf_len_nt - 3)
                records.append(
                    CDSRecord(
                        torf_id=torf_id,
                        gene_id=str(row.get("gene_id", "")),
                        gene_name=str(row.get("gene_name", "")),
                        transcript_id=str(row.get("transcript_id", "")),
                        chrom=str(row.get("chrom", "")),
                        strand=strand,
                        blocks=blocks,
                        flank5=parse_int_list_1based_closed(row.get("flank5_genomic_starts1", ""), row.get("flank5_genomic_ends1", "")),
                        flank3=parse_int_list_1based_closed(row.get("flank3_genomic_starts1", ""), row.get("flank3_genomic_ends1", "")),
                        start_anchor=start_anchor,
                        stop_anchor=stop_anchor,
                    )
                )
    stats["mane_cds_loaded"] = len(records)
    return records, dict(stats)


def bam_strand(read: Any) -> str:
    return "-" if read.is_reverse else "+"


def read_has_rl(read: Any, rl: Optional[int], rl_source: str) -> bool:
    if rl is None:
        return True
    if rl_source == "tag":
        try:
            return int(read.get_tag("RL")) == rl
        except KeyError:
            return False
    if rl_source == "query_length":
        return int(read.query_length or 0) == rl
    try:
        return int(read.get_tag("RL")) == rl
    except KeyError:
        return int(read.query_length or 0) == rl


def validate_bam_and_rl(bam_path: str, rl: Optional[int], rl_source: str, max_reads: int = 5000) -> None:
    if not os.path.isfile(bam_path):
        raise SystemExit(f"ERROR: BAM not found: {bam_path}")
    bai_candidates = [f"{bam_path}.bai", os.path.splitext(bam_path)[0] + ".bai"]
    if not any(os.path.isfile(path) for path in bai_candidates):
        raise SystemExit(f"ERROR: BAM index not found. Expected one of: {', '.join(bai_candidates)}")
    bam = _pysam.AlignmentFile(bam_path, "rb")
    try:
        try:
            bam.check_index()
        except Exception as exc:
            raise SystemExit(f"ERROR: BAM index is not readable: {exc}") from exc
        if rl is None:
            return
        seen = 0
        tag_count = 0
        tag_match = 0
        qlens = collections.Counter()
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            seen += 1
            qlen = int(read.query_length or 0)
            qlens[qlen] += 1
            try:
                value = int(read.get_tag("RL"))
                tag_count += 1
                if value == rl:
                    tag_match += 1
            except KeyError:
                pass
            if seen >= max_reads:
                break
        if seen == 0:
            raise SystemExit("ERROR: BAM contains no mapped alignments")
        qlen_match = qlens.get(int(rl), 0)
        summary = ",".join(f"{k}:{v}" for k, v in qlens.most_common(8))
        sys.stderr.write(
            f"[rl] source={rl_source}, requested={rl}, sampled_reads={seen:,}, "
            f"RL_tags={tag_count:,}, RL_matches={tag_match:,}, query_length_matches={qlen_match:,}, query_lengths={summary}\n"
        )
        if rl_source == "tag" and tag_count == 0:
            raise SystemExit("ERROR: --rl-source tag requested, but sampled reads have no RL tag")
        if rl_source == "query_length" and qlen_match == 0:
            raise SystemExit("ERROR: --rl-source query_length requested, but sampled reads do not match --rl")
        if rl_source == "auto" and tag_match == 0 and qlen_match == 0:
            if tag_count == 0 and set(qlens).issubset({0, 1}) and rl != 1:
                raise SystemExit(
                    "ERROR: original RPF length cannot be recovered from this apparent 1-nt P-site BAM: "
                    "reads lack RL tags and query_length is 1. Omit --rl or recreate the P-site BAM with RL tags."
                )
            raise SystemExit("ERROR: no sampled BAM reads match --rl under --rl-source auto")
    finally:
        bam.close()


def merge_pieces(pieces: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    items = sorted((s, e) for s, e in pieces if e > s)
    if not items:
        return []
    out = [list(items[0])]
    for s, e in items[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def take_first_w(pieces: Sequence[Tuple[int, int]], w: int, strand: str) -> List[Tuple[int, int]]:
    if w <= 0:
        return []
    out: List[Tuple[int, int]] = []
    need = w
    iterator = iter(pieces) if strand == "+" else reversed(pieces)
    for s, e in iterator:
        take = min(need, e - s)
        if take <= 0:
            continue
        out.append((s, s + take) if strand == "+" else (e - take, e))
        need -= take
        if need <= 0:
            break
    return merge_pieces(out)


def take_last_w(pieces: Sequence[Tuple[int, int]], w: int, strand: str) -> List[Tuple[int, int]]:
    if w <= 0:
        return []
    out: List[Tuple[int, int]] = []
    need = w
    iterator = reversed(pieces) if strand == "+" else iter(pieces)
    for s, e in iterator:
        take = min(need, e - s)
        if take <= 0:
            continue
        out.append((e - take, e) if strand == "+" else (s, s + take))
        need -= take
        if need <= 0:
            break
    return merge_pieces(out)


def clip_transcript_end(pieces: Sequence[Tuple[int, int]], nt: int, strand: str) -> List[Tuple[int, int]]:
    out = [list(x) for x in merge_pieces(pieces)]
    remaining = max(0, int(nt))
    if strand == "+":
        while remaining > 0 and out:
            seg = out[-1][1] - out[-1][0]
            if seg <= remaining:
                remaining -= seg
                out.pop()
            else:
                out[-1][1] -= remaining
                remaining = 0
    else:
        while remaining > 0 and out:
            seg = out[0][1] - out[0][0]
            if seg <= remaining:
                remaining -= seg
                out.pop(0)
            else:
                out[0][0] += remaining
                remaining = 0
    return [(s, e) for s, e in out if e > s]


def gpos_at_orf_tindex(
    blocks: Sequence[Tuple[int, int]],
    strand: str,
    tidx: int,
) -> int:
    """Map a 0-based spliced ORF-relative nucleotide index to genomic position.

    ``blocks`` use 0-based half-open genomic coordinates. The returned genomic
    coordinate is also 0-based. This helper correctly handles negative-strand
    ORFs and the rare case where an anchor codon spans exon junctions.
    """
    if tidx < 0:
        raise ValueError(f"Negative ORF-relative index: {tidx}")
    ordered = blocks if strand == "+" else tuple(reversed(blocks))
    remaining = tidx
    for start, end in ordered:
        seg_len = end - start
        if remaining < seg_len:
            return start + remaining if strand == "+" else end - 1 - remaining
        remaining -= seg_len
    raise ValueError(
        f"ORF-relative index {tidx} exceeds ORF length "
        f"{sum(end - start for start, end in blocks)}"
    )


def tindex_of_gpos(rec: CDSRecord, gpos: int) -> Optional[int]:
    offset = 0
    blocks = rec.blocks if rec.strand == "+" else tuple(reversed(rec.blocks))
    for s, e in blocks:
        if s <= gpos < e:
            return offset + ((gpos - s) if rec.strand == "+" else (e - 1 - gpos))
        offset += e - s
    return None


def spliced_index_of_gpos(
    pieces: Sequence[Tuple[int, int]],
    strand: str,
    gpos: int,
) -> Optional[int]:
    """Return the 0-based transcript-order index of ``gpos`` in ``pieces``."""
    offset = 0
    ordered = pieces if strand == "+" else tuple(reversed(pieces))
    for start, end in ordered:
        if start <= gpos < end:
            return offset + ((gpos - start) if strand == "+" else (end - 1 - gpos))
        offset += end - start
    return None


def metagene_relative_index(rec: CDSRecord, gpos: int, which: str) -> Optional[int]:
    """Map a genomic P-site to a spliced start/stop-relative coordinate."""
    cds_index = tindex_of_gpos(rec, gpos)
    if cds_index is not None:
        anchor_index = 0 if which == "start" else rec.orf_len_nt - 3
        return cds_index - anchor_index

    if which == "start":
        flank_index = spliced_index_of_gpos(rec.flank5, rec.strand, gpos)
        if flank_index is not None:
            flank_len = sum(end - start for start, end in rec.flank5)
            return flank_index - flank_len
    else:
        flank_index = spliced_index_of_gpos(rec.flank3, rec.strand, gpos)
        if flank_index is not None:
            # The stop anchor is the first nucleotide of the stop codon, so the
            # first 3' flank nucleotide has relative position +3.
            return flank_index + 3
    return None


def count_tindices(bam: Any, rec: CDSRecord, rl: Optional[int], rl_source: str) -> List[int]:
    out: List[int] = []
    for s, e in rec.blocks:
        for read in bam.fetch(rec.chrom, s, e):
            if read.is_unmapped or bam_strand(read) != rec.strand or not read_has_rl(read, rl, rl_source):
                continue
            pos = read.reference_start
            if not (s <= pos < e):
                continue
            tidx = tindex_of_gpos(rec, pos)
            if tidx is not None and 0 <= tidx < rec.orf_len_nt:
                out.append(tidx)
    return out


def count_piece_reads(bam: Any, rec: CDSRecord, pieces: Sequence[Tuple[int, int]], rl: Optional[int], rl_source: str) -> int:
    total = 0
    for s, e in pieces:
        for read in bam.fetch(rec.chrom, s, e):
            if read.is_unmapped or bam_strand(read) != rec.strand or not read_has_rl(read, rl, rl_source):
                continue
            if s <= read.reference_start < e:
                total += 1
    return total


def relative_hist(bam: Any, rec: CDSRecord, which: str, window: int, rl: Optional[int], rl_source: str) -> collections.Counter:
    out = collections.Counter()
    if which == "start":
        before = take_last_w(rec.flank5, window, rec.strand)
        after = take_first_w(rec.blocks, window + 1, rec.strand)
    else:
        # Include -window through the complete three-nucleotide stop codon.
        before = take_last_w(rec.blocks, window + 3, rec.strand)
        # Relative +1/+2 are in the stop codon; the 3' flank begins at +3.
        after = take_first_w(rec.flank3, max(0, window - 2), rec.strand)
    for pieces in (before, after):
        for s, e in pieces:
            for read in bam.fetch(rec.chrom, s, e):
                if read.is_unmapped or bam_strand(read) != rec.strand or not read_has_rl(read, rl, rl_source):
                    continue
                pos = read.reference_start
                if not (s <= pos < e):
                    continue
                rel = metagene_relative_index(rec, pos, which)
                if rel is not None and -window <= rel <= window:
                    out[rel] += 1
    return out


def codon_patterns_from_tindices(tindices: Sequence[int], orf_len_nt: int, trim_start_nt: int, trim_stop_nt: int, min_codon_reads: int) -> List[Tuple[int, int, int]]:
    n_codons = orf_len_nt // 3
    matrix = [[0, 0, 0] for _ in range(n_codons)]
    for tidx in tindices:
        codon = tidx // 3
        frame = tidx % 3
        if 0 <= codon < n_codons:
            matrix[codon][frame] += 1
    first = int(math.ceil(trim_start_nt / 3.0))
    last = int(math.floor((orf_len_nt - trim_stop_nt) / 3.0))
    if last <= first:
        return []
    return [tuple(row) for row in matrix[first:last] if sum(row) >= min_codon_reads]


def cds_metrics(rec: CDSRecord, bam: Any, tindices: Sequence[int], args: argparse.Namespace) -> Dict[str, object]:
    total = len(tindices)
    codon_count = max(0, rec.orf_len_nt // 3)
    frame_counts = [0, 0, 0]
    codon_hits = [0] * codon_count
    codon_frame0 = [0] * codon_count
    for tidx in tindices:
        frame = tidx % 3
        frame_counts[frame] += 1
        codon = tidx // 3
        if 0 <= codon < codon_count:
            codon_hits[codon] += 1
            if frame == 0:
                codon_frame0[codon] += 1
    pif = frame_counts[0] / total if total else 0.0
    codons_with_reads = sum(x > 0 for x in codon_hits)
    codons_passing = sum(h > 0 and (codon_frame0[i] / h) > args.uniformity_thr for i, h in enumerate(codon_hits))
    uniformity = codons_passing / codon_count if codon_count else 0.0

    w = 15
    start_before = count_piece_reads(bam, rec, take_last_w(rec.flank5, w, rec.strand), args.rl, args.rl_source)
    start_after = count_piece_reads(bam, rec, take_first_w(rec.blocks, w, rec.strand), args.rl, args.rl_source)
    start_increase = start_after / (start_before + start_after) if (start_before + start_after) else 0.0

    coding_without_stop = clip_transcript_end(rec.blocks, 3, rec.strand)
    stop_before = count_piece_reads(bam, rec, take_last_w(coding_without_stop, w, rec.strand), args.rl, args.rl_source)
    stop_after = count_piece_reads(bam, rec, take_first_w(rec.flank3, w, rec.strand), args.rl, args.rl_source)
    dropoff = stop_before / (stop_before + stop_after) if (stop_before + stop_after) else 0.0

    tis_reads = sum(tidx == 0 for tidx in tindices)
    mean_frame0 = frame_counts[0] / codon_count if codon_count else 0.0
    tis_global = tis_reads / mean_frame0 if mean_frame0 > 0 else 0.0
    local_bg = [sum(tidx == pos for tidx in tindices) for pos in range(3, min(rec.orf_len_nt, 33), 3)]
    local_median = float(_np.median(local_bg)) if local_bg else 0.0
    tis_local = (tis_reads + 1.0) / (local_median + 1.0)

    return {
        "torf_id": rec.torf_id,
        "gene_id": rec.gene_id,
        "gene_name": rec.gene_name,
        "transcript_id": rec.transcript_id,
        "chrom": rec.chrom,
        "strand": rec.strand,
        "cds_nt": rec.orf_len_nt,
        "total_reads": total,
        "read_density": total / rec.orf_len_nt if rec.orf_len_nt else 0.0,
        "PIF": pif,
        "Uniformity": uniformity,
        "StartIncrease": start_increase,
        "DropOff": dropoff,
        "TIS_global_enrich": tis_global,
        "TIS_local_score": tis_local,
        "codon_count": codon_count,
        "codons_with_reads": codons_with_reads,
        "start_before": start_before,
        "start_after": start_after,
        "stop_before": stop_before,
        "stop_after": stop_after,
    }


def process_chrom(task: Tuple[str, List[CDSRecord], str, argparse.Namespace]) -> Dict[str, object]:
    require_worker_dependencies()
    chrom, records, bam_path, args = task
    bam = _pysam.AlignmentFile(bam_path, "rb")
    start_hist = collections.Counter()
    stop_hist = collections.Counter()
    frame_counts = [0, 0, 0]
    profile_sum = [0.0] * args.bins
    profile_frame_sum = [[0.0] * args.bins for _ in range(3)]
    profile_cds_n = 0
    metrics_rows: List[Dict[str, object]] = []
    template = TemplateStats.empty()
    kept = 0
    try:
        for rec in records:
            tindices = count_tindices(bam, rec, args.rl, args.rl_source)
            density = len(tindices) / float(max(rec.orf_len_nt, 1))
            if density < args.min_read_density:
                continue
            kept += 1
            start_hist.update(relative_hist(bam, rec, "start", args.window, args.rl, args.rl_source))
            stop_hist.update(relative_hist(bam, rec, "stop", args.window, args.rl, args.rl_source))
            metrics_rows.append(cds_metrics(rec, bam, tindices, args))

            if len(tindices) >= args.min_reads:
                profile_cds_n += 1
                bins_total = [0] * args.bins
                bins_frame = [[0] * args.bins for _ in range(3)]
                for tidx in tindices:
                    frame = tidx % 3
                    frame_counts[frame] += 1
                    b = min(args.bins - 1, max(0, int(math.floor(tidx * args.bins / rec.orf_len_nt))))
                    bins_total[b] += 1
                    bins_frame[frame][b] += 1
                for b in range(args.bins):
                    profile_sum[b] += bins_total[b]
                    for frame in range(3):
                        profile_frame_sum[frame][b] += bins_frame[frame][b]

            patterns = codon_patterns_from_tindices(
                tindices,
                rec.orf_len_nt,
                args.template_trim_start_nt,
                args.template_trim_stop_nt,
                args.template_min_codon_reads,
            )
            if patterns:
                template.n_cds_kept += 1
                for c0, c1, c2 in patterns:
                    template.add_pattern(c0, c1, c2)
    finally:
        bam.close()
    sys.stderr.write(f"[worker] chrom={chrom} CDS={len(records):,}, kept={kept:,}\n")
    return {
        "chrom": chrom,
        "records": len(records),
        "kept": kept,
        "start_hist": start_hist,
        "stop_hist": stop_hist,
        "frame_counts": frame_counts,
        "profile_sum": profile_sum,
        "profile_frame_sum": profile_frame_sum,
        "profile_cds_n": profile_cds_n,
        "metrics": metrics_rows,
        "template": template,
    }


def dm_logpmf_pattern(counts: Sequence[int], alpha: Sequence[float]) -> float:
    n = int(sum(counts))
    a_sum = float(sum(alpha))
    return float(
        _gammaln(n + 1)
        - sum(_gammaln(c + 1) for c in counts)
        + _gammaln(a_sum)
        - _gammaln(n + a_sum)
        + sum(_gammaln(c + a) - _gammaln(a) for c, a in zip(counts, alpha))
    )


def estimate_a0_mle(patterns: Mapping[Tuple[int, int, int], int], pi: Sequence[float]) -> float:
    pi_arr = _np.asarray(pi, dtype=float)
    pi_arr = _np.clip(pi_arr, 1e-9, 1.0)
    pi_arr /= pi_arr.sum()
    unique = [(key, int(freq)) for key, freq in patterns.items()]
    if not unique:
        return 20.0

    def objective(log_a: float) -> float:
        a0 = float(math.exp(log_a))
        alpha = a0 * pi_arr
        return -sum(freq * dm_logpmf_pattern(key, alpha) for key, freq in unique)

    fit = _minimize_scalar(objective, bounds=(math.log(0.02), math.log(1e5)), method="bounded", options={"xatol": 1e-5, "maxiter": 300})
    if not fit.success:
        sys.stderr.write(f"WARN: A0 optimization did not fully converge: {fit.message}\n")
    return float(math.exp(fit.x))


def estimate_template(stats: TemplateStats, method: str) -> Dict[str, object]:
    if stats.n_codons <= 0 or sum(stats.sum_counts) <= 0:
        raise ValueError("no template codons survived filters")
    pi_read = _np.asarray(stats.sum_counts, dtype=float)
    pi_read /= pi_read.sum()
    pi_equal = _np.asarray(stats.sum_props, dtype=float) / float(stats.n_codons)
    pi_equal /= pi_equal.sum()
    pi_used = pi_read if method == "read_weighted" else pi_equal
    a0 = estimate_a0_mle(stats.patterns, pi_used)
    return {
        "pi_read_weighted": pi_read,
        "pi_codon_equal": pi_equal,
        "pi_used": pi_used,
        "pi_method": method,
        "A0": a0,
        "alpha0": a0 * pi_used,
        "template_codons": stats.n_codons,
        "template_unique_patterns": len(stats.patterns),
        "template_cds_kept": stats.n_cds_kept,
    }


def write_key_value(path: str, rows: Sequence[Tuple[str, object]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in rows:
            if isinstance(value, float):
                value = f"{value:.10g}"
            handle.write(f"{key}\t{value}\n")


def write_counts(path: str, counter: Mapping[int, int], window: int) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("rel_nt\tcount\n")
        for pos in range(-window, window + 1):
            handle.write(f"{pos}\t{counter.get(pos, 0)}\n")


def codon_ternary_grid_key(
    counts: Sequence[int],
    resolution: int = CODON_TERNARY_RESOLUTION,
) -> Tuple[int, int, int]:
    """Map one observed codon composition to a deterministic ternary lattice."""
    if resolution <= 0:
        raise ValueError("ternary resolution must be positive")
    if len(counts) != 3:
        raise ValueError("codon composition must contain exactly three frame counts")
    values = [int(value) for value in counts]
    total = sum(values)
    if total <= 0 or min(values) < 0:
        raise ValueError("codon frame counts must be non-negative with a positive total")
    scaled = [value * resolution / float(total) for value in values]
    base = [int(math.floor(value)) for value in scaled]
    missing = resolution - sum(base)
    if missing > 0:
        order = sorted(
            range(3),
            key=lambda idx: (-(scaled[idx] - base[idx]), idx),
        )
        for idx in order[:missing]:
            base[idx] += 1
    key = tuple(base)
    if sum(key) != resolution or min(key) < 0:
        raise RuntimeError(f"invalid ternary lattice key: {key}")
    return key  # type: ignore[return-value]


def summarize_codon_ternary(
    stats: TemplateStats,
    resolution: int = CODON_TERNARY_RESOLUTION,
) -> List[Dict[str, object]]:
    """Aggregate A0-input codons on a fixed frame-composition grid."""
    grid = collections.Counter()
    for counts, frequency in stats.patterns.items():
        grid[codon_ternary_grid_key(counts, resolution)] += int(frequency)
    total = int(sum(grid.values()))
    rows: List[Dict[str, object]] = []
    for key, count in grid.items():
        p0, p1, p2 = (value / float(resolution) for value in key)
        rows.append({
            "grid0": int(key[0]),
            "grid1": int(key[1]),
            "grid2": int(key[2]),
            "p0": p0,
            "p1": p1,
            "p2": p2,
            "codon_count": int(count),
            "total_codons": total,
            "sample_percentage": 100.0 * count / total if total else 0.0,
        })
    rows.sort(
        key=lambda row: (
            -float(row["sample_percentage"]),
            -int(row["grid0"]),
            -int(row["grid1"]),
            -int(row["grid2"]),
        )
    )
    return rows


def write_codon_ternary(path: str, rows: Sequence[Mapping[str, object]]) -> None:
    columns = [
        "rank", "grid_percent", "grid_P0_percent", "grid_P1_percent",
        "grid_P2_percent", "codon_count", "total_codons", "sample_percentage",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            writer.writerow({
                "rank": rank,
                "grid_percent": CODON_TERNARY_GRID_PERCENT,
                "grid_P0_percent": 100.0 * float(row["p0"]),
                "grid_P1_percent": 100.0 * float(row["p1"]),
                "grid_P2_percent": 100.0 * float(row["p2"]),
                "codon_count": int(row["codon_count"]),
                "total_codons": int(row["total_codons"]),
                "sample_percentage": float(row["sample_percentage"]),
            })


def codon_ternary_vertex_percentages(
    rows: Sequence[Mapping[str, object]],
    resolution: int = CODON_TERNARY_RESOLUTION,
) -> Tuple[float, float, float]:
    by_key = {
        (int(row["grid0"]), int(row["grid1"]), int(row["grid2"])):
        float(row["sample_percentage"])
        for row in rows
    }
    return (
        by_key.get((resolution, 0, 0), 0.0),
        by_key.get((0, resolution, 0), 0.0),
        by_key.get((0, 0, resolution), 0.0),
    )


def codon_ternary_mean_percentages(stats: TemplateStats) -> Tuple[float, float, float]:
    """Return the exact codon-equal mean P0/P1/P2 composition in percent."""
    weighted = [0.0, 0.0, 0.0]
    total_codons = 0
    for counts, frequency in stats.patterns.items():
        depth = sum(counts)
        if depth <= 0 or frequency <= 0:
            continue
        for phase in range(3):
            weighted[phase] += frequency * counts[phase] / depth
        total_codons += frequency
    if total_codons <= 0:
        return 0.0, 0.0, 0.0
    return tuple(100.0 * value / total_codons for value in weighted)  # type: ignore[return-value]


def ternary_xy(p0: float, p1: float, p2: float) -> Tuple[float, float]:
    del p1
    return 0.5 * p0 + p2, math.sqrt(3.0) * p0 / 2.0


def plot_codon_ternary(
    prefix: str,
    rows: Sequence[Mapping[str, object]],
    a0: Optional[float],
    template_status: str,
    mean_percentages: Optional[Tuple[float, float, float]] = None,
) -> None:
    """Write the default codon-level QC as a normalized linear-area bubble plot."""
    if not rows:
        return
    sample_name = os.path.basename(os.path.normpath(prefix)) or "sample"
    total_codons = int(rows[0]["total_codons"])
    vertex = codon_ternary_vertex_percentages(rows)
    if mean_percentages is None:
        total_share = sum(float(row["sample_percentage"]) for row in rows)
        if total_share > 0:
            mean_percentages = tuple(
                100.0 * sum(
                    float(row[f"p{phase}"]) * float(row["sample_percentage"])
                    for row in rows
                ) / total_share
                for phase in range(3)
            )  # type: ignore[assignment]
        else:
            mean_percentages = (0.0, 0.0, 0.0)
    height = math.sqrt(3.0) / 2.0
    fig = _plt.figure(figsize=(7.2, 7.4), dpi=160)
    ax = fig.add_axes([0.10, 0.245, 0.80, 0.60])

    triangle = [(0.0, 0.0), (1.0, 0.0), (0.5, height), (0.0, 0.0)]
    ax.plot(
        [point[0] for point in triangle],
        [point[1] for point in triangle],
        color="#253044",
        linewidth=1.35,
        zorder=1,
    )
    for fraction in (0.25, 0.50, 0.75):
        line_sets = [
            (ternary_xy(fraction, 1-fraction, 0), ternary_xy(fraction, 0, 1-fraction)),
            (ternary_xy(1-fraction, fraction, 0), ternary_xy(0, fraction, 1-fraction)),
            (ternary_xy(1-fraction, 0, fraction), ternary_xy(0, 1-fraction, fraction)),
        ]
        for left, right in line_sets:
            ax.plot(
                [left[0], right[0]],
                [left[1], right[1]],
                color="#d8dee8",
                linewidth=0.6,
                zorder=0,
            )

    ordered = list(reversed(rows))
    x_values = []
    y_values = []
    sizes = []
    for row in ordered:
        x_value, y_value = ternary_xy(float(row["p0"]), float(row["p1"]), float(row["p2"]))
        x_values.append(x_value)
        y_values.append(y_value)
        sizes.append(max(0.25, CODON_TERNARY_AREA_PER_PERCENT * float(row["sample_percentage"])))
    ax.scatter(
        x_values,
        y_values,
        s=sizes,
        color="#258b85",
        alpha=0.82,
        edgecolors="none",
        rasterized=True,
        zorder=2,
    )
    emphasized = [row for row in rows if float(row["sample_percentage"]) >= 0.5]
    if emphasized:
        ex, ey, es = [], [], []
        for row in emphasized:
            x_value, y_value = ternary_xy(float(row["p0"]), float(row["p1"]), float(row["p2"]))
            ex.append(x_value)
            ey.append(y_value)
            es.append(CODON_TERNARY_AREA_PER_PERCENT * float(row["sample_percentage"]))
        ax.scatter(
            ex,
            ey,
            s=es,
            color="#258b85",
            alpha=0.92,
            edgecolors="white",
            linewidths=0.55,
            rasterized=True,
            zorder=3,
        )

    mean_x, mean_y = ternary_xy(
        mean_percentages[0] / 100.0,
        mean_percentages[1] / 100.0,
        mean_percentages[2] / 100.0,
    )
    ax.scatter(
        [mean_x],
        [mean_y],
        marker="*",
        s=135,
        facecolor="white",
        edgecolor="#111827",
        linewidth=1.1,
        zorder=8,
    )

    ax.text(0.5, height + 0.092, "P0", ha="center", va="bottom", fontsize=11, weight="bold")
    ax.text(-0.035, -0.022, "P1", ha="right", va="top", fontsize=11, weight="bold")
    ax.text(1.035, -0.022, "P2", ha="left", va="top", fontsize=11, weight="bold")
    ax.set_xlim(-0.10, 1.10)
    ax.set_ylim(-0.075, height + 0.125)
    ax.set_aspect("equal")
    ax.axis("off")

    a0_label = "NA" if a0 is None else f"{a0:.3f}"
    status_label = template_status.replace("_", " ")
    fig.suptitle("Codon-level phase composition", fontsize=15, weight="bold", y=0.975)
    fig.text(
        0.5,
        0.929,
        sample_name,
        ha="center",
        va="top",
        fontsize=11,
        weight="semibold",
        color="#253044",
    )
    fig.text(
        0.5,
        0.897,
        f"N = {total_codons:,} template codons   ·   A0 = {a0_label}   ·   {status_label}",
        ha="center",
        va="top",
        fontsize=9.5,
        color="#4b5563",
    )
    fig.text(
        0.5,
        0.866,
        "★  Codon-equal mean (%)   "
        f"P0 {mean_percentages[0]:.1f}   ·   "
        f"P1 {mean_percentages[1]:.1f}   ·   "
        f"P2 {mean_percentages[2]:.1f}",
        ha="center",
        va="top",
        fontsize=9.5,
        color="#253044",
    )
    fig.text(
        0.5,
        0.186,
        f"5% grid · exact-vertex share (%)   P0 {vertex[0]:.1f}   ·   P1 {vertex[1]:.1f}   ·   P2 {vertex[2]:.1f}",
        ha="center",
        va="center",
        fontsize=9.5,
        color="#253044",
    )
    fig.text(
        0.5,
        0.115,
        "Linear bubble area = share of template codons",
        ha="center",
        va="center",
        fontsize=9.5,
        color="#253044",
    )
    legend_values = (1.0, 10.0, 50.0)
    legend_handles = [
        ax.scatter(
            [], [],
            s=CODON_TERNARY_AREA_PER_PERCENT * value,
            color="#258b85",
            alpha=0.85,
            edgecolors="white",
            linewidths=0.5,
        )
        for value in legend_values
    ]
    fig.legend(
        legend_handles,
        [f"{value:g}%" for value in legend_values],
        ncol=3,
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.025),
        fontsize=9,
        columnspacing=2.4,
        handletextpad=0.8,
    )
    fig.savefig(
        f"{prefix}.codon_frame_ternary.png",
        dpi=240,
        bbox_inches="tight",
        pad_inches=0.12,
        facecolor="white",
    )
    fig.savefig(
        f"{prefix}.codon_frame_ternary.pdf",
        bbox_inches="tight",
        pad_inches=0.12,
        facecolor="white",
    )
    _plt.close(fig)


def clear_previous_outputs(prefix: str) -> None:
    """Remove known outputs so a failed run cannot be mixed with older files."""
    suffixes = (
        ".start.meta.tsv", ".stop.meta.tsv", ".frame.tsv", ".cds_profile.tsv",
        ".cds_profile.f0.tsv", ".cds_profile.f1.tsv", ".cds_profile.f2.tsv",
        ".cds_metrics.tsv", ".template.tsv", ".dm_background.tsv",
        ".qc_summary.tsv", ".metagene.png", ".metagene.pdf",
        ".codon_frame_ternary.tsv", ".codon_frame_ternary.png",
        ".codon_frame_ternary.pdf",
    )
    for suffix in suffixes:
        path = f"{prefix}{suffix}"
        try:
            os.remove(path)
        except FileNotFoundError:
            pass


def write_frame(path: str, counts: Sequence[int]) -> None:
    total = sum(counts)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("frame\tcount\tproportion\ttotal\n")
        for frame, count in enumerate(counts):
            prop = count / total if total else 0.0
            handle.write(f"{frame}\t{count}\t{prop:.10g}\t{total}\n")


def write_profile(path: str, values: Sequence[float]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("bin\tmean_count_per_cds\n")
        for idx, value in enumerate(values):
            handle.write(f"{idx}\t{value:.10g}\n")


def write_metrics(path: str, rows: Sequence[Mapping[str, object]]) -> None:
    columns = [
        "torf_id", "gene_id", "gene_name", "transcript_id", "chrom", "strand", "cds_nt", "total_reads", "read_density",
        "PIF", "Uniformity", "StartIncrease", "DropOff", "TIS_global_enrich", "TIS_local_score", "codon_count",
        "codons_with_reads", "start_before", "start_after", "stop_before", "stop_after",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_metagene(prefix: str, start: Mapping[int, int], stop: Mapping[int, int], frames: Sequence[int], profile_frames: Sequence[Sequence[float]], window: int) -> None:
    colors = ["#ff66b3", "#33cc66", "#3399ff"]
    total = sum(frames)
    props = [x / total if total else 0.0 for x in frames]
    x_start = list(range(-window, window + 1))
    y_start = [start.get(x, 0) for x in x_start]
    x_stop = list(range(-window, window + 1))
    y_stop = [stop.get(x, 0) for x in x_stop]
    fig = _plt.figure(figsize=(14, 3.2), dpi=150)
    ax1 = fig.add_subplot(1, 4, 1)
    ax1.bar([0, 1, 2], props, color=colors)
    ax1.set_xticks([0, 1, 2])
    ax1.set_ylim(0, 1)
    ax1.set_title("Frame distribution")
    ax1.set_xlabel("Frame")
    for idx, value in enumerate(props):
        ax1.text(idx, value + 0.02, f"{value:.2f}", ha="center", va="bottom", fontsize=8)
    ax2 = fig.add_subplot(1, 4, 2)
    ax2.bar(x_start, y_start, width=1.0, color=[colors[x % 3] for x in x_start])
    ax2.axvline(0, linestyle="--", linewidth=1)
    ax2.set_xlim(-window, window)
    ax2.set_title("P-sites near CDS start")
    ax2.set_xlabel("Relative nt")
    ax3 = fig.add_subplot(1, 4, 3)
    ax3.bar(x_stop, y_stop, width=1.0, color=[colors[x % 3] for x in x_stop])
    ax3.axvline(0, linestyle="--", linewidth=1)
    ax3.set_xlim(-window, window)
    ax3.set_title("P-sites near CDS stop")
    ax3.set_xlabel("Relative nt")
    ax4 = fig.add_subplot(1, 4, 4)
    for frame in range(3):
        ax4.plot(range(len(profile_frames[frame])), profile_frames[frame], color=colors[frame], linewidth=1.8, label=f"frame{frame}")
    ax4.set_xlim(0, max(0, len(profile_frames[0]) - 1))
    ax4.set_title("Scaled CDS profile")
    ax4.set_xlabel("Scaled CDS bin")
    ax4.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{prefix}.metagene.png", bbox_inches="tight")
    fig.savefig(f"{prefix}.metagene.pdf", bbox_inches="tight")
    _plt.close(fig)


def summarize_metric(rows: Sequence[Mapping[str, object]], key: str) -> Tuple[float, float, int]:
    values = _np.asarray([float(row[key]) for row in rows], dtype=float)
    if values.size == 0:
        return 0.0, 0.0, 0
    return float(_np.nanmean(values)), float(_np.nanmedian(values)), int(values.size)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    require_dependencies()
    if args.min_a0 <= 0:
        raise SystemExit("ERROR: --min-a0 must be > 0")
    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)), exist_ok=True)
    clear_previous_outputs(args.out_prefix)
    validate_bam_and_rl(args.psite_bam, args.rl, args.rl_source)
    tsvs = list_input_torfs(args.torf)
    sys.stderr.write(f"[init] input tORF files: {len(tsvs):,}\n")
    records, load_stats = load_mane_cds(tsvs, args.min_cds_nt)
    if not records:
        raise SystemExit("ERROR: no annotated MANE CDS records loaded")
    sys.stderr.write(f"[init] MANE CDS loaded: {len(records):,}\n")
    by_chrom: Dict[str, List[CDSRecord]] = collections.defaultdict(list)
    for rec in records:
        by_chrom[rec.chrom].append(rec)
    tasks = [(chrom, rows, args.psite_bam, args) for chrom, rows in sorted(by_chrom.items())]
    workers = min(max(1, args.workers), len(tasks))
    sys.stderr.write(f"[mp] chromosome tasks={len(tasks):,}, workers={workers}\n")
    if workers <= 1:
        results = [process_chrom(task) for task in tasks]
    else:
        with Pool(processes=workers, initializer=require_worker_dependencies) as pool:
            results = list(pool.imap_unordered(process_chrom, tasks, chunksize=max(1, args.mp_chunksize)))

    start_hist = collections.Counter()
    stop_hist = collections.Counter()
    frame_counts = [0, 0, 0]
    profile_sum = [0.0] * args.bins
    profile_frame_sum = [[0.0] * args.bins for _ in range(3)]
    profile_cds_n = 0
    metrics: List[Dict[str, object]] = []
    template_stats = TemplateStats.empty()
    kept = 0
    for result in results:
        start_hist.update(result["start_hist"])
        stop_hist.update(result["stop_hist"])
        for idx in range(3):
            frame_counts[idx] += result["frame_counts"][idx]
        for b in range(args.bins):
            profile_sum[b] += result["profile_sum"][b]
            for frame in range(3):
                profile_frame_sum[frame][b] += result["profile_frame_sum"][frame][b]
        profile_cds_n += int(result["profile_cds_n"])
        metrics.extend(result["metrics"])
        template_stats.merge(result["template"])
        kept += int(result["kept"])

    profile_total = [x / profile_cds_n if profile_cds_n else 0.0 for x in profile_sum]
    profile_frames = [[x / profile_cds_n if profile_cds_n else 0.0 for x in profile_frame_sum[frame]] for frame in range(3)]
    frame_total = sum(frame_counts)

    # Basic QC outputs do not depend on successful template estimation.
    write_counts(f"{args.out_prefix}.start.meta.tsv", start_hist, args.window)
    write_counts(f"{args.out_prefix}.stop.meta.tsv", stop_hist, args.window)
    write_frame(f"{args.out_prefix}.frame.tsv", frame_counts)
    write_profile(f"{args.out_prefix}.cds_profile.tsv", profile_total)
    for frame in range(3):
        write_profile(f"{args.out_prefix}.cds_profile.f{frame}.tsv", profile_frames[frame])
    write_metrics(f"{args.out_prefix}.cds_metrics.tsv", metrics)

    template: Optional[Dict[str, object]] = None
    template_error = ""
    try:
        template = estimate_template(template_stats, args.pi_method)
    except ValueError as exc:
        template_error = str(exc)

    template_status = "failed"
    if template is not None:
        pi_used = template["pi_used"]
        alpha0 = template["alpha0"]
        a0 = float(template["A0"])
        if a0 < args.min_a0:
            template_status = "failed_low_a0"
            template_error = f"A0={a0:.10g} is below --min-a0={args.min_a0:.10g}"
        else:
            template_status = "success"
        template_rows: List[Tuple[str, object]] = [
            ("program", PROGRAM), ("version", VERSION), ("author", AUTHOR), ("contact", CONTACT),
            ("template_status", template_status), ("template_error", template_error),
            ("pi_method", template["pi_method"]),
            ("pi_read_weighted_0", float(template["pi_read_weighted"][0])),
            ("pi_read_weighted_1", float(template["pi_read_weighted"][1])),
            ("pi_read_weighted_2", float(template["pi_read_weighted"][2])),
            ("pi_codon_equal_0", float(template["pi_codon_equal"][0])),
            ("pi_codon_equal_1", float(template["pi_codon_equal"][1])),
            ("pi_codon_equal_2", float(template["pi_codon_equal"][2])),
            ("pi_used_0", float(pi_used[0])), ("pi_used_1", float(pi_used[1])), ("pi_used_2", float(pi_used[2])),
            ("A0", a0), ("min_a0", args.min_a0),
            ("alpha0_0", float(alpha0[0])), ("alpha0_1", float(alpha0[1])), ("alpha0_2", float(alpha0[2])),
            ("template_cds_kept", template["template_cds_kept"]),
            ("template_codons", template["template_codons"]),
            ("template_unique_patterns", template["template_unique_patterns"]),
            ("template_trim_start_nt", args.template_trim_start_nt),
            ("template_trim_stop_nt", args.template_trim_stop_nt),
            ("template_min_codon_reads", args.template_min_codon_reads),
        ]
        write_key_value(f"{args.out_prefix}.template.tsv", template_rows)
        if template_status == "success":
            write_key_value(
                f"{args.out_prefix}.dm_background.tsv",
                [("program", PROGRAM), ("version", VERSION), ("template_status", template_status),
                 ("pi_method", template["pi_method"]),
                 ("pi0", float(pi_used[0])), ("pi1", float(pi_used[1])), ("pi2", float(pi_used[2])),
                 ("pi_read_weighted_0", float(template["pi_read_weighted"][0])),
                 ("pi_read_weighted_1", float(template["pi_read_weighted"][1])),
                 ("pi_read_weighted_2", float(template["pi_read_weighted"][2])),
                 ("pi_codon_equal_0", float(template["pi_codon_equal"][0])),
                 ("pi_codon_equal_1", float(template["pi_codon_equal"][1])),
                 ("pi_codon_equal_2", float(template["pi_codon_equal"][2])),
                 ("A", a0), ("min_a0", args.min_a0),
                 ("alpha0", float(alpha0[0])), ("alpha1", float(alpha0[1])), ("alpha2", float(alpha0[2])),
                 ("template_cds_kept", template["template_cds_kept"]),
                 ("template_codons", template["template_codons"]),
                 ("template_unique_patterns", template["template_unique_patterns"]),
                 ("template_trim_start_nt", args.template_trim_start_nt),
                 ("template_trim_stop_nt", args.template_trim_stop_nt),
                 ("template_min_codon_reads", args.template_min_codon_reads),
                 ("template_min_density", args.min_read_density),
                 ("psite_bam", os.path.abspath(args.psite_bam)),
                 ("torf_spec", os.path.abspath(args.torf)),
                 ("rl_filter", args.rl if args.rl is not None else "pooled"),
                 ("rl_source", args.rl_source)],
            )
        else:
            sys.stderr.write(f"WARN: {template_error}; DM background was not written\n")
    else:
        write_key_value(
            f"{args.out_prefix}.template.tsv",
            [("program", PROGRAM), ("version", VERSION), ("author", AUTHOR), ("contact", CONTACT),
             ("template_status", "failed"), ("template_error", template_error),
             ("pi_method", args.pi_method), ("min_a0", args.min_a0),
             ("template_cds_kept", template_stats.n_cds_kept),
             ("template_codons", template_stats.n_codons),
             ("template_unique_patterns", len(template_stats.patterns)),
             ("template_trim_start_nt", args.template_trim_start_nt),
             ("template_trim_stop_nt", args.template_trim_stop_nt),
             ("template_min_codon_reads", args.template_min_codon_reads)],
        )
        sys.stderr.write(f"WARN: template estimation failed: {template_error}; basic QC outputs were still written\n")

    ternary_rows = summarize_codon_ternary(template_stats)
    ternary_vertex = codon_ternary_vertex_percentages(ternary_rows)
    ternary_mean = codon_ternary_mean_percentages(template_stats)
    if ternary_rows:
        write_codon_ternary(f"{args.out_prefix}.codon_frame_ternary.tsv", ternary_rows)
        plot_codon_ternary(
            args.out_prefix,
            ternary_rows,
            float(template["A0"]) if template is not None else None,
            template_status,
            ternary_mean,
        )

    summary_rows: List[Tuple[str, object]] = [
        ("program", PROGRAM), ("version", VERSION), ("author", AUTHOR), ("contact", CONTACT),
        ("input_torf_files", len(tsvs)), ("mane_cds_loaded", len(records)), ("mane_cds_kept_by_density", kept),
        ("min_read_density", args.min_read_density), ("profile_cds_n", profile_cds_n),
        ("frame0_count", frame_counts[0]), ("frame1_count", frame_counts[1]), ("frame2_count", frame_counts[2]),
        ("frame0_prop", frame_counts[0] / frame_total if frame_total else 0.0),
        ("frame1_prop", frame_counts[1] / frame_total if frame_total else 0.0),
        ("frame2_prop", frame_counts[2] / frame_total if frame_total else 0.0),
        ("template_status", template_status), ("min_a0", args.min_a0),
        ("codon_ternary_grid_percent", CODON_TERNARY_GRID_PERCENT),
        ("codon_ternary_total_codons", template_stats.n_codons),
        ("codon_ternary_occupied_points", len(ternary_rows)),
        ("codon_ternary_P0_vertex_percent", ternary_vertex[0]),
        ("codon_ternary_P1_vertex_percent", ternary_vertex[1]),
        ("codon_ternary_P2_vertex_percent", ternary_vertex[2]),
        ("codon_ternary_mean_P0_percent", ternary_mean[0]),
        ("codon_ternary_mean_P1_percent", ternary_mean[1]),
        ("codon_ternary_mean_P2_percent", ternary_mean[2]),
        ("rl", args.rl if args.rl is not None else "pooled"), ("rl_source", args.rl_source),
    ]
    if template is not None:
        summary_rows.extend([
            ("pi_trans_0", float(pi_used[0])), ("pi_trans_1", float(pi_used[1])), ("pi_trans_2", float(pi_used[2])),
            ("A0", float(template["A0"])),
        ])
    else:
        summary_rows.append(("template_error", template_error))
    if template_status != "success" and template_error and template is not None:
        summary_rows.append(("template_error", template_error))
    for key in ("PIF", "Uniformity", "StartIncrease", "DropOff", "TIS_global_enrich", "TIS_local_score"):
        mean, median, n = summarize_metric(metrics, key)
        summary_rows.extend([(f"{key}_mean", mean), (f"{key}_median", median), (f"{key}_n", n)])
    for key, value in sorted(load_stats.items()):
        summary_rows.append((f"load::{key}", value))
    write_key_value(f"{args.out_prefix}.qc_summary.tsv", summary_rows)
    plot_metagene(args.out_prefix, start_hist, stop_hist, frame_counts, profile_frames, args.window)
    sys.stderr.write(f"[ok] wrote {args.out_prefix}.*\n")
    sys.stderr.write(f"[filter] kept {kept:,}/{len(records):,} MANE CDS at density >= {args.min_read_density}\n")
    if template is not None:
        sys.stderr.write(
            f"[template] status={template_status}, pi_trans=({pi_used[0]:.6f}, {pi_used[1]:.6f}, {pi_used[2]:.6f}), "
            f"A0={template['A0']:.6f}, min_A0={args.min_a0:.6f}\n"
        )
    if template_status != "success":
        sys.stderr.write(f"[qc] background gate failed: {template_error}\n")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
