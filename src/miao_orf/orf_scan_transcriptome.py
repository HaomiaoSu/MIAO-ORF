#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transcriptome-wide ORF scanner with consistent tORF-to-gORF collapsing.

Author: Haomiao Su
Contact: suhaomiao@csu.edu.cn
Version: 1.0.0

Overview
--------
This command-line tool scans transcript sequences reconstructed from a genome
FASTA and a GTF annotation, reports transcript-level ORFs (tORFs), collapses
them into stable genomic ORFs (gORFs), and writes an explicit tORF-to-gORF
membership table. It supports chromosome-wise processing for large annotations
while rebuilding one canonical global tORF table at the end of the run.

Key features
------------
- IGV-friendly genomic coordinates: 1-based, closed intervals.
- Transcript model flags and ORF biotyping:
  CDS, intORF_inframe, intORF_altframe, uoORF, doORF, uORF, dORF, other, ncORF.
- Stable gORF identifiers derived from the complete grouping signature.
- tORF-to-gORF membership backfill into the final <out>.torf.tsv.
- Multiprocessing over transcripts with main-process streaming writes.
- Real-time gORF ID collision detection during chromosome-wise appends.
- Final cross-file validation of tORF, gORF, and membership outputs.

Outputs
-------
  <out>.torf.tsv
  <out>.torf.faa
  <out>.gorf.tsv
  <out>.gorf_members.tsv
  <out>.gorf.faa
  <out>.gorf_validation.tsv

Length-filtering invariant
--------------------------
--min-aa is the only ORF-length cutoff. It is applied before a tORF is written,
and the exact emitted tORF set is collapsed into gORFs. No second independent
length filter is applied during collapsing.

Dependencies
------------
pyfaidx
"""
from __future__ import annotations
import argparse, atexit, os, re, sys, csv, hashlib, json
from typing import Any, Dict, List, Optional, Sequence, Tuple, TYPE_CHECKING
from collections import defaultdict, Counter
from multiprocessing import Pool

# Keep pyfaidx optional at module import time so --version and --help remain
# usable before dependencies are installed. Static type checkers still see the
# concrete Fasta type, while runtime annotation resolution remains safe.
if TYPE_CHECKING:
    from pyfaidx import Fasta
else:
    Fasta = Any

__author__ = "Haomiao Su"
__contact__ = "suhaomiao@csu.edu.cn"
__version__ = "1.0.0"
__program__ = "miao-orf-orfscan"

_PROCESS_FASTA = None
_PROCESS_FASTA_PATH: Optional[str] = None


def close_process_fasta() -> None:
    """Close the FASTA handle owned by the current process, if any."""
    global _PROCESS_FASTA, _PROCESS_FASTA_PATH
    if _PROCESS_FASTA is not None:
        close = getattr(_PROCESS_FASTA, "close", None)
        if callable(close):
            close()
    _PROCESS_FASTA = None
    _PROCESS_FASTA_PATH = None


def init_fasta_worker(fa_path: str) -> None:
    """Open one pyfaidx handle per process for reuse across transcripts."""
    global _PROCESS_FASTA, _PROCESS_FASTA_PATH
    normalized_path = os.path.abspath(os.path.expanduser(fa_path))
    if _PROCESS_FASTA is not None and _PROCESS_FASTA_PATH == normalized_path:
        return
    close_process_fasta()
    try:
        from pyfaidx import Fasta
    except ImportError as exc:
        raise RuntimeError(
            "Missing dependency 'pyfaidx'. Install it in the active Python "
            "environment before running the scanner."
        ) from exc
    _PROCESS_FASTA = Fasta(normalized_path, as_raw=True, sequence_always_upper=True)
    _PROCESS_FASTA_PATH = normalized_path


def get_process_fasta(fa_path: str):
    """Return the reusable FASTA handle local to the current process."""
    normalized_path = os.path.abspath(os.path.expanduser(fa_path))
    if _PROCESS_FASTA is None or _PROCESS_FASTA_PATH != normalized_path:
        init_fasta_worker(normalized_path)
    return _PROCESS_FASTA


atexit.register(close_process_fasta)

# --------------------------- CLI ---------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        prog=__program__,
        description=(
            "MIAO transcriptome-wide ORF scanner with consistent tORF-to-gORF "
            "collapsing and validation."
        ),
        epilog=f"Author: {__author__} | Contact: {__contact__}",
    )
    ap.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    ap.add_argument("--gtf", required=True, help="GTF (GENCODE/Ensembl-like)")
    ap.add_argument("--fa", required=True, help="Genome FASTA (indexed with .fai)")
    ap.add_argument("--out-prefix", required=True, help="Output prefix")

    ap.add_argument("--min-aa", type=int, default=15, help="Minimum ORF length (aa) applied before any tORF is written. Default: 15 (45 nt). This is the single ORF-length cutoff used by both tORF output and gORF collapsing.")
    ap.add_argument("--start", nargs="+", default=["ATG"], help="Start codons for scanning. Default: ATG")
    ap.add_argument("--near-cognate", action="store_true", help="Also allow CTG/TTG/GTG as start (first codon→M)")
    ap.add_argument("--primary-only", action="store_true", help="Keep primary assembly contigs only")
    ap.add_argument("--max-transcripts", type=int, default=0, help="Process at most N transcripts (debug)")
    ap.add_argument("--flank-nt", type=int, default=30, help="Transcript flank length to export for 5'/3' (default 30nt)")
    ap.add_argument(
        "--skip-peptide-fasta",
        action="store_true",
        help=(
            "Do not write redundant *.torf.faa/*.gorf.faa files. Peptide sequences remain in the TSV outputs; "
            "tORF-to-gORF collapsing and validation are unchanged. Useful for storage-limited temporary runs."
        ),
    )

    # collapsing options
    ap.add_argument("--no-collapse", action="store_true", help="Skip tORF→gORF collapsing")
    ap.add_argument("--collapse-merge-by", choices=["blocks","stop"], default="blocks",
                    help="Grouping key for collapsing: exact genomic blocks (default) or genomic stop position+genomic frame ('stop').")

    # multiprocessing
    ap.add_argument("--workers", type=int, default=1, help="Worker processes (>=2 to enable multiprocessing)")
    ap.add_argument("--mp-chunksize", type=int, default=50, help="Transcripts per task for ordered multiprocessing imap")
    
    ap.add_argument("--by-chrom", action="store_true",
                help="Process per chromosome: write per-chrom tORF/FAA (and per-chrom updated tORF with gorf_id); append gORFs to one global set of files.")
    ap.add_argument("--chroms", nargs="+", default=None,
                help="Optional subset of chromosomes to process (e.g. chr1 chr2). Defaults to all kept by FASTA whitelist.")
    ap.add_argument("--perchrom-outdir", default="",
                help="Optional output directory for per-chrom files. If empty, files are written next to --out-prefix.")

    return ap.parse_args()

# -------------------- FASTA index helpers --------------------

def read_fai_lengths(fa_path: str):
    fai_path = fa_path + ".fai"
    if not os.path.exists(fai_path):
        return None
    chroms = {}
    with open(fai_path, "r") as f:
        for line in f:
            if not line.strip(): continue
            parts = line.split("\t")
            chroms[parts[0]] = int(parts[1])
    return chroms

_PAT_ALT = re.compile(r"(_alt|_fix|_ALT|_FIX|_random|_decoy|_PATCH|_patch|chrUn|Un_|_UNLOCALIZED)", re.I)
_PAT_UNLOC = re.compile(r"(chrUn|_random|unlocalized|_UNLOCALIZED)", re.I)

def classify_contig(name: str) -> str:
    n = name
    if _PAT_ALT.search(n):
        if "patch" in n.lower() or "fix" in n.lower():
            return "patch"
        return "alt"
    if _PAT_UNLOC.search(n):
        return "unlocalized"
    base = n[3:] if n.startswith("chr") else n
    if base in [str(i) for i in range(1,23)] + ["X","Y","M","MT"]:
        return "primary"
    return "other"

def build_chrom_whitelist(fa_path: str, primary_only: bool):
    fai = read_fai_lengths(fa_path)
    if fai is None:
        sys.exit("ERROR: FASTA index (.fai) not found. Run: samtools faidx <fa>")
    chroms = {nm: {"length": ln, "class": classify_contig(nm)} for nm, ln in fai.items()}
    keep = {k for k, v in chroms.items()} if not primary_only else {k for k, v in chroms.items() if v["class"]=="primary"}
    if not keep:
        sys.exit("ERROR: no chromosomes kept after filtering. Check --primary-only and contig names.")
    return keep

# ------------------------ GTF parsing ------------------------

_ATTR_KV = re.compile(r'\s*([A-Za-z0-9_]+)\s+"([^"]*)"')

def parse_attr(s: str) -> Dict[str, str]:
    d = {}
    for m in _ATTR_KV.finditer(s):
        d[m.group(1)] = m.group(2)
    return d

def parse_tags_from_attr(s: str):
    return [m.group(1) for m in re.finditer(r'\btag\s+"([^\"]+)"', s)]

def gtf_iter(path):
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("#"): continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            yield chrom, feature, int(start)-1, int(end), strand, attrs, parse_attr(attrs)

class TxRec:
    __slots__ = ("gene_id","gene_name","gene_biotype","transcript_id","transcript_biotype",
                 "transcript_name","chrom","strand","exons","cds","start_codon","stop_codon","tags")
    def __init__(self, gid, gname, gbi, tid, tbi, tname, chrom, strand):
        self.gene_id = gid
        self.gene_name = gname or ""
        self.gene_biotype = gbi or ""
        self.transcript_id = tid
        self.transcript_biotype = tbi or ""
        self.transcript_name = tname or ""
        self.chrom = chrom
        self.strand = strand
        self.exons: List[Tuple[int,int]] = []
        self.cds:   List[Tuple[int,int]] = []
        self.start_codon: List[Tuple[int,int]] = []
        self.stop_codon:  List[Tuple[int,int]] = []
        self.tags: List[str] = []

def build_indices(gtf_path: str, chrom_whitelist: set):
    gene_dict = {}
    tx_map: Dict[str, TxRec] = {}
    for chrom, feature, gstart, gend, strand, raw_attr, a in gtf_iter(gtf_path):
        if chrom not in chrom_whitelist: continue
        gid = a.get("gene_id"); tid = a.get("transcript_id")
        if feature == "gene":
            gname = a.get("gene_name",""); gbi = a.get("gene_biotype", a.get("gene_type",""))
            if gid and gid not in gene_dict:
                gene_dict[gid] = {"gene_name": gname, "gene_biotype": gbi, "chrom": chrom, "strand": strand, "transcripts": []}
            continue
        if gid is None or tid is None: continue
        if gid not in gene_dict:
            gene_dict[gid] = {"gene_name": a.get("gene_name",""), "gene_biotype": a.get("gene_biotype", a.get("gene_type","")),
                              "chrom": chrom, "strand": strand, "transcripts": []}
        if feature == "transcript" and tid not in tx_map:
            tbi = a.get("transcript_biotype", a.get("transcript_type",""))
            tname = a.get("transcript_name", "")
            tx_map[tid] = TxRec(gid, gene_dict[gid]["gene_name"], gene_dict[gid]["gene_biotype"],
                                 tid, tbi, tname, chrom, strand)
            gene_dict[gid]["transcripts"].append(tid)
            tags = parse_tags_from_attr(raw_attr); tx_map[tid].tags.extend(tags)
        if tid not in tx_map:
            tbi = a.get("transcript_biotype", a.get("transcript_type",""))
            tname = a.get("transcript_name", "")
            tx_map[tid] = TxRec(gid, gene_dict[gid]["gene_name"], gene_dict[gid]["gene_biotype"],
                                 tid, tbi, tname, chrom, strand)
            gene_dict[gid]["transcripts"].append(tid)
        tx = tx_map[tid]
        if feature == "transcript":
            tags = parse_tags_from_attr(raw_attr)
            if tags: tx.tags.extend([t for t in tags if t not in tx.tags])
        elif feature == "exon":
            tx.exons.append((gstart, gend))
        elif feature == "CDS":
            tx.cds.append((gstart, gend))
        elif feature == "start_codon":
            tx.start_codon.append((gstart, gend))
        elif feature == "stop_codon":
            tx.stop_codon.append((gstart, gend))
    tx_map = {tid: tx for tid, tx in tx_map.items() if tx.exons}
    for gid in list(gene_dict.keys()):
        gene_dict[gid]["transcripts"] = [tid for tid in gene_dict[gid]["transcripts"] if tid in tx_map]
    return gene_dict, tx_map


def genomic_stop_signature(chrom: str, strand: str, positions: Sequence[int]):
    """Return an exact genomic signature for one complete 3-nt stop codon."""
    unique_positions = tuple(sorted(set(int(p) for p in positions)))
    if strand not in {"+", "-"} or len(unique_positions) != 3:
        return None
    return (chrom, strand, unique_positions)


def build_annotated_stop_indexes(tx_map: Dict[str, TxRec]):
    """Index explicit GTF stop_codon features from transcripts carrying a CDS."""
    by_gene = defaultdict(set)
    all_stops = set()
    invalid_stop_transcripts = 0
    for tx in tx_map.values():
        if not tx.cds or not tx.stop_codon:
            continue
        positions = [p for start, end in tx.stop_codon for p in range(start, end)]
        signature = genomic_stop_signature(tx.chrom, tx.strand, positions)
        if signature is None:
            invalid_stop_transcripts += 1
            continue
        by_gene[tx.gene_id].add(signature)
        all_stops.add(signature)
    return dict(by_gene), all_stops, invalid_stop_transcripts


def annotate_intorf_altframe_stop_matches(
    rows: List[Dict[str, object]],
    annotated_stops_by_gene: Dict[str, set],
    all_annotated_stops: set,
) -> None:
    """Annotate only intORF_altframe rows against explicit annotated CDS stops."""
    for row in rows:
        row["intorf_shares_annotated_stop_same_gene"] = 0
        row["intorf_shares_annotated_stop_any_gene"] = 0
        row["intorf_annotated_stop_confounded"] = 0
        signature = row.pop("_intorf_stop_signature", None)
        if row.get("orf_biotype") != "intORF_altframe" or signature is None:
            continue
        same_gene = signature in annotated_stops_by_gene.get(str(row.get("gene_id", "")), set())
        any_gene = signature in all_annotated_stops
        row["intorf_shares_annotated_stop_same_gene"] = int(same_gene)
        row["intorf_shares_annotated_stop_any_gene"] = int(any_gene)
        row["intorf_annotated_stop_confounded"] = int(any_gene)

# ---------------- Transcript utilities & mapping ----------------

def build_transcript_pos_list(exons_genomic: List[Tuple[int,int]], strand: str) -> List[int]:
    if not exons_genomic: return []
    ex_sorted = sorted(exons_genomic, key=lambda x: x[0])
    pos_list: List[int] = []
    for gs, ge in ex_sorted: pos_list.extend(range(gs, ge))
    if strand == "-": pos_list.reverse()
    return pos_list

def build_pos_index(pos_list: List[int]) -> Dict[int,int]:
    return {gpos: idx for idx, gpos in enumerate(pos_list)}

def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def splice_transcript_sequence(fa: Fasta, chrom: str, strand: str, exons_genomic: List[Tuple[int,int]]) -> str:
    ex_sorted = sorted(exons_genomic, key=lambda x: x[0])
    seq_parts = []
    if strand == "+":
        for gs, ge in ex_sorted: seq_parts.append(str(fa[chrom][gs:ge]))
    else:
        for gs, ge in reversed(ex_sorted): seq_parts.append(reverse_complement(str(fa[chrom][gs:ge])))
    return "".join(seq_parts).upper()

# --------- transcript model, mapping to tx coords, helpers ---------

def _compress_sorted_indices_to_ranges(sorted_idx: List[int]) -> List[Tuple[int,int]]:
    if not sorted_idx: return []
    ranges = []; run_start = sorted_idx[0]; prev = sorted_idx[0]
    for x in sorted_idx[1:]:
        if x == prev + 1: prev = x; continue
        ranges.append((run_start, prev + 1)); run_start = x; prev = x
    ranges.append((run_start, prev + 1)); return ranges

def merge_ranges(ranges: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    if not ranges: return []
    rs = sorted(ranges, key=lambda x: x[0])
    out = [rs[0]]
    for s,e in rs[1:]:
        ps,pe = out[-1]
        if s <= pe: out[-1] = (ps, max(pe, e))
        else: out.append((s,e))
    return out

def map_intervals_by_poslist(intervals: List[Tuple[int,int]], pos2tidx: Dict[int,int]) -> List[Tuple[int,int]]:
    out: List[Tuple[int,int]] = []
    for s, e in intervals:
        if e <= s: continue
        idxs = [pos2tidx[p] for p in range(s, e) if p in pos2tidx]
        if not idxs: continue
        idxs = sorted(set(idxs))
        out.extend(_compress_sorted_indices_to_ranges(idxs))
    return merge_ranges(out)

STOP = {"TAA","TAG","TGA"}
FORWARD_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def scan_orfs(seq: str, start_set: set, min_aa: int) -> List[Tuple[int,int,int,str]]:
    L = len(seq); hits: List[Tuple[int,int,int,str]] = []
    for f in (0,1,2):
        i = f; earliest_start = None
        while i + 3 <= L:
            cod = seq[i:i+3]
            if cod in STOP:
                if earliest_start is not None:
                    aa_len = (i - earliest_start) // 3
                    if aa_len >= min_aa:
                        hits.append((earliest_start, i + 3, f, seq[earliest_start:earliest_start+3]))
                earliest_start = None; i += 3; continue
            if cod in start_set and earliest_start is None:
                earliest_start = i
            i += 3
    return hits

def translate_orf(nuc_seq: str, start_idx: int, end_idx: int, start_set: set) -> str:
    aa = []; first = True; i = start_idx
    while i + 3 <= end_idx:
        cod = nuc_seq[i:i+3]
        if first and cod in start_set: aa.append('M')
        else:
            if cod in STOP: break
            aa.append(FORWARD_TABLE.get(cod, 'X'))
        first = False; i += 3
    return "".join(aa)

# ---------- Fast-path helpers for contiguous CDS ----------

def is_contiguous_union(ranges: List[Tuple[int,int]]) -> bool:
    return bool(ranges) and (len(ranges) == 1)

def span_from_union(ranges: List[Tuple[int,int]]):
    return ranges[0] if is_contiguous_union(ranges) else None

# ---------------------- IGV lifting ----------------------

def lift_orf_to_genome_poslist_igv(chrom: str, strand: str, pos_list: List[int], t_start: int, t_end: int):
    if t_end <= t_start: raise ValueError("Empty ORF range")
    if t_start < 0 or t_end > len(pos_list): raise ValueError("ORF outside transcript")
    step_expect = 1 if strand == "+" else -1
    runs = []; i0 = t_start; prev_g = pos_list[t_start]
    for ti in range(t_start + 1, t_end):
        g = pos_list[ti]
        if g - prev_g != step_expect:
            runs.append((i0, ti-1)); i0 = ti
        prev_g = g
    runs.append((i0, t_end - 1))
    # build 0-based half-open blocks then convert to 1-based closed
    blocks01 = []
    for ts_idx, te_idx in runs:
        g_first = pos_list[ts_idx]; g_last = pos_list[te_idx]
        s0 = min(g_first, g_last); e0 = max(g_first, g_last) + 1
        blocks01.append((s0, e0))
    blocks01.sort(key=lambda x: x[0])
    starts1 = [s+1 for s,e in blocks01]; ends1 = [e for s,e in blocks01]
    chromStart1 = starts1[0]; chromEnd1 = ends1[-1]
    blockSizes = [e-s for s,e in blocks01]
    return {"chrom": chrom, "strand": strand, "chromStart1": chromStart1, "chromEnd1": chromEnd1,
            "genomic_block_starts1": ",".join(map(str, starts1)),
            "genomic_block_ends1": ",".join(map(str, ends1)),
            "blockCount": len(blockSizes),
            "blockSizes": ",".join(map(str, blockSizes))}

def lift_tx_window_to_genomic_pieces_igv(pos_list: List[int], strand: str, t_win_start: int, t_win_end: int) -> List[Tuple[int,int]]:
    if t_win_end <= t_win_start: return []
    t_win_start = max(0, t_win_start); t_win_end = min(len(pos_list), t_win_end)
    if t_win_end <= t_win_start: return []
    step_expect = 1 if strand == "+" else -1
    pieces01: List[Tuple[int,int]] = []
    i0 = t_win_start; prev_g = pos_list[i0]
    for ti in range(t_win_start + 1, t_win_end):
        g = pos_list[ti]
        if g - prev_g != step_expect:
            g_first = pos_list[i0]; g_last = prev_g
            s0 = min(g_first, g_last); e0 = max(g_first, g_last) + 1
            pieces01.append((s0, e0)); i0 = ti
        prev_g = g
    g_first = pos_list[i0]; g_last = pos_list[t_win_end-1]
    s0 = min(g_first, g_last); e0 = max(g_first, g_last) + 1
    pieces01.append((s0, e0))
    pieces01.sort(key=lambda x: x[0])
    return [(s0+1, e0) for (s0, e0) in pieces01]

def pieces_to_csv_lists_igv(pieces1: List[Tuple[int,int]]):
    if not pieces1: return "", ""
    return ",".join(str(s) for s,_ in pieces1), ",".join(str(e) for _,e in pieces1)

# ---------------- transcript model & ORF classification ----------------

def find_inframe_stop_from(tseq: str, start_idx: int) -> int:
    i = start_idx + 3; L = len(tseq)
    while i + 3 <= L:
        if tseq[i:i+3] in STOP: return i + 3
        i += 3
    return -1

def recover_stop_if_missing(tseq: str, cds_t_union: List[Tuple[int,int]], stopcodon_t: List[Tuple[int,int]]):
    if stopcodon_t: return stopcodon_t, False
    if not cds_t_union: return [], False
    t_stop_start = cds_t_union[-1][1]
    if tseq[t_stop_start:t_stop_start+3].upper() in STOP:
        return [(t_stop_start, t_stop_start+3)], True
    return [], False

def derive_transcript_model(tseq, cds_t_union_raw, startcodon_t, stopcodon_t):
    flags = {
        "tx_has_start_annot": 1 if startcodon_t else 0,
        "tx_has_stop_annot":  1 if stopcodon_t else 0,
        "tx_stop_recovered":  0,
        "tx_cds_synthesized": 0,
    }
    if cds_t_union_raw:
        cds_t_union = cds_t_union_raw
        if not stopcodon_t:
            stopcodon_t_rec, rec = recover_stop_if_missing(tseq, cds_t_union, stopcodon_t)
            if rec:
                flags["tx_stop_recovered"] = 1
                stopcodon_t = stopcodon_t_rec
        has_start = bool(startcodon_t); has_stop = bool(stopcodon_t)
        if has_start and has_stop: model = "coding"
        elif not has_start:        model = "non_coding"
        else:                      model = "no_stop"
        return model, cds_t_union, startcodon_t, stopcodon_t, flags

    if startcodon_t:
        flags["tx_cds_synthesized"] = 1
        s0 = min(s for (s, e) in startcodon_t)
        stop_end = find_inframe_stop_from(tseq, s0)
        if stop_end != -1:
            cds_t_union = [(s0, stop_end - 3)]
            stopcodon_t = [(stop_end - 3, stop_end)]
            model = "coding"
        else:
            cds_t_union = [(s0, len(tseq))]
            stopcodon_t = []
            model = "no_stop"
        return model, cds_t_union, startcodon_t, stopcodon_t, flags

    return "non_coding", [], startcodon_t, stopcodon_t, flags

def annotated_cds_trange(cds_t_union, startcodon_t, stopcodon_t):
    if startcodon_t and stopcodon_t:
        t_start = min(s for s, e in startcodon_t)
        t_end   = max(e for s, e in stopcodon_t)
        if t_end > t_start:
            return (t_start, t_end)
    return (None, None)

def interval_overlap_len(s1,e1,s2,e2) -> int:
    return max(0, min(e1,e2) - max(s1,s2))

def sum_overlap_with_union(s: int, e: int, union: List[Tuple[int,int]]) -> int:
    return sum(interval_overlap_len(s,e,us,ue) for us,ue in union)

def is_inside_union(s: int, e: int, union: List[Tuple[int,int]]) -> bool:
    return sum_overlap_with_union(s, e, union) == (e - s)

# ---------------------- Worker ----------------------

def kozak_context(seq: str, start_idx: int):
    left = seq[max(0, start_idx-6): start_idx]
    codon = seq[start_idx:start_idx+3]
    right = seq[start_idx+3: start_idx+3+5]
    left = (left[-6:] if len(left) >= 6 else ("N"*(6-len(left)) + left))
    right = (right[:5] if len(right) >= 5 else (right + "N"*(5-len(right))))
    kseq = left + codon + right
    base_minus3 = kseq[3] if len(kseq) >= 7 else "N"
    base_plus4  = kseq[9] if len(kseq) >= 10 else "N"
    strong = (base_minus3 in "AG") and (base_plus4 == "G")
    moderate = (base_minus3 in "AG") or (base_plus4 == "G")
    klass = "Strong" if strong else ("Moderate" if moderate else "Weak")
    return kseq, klass

# bundle minimal tx payload so worker can re-open FASTA locally
class TxPayload:
    __slots__ = ("tid","chrom","strand","exons","cds","start_codon","stop_codon",
                 "gene_id","gene_name","gene_biotype","transcript_name","transcript_biotype","tags")
    def __init__(self, tid, tx: TxRec):
        self.tid = tid
        self.chrom = tx.chrom; self.strand = tx.strand
        self.exons = tx.exons; self.cds = tx.cds
        self.start_codon = tx.start_codon; self.stop_codon = tx.stop_codon
        self.gene_id = tx.gene_id; self.gene_name = tx.gene_name; self.gene_biotype = tx.gene_biotype
        self.transcript_name = tx.transcript_name; self.transcript_biotype = tx.transcript_biotype
        self.tags = list(tx.tags)

def process_tx(args: Tuple[str, "TxPayload", set, int, int]) -> Tuple[List[Dict[str, str]], List[Tuple[str,str]]]:
    fa_path, payload, start_set, min_aa, flank_nt = args
    fa = get_process_fasta(fa_path)

    tid = payload.tid
    pos_list = build_transcript_pos_list(payload.exons, payload.strand)
    if not pos_list: return [], []
    tx_len = len(pos_list)

    try:
        # splice_transcript_sequence() sorts exons internally and handles strand.
        tseq = splice_transcript_sequence(fa, payload.chrom, payload.strand, payload.exons)
    except KeyError:
        return [], []
    if len(tseq) != tx_len or tx_len < 3:
        return [], []

    # map annotations
    pos2tidx = build_pos_index(pos_list)
    cds_t_union_raw = map_intervals_by_poslist(payload.cds, pos2tidx)
    startcodon_t    = map_intervals_by_poslist(payload.start_codon, pos2tidx) if payload.start_codon else []
    stopcodon_t     = map_intervals_by_poslist(payload.stop_codon,  pos2tidx) if payload.stop_codon  else []

    # normalize model
    tx_model, cds_t_union, startcodon_t, stopcodon_t, tx_flags = derive_transcript_model(
        tseq, cds_t_union_raw, startcodon_t, stopcodon_t
    )
    ann_cds_t_start, ann_cds_t_end = annotated_cds_trange(cds_t_union, startcodon_t, stopcodon_t)
    tx_complete_cds = 1 if (cds_t_union and startcodon_t and stopcodon_t) else 0

    cds_span = span_from_union(cds_t_union)
    if cds_span is not None:
        cds_first_t, cds_last_t = cds_span[0], cds_span[1]
    else:
        if cds_t_union:
            cds_first_t = cds_t_union[0][0]
            cds_last_t  = cds_t_union[-1][1]
        else:
            cds_first_t = None; cds_last_t = None

    tags_str = ";".join(payload.tags) if payload.tags else ""
    has_mane = 1 if any(t == "MANE_Select" for t in payload.tags) else 0

    rows: List[Dict[str,str]] = []
    faa: List[Tuple[str,str]] = []

    for (t_start, t_end, frame, _sc) in scan_orfs(tseq, start_set, min_aa):
        t_start_out, t_end_out = t_start, t_end
        biotype_override = None

        # force to annotated CDS if valid and shares annotated STOP
        force_to_ann = (
            tx_complete_cds == 1 and ann_cds_t_start is not None and ann_cds_t_end is not None and ann_cds_t_end > ann_cds_t_start and t_end == ann_cds_t_end
        )
        if force_to_ann:
            t_start_out = ann_cds_t_start
            t_end_out   = ann_cds_t_end
            biotype_override = "CDS"

        if t_end_out <= t_start_out: continue
        aa_len = (t_end_out - t_start_out)//3 - 1
        if aa_len < min_aa: continue

        pep = translate_orf(tseq, t_start_out, t_end_out, start_set)
        kseq, kclass = kozak_context(tseq, t_start_out)
        start_codon_out = tseq[t_start_out:t_start_out+3]

        # classification
        if cds_span is not None:
            cds_start, cds_end = cds_span
            if tx_model == "non_coding":
                biotype = "ncORF"; ov_cds = ov_u5 = ov_u3 = 0
            else:
                in_cds = (t_start_out >= cds_start and t_end_out <= cds_end)
                if in_cds:
                    rel = (t_start_out - cds_start) % 3
                    biotype = "intORF_inframe" if rel == 0 else "intORF_altframe"
                elif t_end_out <= cds_start:
                    biotype = "uORF"
                elif (ann_cds_t_end is not None) and (t_start_out >= ann_cds_t_end):
                    biotype = "dORF"
                else:
                    biotype = None
                ov_cds = interval_overlap_len(t_start_out, t_end_out, cds_start, cds_end)
                ov_u5  = interval_overlap_len(t_start_out, t_end_out, 0, cds_start)
                ov_u3  = interval_overlap_len(t_start_out, t_end_out, ann_cds_t_end, tx_len) if (ann_cds_t_end is not None) else 0
                if biotype is None:
                    if ov_cds > 0 and ov_u5 > 0: biotype = "uoORF"
                    elif ov_cds > 0 and ov_u3 > 0: biotype = "doORF"
                    else: biotype = "other"
        else:
            if (tx_model == "non_coding") or (not cds_t_union):
                biotype = "ncORF"; ov_cds = ov_u5 = ov_u3 = 0
            else:
                cds_first = cds_t_union[0][0]
                if is_inside_union(t_start_out, t_end_out, cds_t_union):
                    rel = (t_start_out - cds_first) % 3
                    biotype = "intORF_inframe" if rel == 0 else "intORF_altframe"
                elif t_end_out <= cds_first:
                    biotype = "uORF"
                elif (ann_cds_t_end is not None) and (t_start_out >= ann_cds_t_end):
                    biotype = "dORF"
                else:
                    biotype = None
                ov_cds = sum_overlap_with_union(t_start_out, t_end_out, cds_t_union)
                ov_u5  = interval_overlap_len(t_start_out, t_end_out, 0, cds_first)
                ov_u3  = interval_overlap_len(t_start_out, t_end_out, ann_cds_t_end, tx_len) if (ann_cds_t_end is not None) else 0
                if biotype is None:
                    if ov_cds > 0 and ov_u5 > 0: biotype = "uoORF"
                    elif ov_cds > 0 and ov_u3 > 0: biotype = "doORF"
                    else: biotype = "other"

        if biotype_override is not None:
            biotype = biotype_override

        intorf_stop_signature = None
        if biotype == "intORF_altframe":
            intorf_stop_signature = genomic_stop_signature(
                payload.chrom,
                payload.strand,
                pos_list[t_end_out-3:t_end_out],
            )

        # annotations presence
        start_sites = {s for (s, e) in startcodon_t}
        stop_ends   = {e for (s, e) in stopcodon_t}
        is_ann = 1 if (ann_cds_t_start is not None and ann_cds_t_end is not None and t_start_out == ann_cds_t_start and t_end_out == ann_cds_t_end) else 0
        ann_id = f"cds|{tid}|{t_start_out}|{t_end_out}" if is_ann else ""
        orf_has_upstream_start_annot = 1 if any(s < t_start_out for s in start_sites) else 0
        orf_shares_annotated_stop    = 1 if (t_end_out in stop_ends) else 0

        lift = lift_orf_to_genome_poslist_igv(payload.chrom, payload.strand, pos_list, t_start_out, t_end_out)
        torf_id = f"t|{tid}|{t_start_out}|{t_end_out}|f{frame}"

        # flanks (1-based closed)
        t5_s = max(0, t_start_out - max(0, int(flank_nt))); t5_e = t_start_out
        pieces5 = lift_tx_window_to_genomic_pieces_igv(pos_list, payload.strand, t5_s, t5_e)
        f5_starts, f5_ends = pieces_to_csv_lists_igv(pieces5)
        flank5_tx_len_requested = max(0, int(flank_nt))
        flank5_tx_len_actual = max(0, t5_e - t5_s)
        flank5_complete = 1 if flank5_tx_len_actual >= flank5_tx_len_requested else 0
        flank5_piece_count = len(pieces5)

        t3_s = t_end_out; t3_e = min(tx_len, t_end_out + max(0, int(flank_nt)))
        pieces3 = lift_tx_window_to_genomic_pieces_igv(pos_list, payload.strand, t3_s, t3_e)
        f3_starts, f3_ends = pieces_to_csv_lists_igv(pieces3)
        flank3_tx_len_requested = max(0, int(flank_nt))
        flank3_tx_len_actual = max(0, t3_e - t3_s)
        flank3_complete = 1 if flank3_tx_len_actual >= flank3_tx_len_requested else 0
        flank3_piece_count = len(pieces3)

        row = {
            "torf_id": torf_id,
            "gene_id": payload.gene_id, "gene_name": payload.gene_name, "gene_biotype": payload.gene_biotype,
            "transcript_id": tid, "transcript_name": payload.transcript_name,
            "transcript_biotype": payload.transcript_biotype, "tx_tags": tags_str, "has_MANE_Select_tag": has_mane,
            "chrom": payload.chrom, "strand": payload.strand,
            "tx_model": tx_model, "tx_complete_cds": tx_complete_cds,
            "tx_has_start_annot": tx_flags['tx_has_start_annot'], "tx_has_stop_annot": tx_flags['tx_has_stop_annot'],
            "tx_stop_recovered": tx_flags['tx_stop_recovered'], "tx_cds_synthesized": tx_flags['tx_cds_synthesized'],
            "t_start": t_start_out, "t_end": t_end_out, "frame": frame, "aa_len": aa_len,
            "start_codon": start_codon_out, "kozak_seq": kseq, "kozak_class": kclass,
            "orf_biotype": biotype, "overlap_cds_bp": ov_cds, "overlap_utr5_bp": ov_u5, "overlap_utr3_bp": ov_u3,
            "cds_first_t": (cds_span[0] if cds_span else (cds_t_union[0][0] if cds_t_union else "")),
            "cds_last_t":  (cds_span[1] if cds_span else (cds_t_union[-1][1] if cds_t_union else "")),
            "is_annotated_cds": is_ann, "annotated_cds_id": ann_id,
            "orf_has_upstream_start_annot": orf_has_upstream_start_annot,
            "orf_shares_annotated_stop": orf_shares_annotated_stop,
            "intorf_shares_annotated_stop_same_gene": 0,
            "intorf_shares_annotated_stop_any_gene": 0,
            "intorf_annotated_stop_confounded": 0,
            "_intorf_stop_signature": intorf_stop_signature,
            # IGV 1-based closed
            "chromStart1": lift["chromStart1"], "chromEnd1": lift["chromEnd1"],
            "blockCount": lift["blockCount"], "blockSizes": lift["blockSizes"],
            "genomic_block_starts1": lift["genomic_block_starts1"],
            "genomic_block_ends1": lift["genomic_block_ends1"],
            # flanks
            "flank5_genomic_starts1": f5_starts, "flank5_genomic_ends1": f5_ends,
            "flank3_genomic_starts1": f3_starts, "flank3_genomic_ends1": f3_ends,
            "flank5_tx_len_requested": flank5_tx_len_requested, "flank5_tx_len_actual": flank5_tx_len_actual,
            "flank5_complete": flank5_complete, "flank5_piece_count": flank5_piece_count,
            "flank3_tx_len_requested": flank3_tx_len_requested, "flank3_tx_len_actual": flank3_tx_len_actual,
            "flank3_complete": flank3_complete, "flank3_piece_count": flank3_piece_count,
            "peptide_len": len(pep), "peptide": pep,
        }
        rows.append(row)
        faa.append((torf_id, pep))

    return rows, faa

# ---------------------- Collapse helpers ----------------------

BIOTYPE_PRIORITY = [
    "CDS",
    "intORF_inframe",
    "intORF_altframe",
    "uoORF",
    "doORF",
    "uORF",
    "dORF",
    "other",
    "ncORF",
]
PRIO = {b:i for i,b in enumerate(BIOTYPE_PRIORITY)}

def biotype_rank(bt: str) -> int: return PRIO.get(bt, len(BIOTYPE_PRIORITY) + 10)

def safe_int(x, default=None):
    try: return int(x)
    except Exception: return default


def parse_blocks1(block_sizes: str, genomic_starts1: str):
    sizes = tuple(int(x) for x in block_sizes.split(",") if x != "")
    starts1= tuple(int(x) for x in genomic_starts1.split(",") if x != "")
    return sizes, starts1


def torf_blocks_key_igv(row: dict):
    sizes, starts1 = parse_blocks1(row["blockSizes"], row["genomic_block_starts1"])
    return (row["chrom"], row["strand"], sizes, starts1)


def torf_stop_key_igv(row: dict):
    """Stable genomic stop key for optional stop-based collapsing.

    The transcript-level ORF frame stored in row["frame"] is relative to the
    spliced transcript sequence and can differ across transcript isoforms with
    different 5' UTR lengths. It must therefore not be used for genomic ORF
    collapsing. Use the genomic stop coordinate to derive a stable genomic
    frame instead.
    """
    strand = row["strand"]
    chromStart1 = safe_int(row.get("chromStart1", None), None)
    chromEnd1   = safe_int(row.get("chromEnd1", None), None)
    if chromStart1 is None or chromEnd1 is None:
        return None
    stop_pos1 = chromStart1 if strand == "-" else chromEnd1
    stop_pos0 = stop_pos1 - 1
    genomic_stop_frame = stop_pos0 % 3
    return (row["chrom"], strand, stop_pos1, genomic_stop_frame)


def best_record_by_priority(rows: List[dict]) -> dict:
    def key(row):
        bt = row.get("orf_biotype","other")
        pep_len = safe_int(row.get("peptide_len", None), None)
        if pep_len is None: pep_len = len(row.get("peptide",""))
        aa_len = safe_int(row.get("aa_len", 0), 0)
        return (biotype_rank(bt), -pep_len, -aa_len, row.get("torf_id",""))
    return sorted(rows, key=key)[0]


def _canonical_group_signature(merge_by: str, key: tuple) -> str:
    """Stable, complete signature for one gORF grouping key."""
    return json.dumps(
        {"merge_by": merge_by, "key": key},
        sort_keys=True, separators=(",", ":"), ensure_ascii=True,
    )


def _stable_gorf_id(merge_by: str, key: tuple, rep: dict) -> Tuple[str, str]:
    """Create a compact stable gORF ID from the complete grouping key."""
    signature = _canonical_group_signature(merge_by, key)
    digest = hashlib.sha1(signature.encode("utf-8")).hexdigest()[:16]
    chrom = str(rep.get("chrom", ""))
    strand = str(rep.get("strand", ""))
    return f"g|{chrom}|{strand}|{merge_by}|h{digest}", signature


def collapse_torf_rows_igv(rows: List[dict], merge_by: str):
    """Collapse the complete emitted tORF set into gORFs.

    Length filtering is intentionally absent here. The single ORF-length cutoff
    is --min-aa and is applied during transcript scanning before rows are written.
    Therefore every emitted tORF must receive exactly one gORF membership.
    """

    groups = defaultdict(list)
    if merge_by == "stop":
        for row in rows:
            key0 = torf_stop_key_igv(row)
            key = (("stop",) + key0) if key0 is not None else (("fallback_blocks",) + torf_blocks_key_igv(row))
            groups[key].append(row)
    else:
        for row in rows:
            key = ("blocks",) + torf_blocks_key_igv(row)
            groups[key].append(row)

    g_list = []
    members_map = []
    faa_list = []
    torf2gorf = {}
    seen_gorf_signatures = {}

    for key, mem in groups.items():
        rep = best_record_by_priority(mem)
        gene_ids   = [m.get("gene_id","") for m in mem if m.get("gene_id")]
        gene_names = [m.get("gene_name","") for m in mem if m.get("gene_name")]
        gene_id = Counter(gene_ids).most_common(1)[0][0] if gene_ids else rep.get("gene_id","")
        gene_name = Counter(gene_names).most_common(1)[0][0] if gene_names else rep.get("gene_name","")

        biotype_counts = Counter(m.get("orf_biotype","other") for m in mem)
        biotype_summary = ",".join(f"{k}:{biotype_counts[k]}" for k in sorted(biotype_counts.keys(), key=biotype_rank))
        intorf_shared_stop_same_gene_any = int(any(
            safe_int(m.get("intorf_shares_annotated_stop_same_gene", 0), 0) == 1 for m in mem
        ))
        intorf_shared_stop_any_gene_any = int(any(
            safe_int(m.get("intorf_shares_annotated_stop_any_gene", 0), 0) == 1 for m in mem
        ))
        intorf_annotated_stop_confounded_any = int(any(
            safe_int(m.get("intorf_annotated_stop_confounded", 0), 0) == 1 for m in mem
        ))

        txs = sorted(set(m.get("transcript_id","") for m in mem if m.get("transcript_id")))
        torf_ids = [str(m["torf_id"]).strip() for m in mem]

        pep = rep.get("peptide","")
        pep_len = safe_int(rep.get("peptide_len", None), None)
        if pep_len is None: pep_len = len(pep)

        chrom = rep["chrom"]; strand = rep["strand"]
        chromStart1 = safe_int(rep["chromStart1"], 0)
        chromEnd1   = safe_int(rep["chromEnd1"],   0)

        gorf_id, signature = _stable_gorf_id(merge_by, key, rep)
        old_sig = seen_gorf_signatures.get(gorf_id)
        if old_sig is not None and old_sig != signature:
            raise RuntimeError(f"gORF hash collision detected for {gorf_id}")
        seen_gorf_signatures[gorf_id] = signature

        g_list.append({
            "gorf_id": gorf_id,
            "gorf_signature": signature,
            "merge_by": merge_by,
            "gene_id": gene_id,
            "gene_name": gene_name,
            "chrom": chrom,
            "strand": strand,
            "chromStart1": chromStart1,
            "chromEnd1": chromEnd1,
            "blockCount": rep["blockCount"],
            "blockSizes": rep["blockSizes"],
            "genomic_block_starts1": rep["genomic_block_starts1"],
            "genomic_block_ends1": rep["genomic_block_ends1"],
            "orf_biotype_rep": rep.get("orf_biotype",""),
            "biotype_summary": biotype_summary,
            "intorf_shared_annotated_stop_same_gene_any": intorf_shared_stop_same_gene_any,
            "intorf_shared_annotated_stop_any_gene_any": intorf_shared_stop_any_gene_any,
            "intorf_annotated_stop_confounded_any": intorf_annotated_stop_confounded_any,
            "peptide_len": pep_len,
            "peptide": pep,
            "n_members": len(mem),
            "n_transcripts": len(txs),
            "transcripts": ",".join(txs),
            "torf_ids": ",".join(torf_ids),
        })
        for tid in torf_ids:
            previous = torf2gorf.get(tid)
            if previous is not None and previous != gorf_id:
                raise RuntimeError(f"tORF {tid} assigned to multiple gORFs: {previous} vs {gorf_id}")
            torf2gorf[tid] = gorf_id
            members_map.append({"gorf_id": gorf_id, "torf_id": tid})
        if pep:
            faa_list.append((gorf_id, pep))

    if len(torf2gorf) != len(rows):
        raise RuntimeError(
            f"Incomplete collapse mapping: input tORFs={len(rows):,}, mapped unique tORFs={len(torf2gorf):,}"
        )

    return g_list, members_map, faa_list, torf2gorf


def write_gorf_outputs(
    out_prefix: str,
    g_list: List[dict],
    members_map: List[dict],
    faa_list: List[Tuple[str,str]],
    write_faa: bool = True,
):
    out_tsv = f"{out_prefix}.gorf.tsv"
    out_mem = f"{out_prefix}.gorf_members.tsv"
    out_faa = f"{out_prefix}.gorf.faa"

    tsv_cols = [
        "gorf_id","gorf_signature","merge_by","gene_id","gene_name","chrom","strand",
        "chromStart1","chromEnd1","blockCount","blockSizes","genomic_block_starts1","genomic_block_ends1",
        "orf_biotype_rep","biotype_summary",
        "intorf_shared_annotated_stop_same_gene_any","intorf_shared_annotated_stop_any_gene_any",
        "intorf_annotated_stop_confounded_any","peptide_len","peptide",
        "n_members","n_transcripts","transcripts","torf_ids"
    ]
    with open(out_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=tsv_cols, delimiter="\t")
        w.writeheader()
        for row in g_list:
            w.writerow(row)

    with open(out_mem, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gorf_id","torf_id"], delimiter="\t")
        w.writeheader()
        for m in members_map:
            w.writerow(m)

    if write_faa:
        with open(out_faa, "w", encoding="utf-8") as f:
            for gid, pep in faa_list:
                f.write(f">{gid}\n{pep}\n")

# ------------------------------- Main -------------------------------



def write_gorf_outputs_append(
    out_prefix: str,
    g_list: List[dict],
    members_map: List[dict],
    faa_list: List[Tuple[str, str]],
    seen_global_gorf_signatures: Dict[str, Tuple[str, str]],
    source_label: str,
    write_faa: bool = True,
) -> None:
    """Append chromosome-level gORFs to the global outputs with early checks.

    Parameters
    ----------
    seen_global_gorf_signatures
        Mutable registry mapping gorf_id to (signature, first_source_label).
        It persists across chromosome shards so hash collisions or duplicate
        identifiers are detected before anything is appended to disk.
    source_label
        Human-readable shard label, normally the chromosome name.
    """
    import csv, os
    out_tsv = f"{out_prefix}.gorf.tsv"
    out_mem = f"{out_prefix}.gorf_members.tsv"
    out_faa = f"{out_prefix}.gorf.faa"

    tsv_cols = [
        "gorf_id","gorf_signature","merge_by","gene_id","gene_name","chrom","strand",
        "chromStart1","chromEnd1","blockCount","blockSizes","genomic_block_starts1","genomic_block_ends1",
        "orf_biotype_rep","biotype_summary",
        "intorf_shared_annotated_stop_same_gene_any","intorf_shared_annotated_stop_any_gene_any",
        "intorf_annotated_stop_confounded_any","peptide_len","peptide",
        "n_members","n_transcripts","transcripts","torf_ids"
    ]

    # Detect conflicts before appending any part of this shard. This provides a
    # clearer error than the final validation alone and prevents partial writes.
    local_seen: Dict[str, str] = {}
    for row in g_list:
        gid = str(row.get("gorf_id", "")).strip()
        signature = str(row.get("gorf_signature", "")).strip()
        if not gid or not signature:
            raise RuntimeError(
                f"Blank gORF identifier/signature in shard {source_label}: "
                f"gorf_id={gid!r}, signature={signature!r}"
            )
        local_previous = local_seen.get(gid)
        if local_previous is not None:
            if local_previous == signature:
                raise RuntimeError(
                    f"Duplicate gorf_id within shard {source_label}: {gid}. "
                    "The same gORF row would be appended twice."
                )
            raise RuntimeError(
                f"gORF hash collision within shard {source_label}: {gid}. "
                f"signature_1={local_previous}; signature_2={signature}"
            )
        local_seen[gid] = signature

        global_previous = seen_global_gorf_signatures.get(gid)
        if global_previous is not None:
            previous_signature, first_source = global_previous
            if previous_signature == signature:
                raise RuntimeError(
                    f"Duplicate gorf_id across chromosome shards: {gid}. "
                    f"first_seen_in={first_source}; seen_again_in={source_label}. "
                    "The same gORF row would be appended twice."
                )
            raise RuntimeError(
                f"gORF hash collision across chromosome shards: {gid}. "
                f"first_seen_in={first_source}; seen_again_in={source_label}; "
                f"signature_1={previous_signature}; signature_2={signature}"
            )

    # Update the persistent registry only after the complete shard passes.
    for gid, signature in local_seen.items():
        seen_global_gorf_signatures[gid] = (signature, source_label)

    # TSV (gorf)
    header_exists = os.path.exists(out_tsv)
    with open(out_tsv, "a", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=tsv_cols, delimiter="\t")
        if not header_exists:
            w.writeheader()
        for row in g_list:
            w.writerow(row)

    # TSV (members)
    mem_header_exists = os.path.exists(out_mem)
    with open(out_mem, "a", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gorf_id","torf_id"], delimiter="\t")
        if not mem_header_exists:
            w.writeheader()
        for m in members_map:
            w.writerow(m)

    # FAA is optional for storage-limited temporary runs. The peptide remains
    # present in the gORF TSV and is not part of the cross-file validation.
    if write_faa:
        with open(out_faa, "a", encoding="utf-8") as f:
            for gid, pep in faa_list:
                f.write(f">{gid}\n{pep}\n")



def update_torf_tsv_with_gorf_id(torf_tsv: str, torf2gorf: Dict[str, str]) -> None:
    """Backfill gorf_id and fail if any emitted tORF lacks a mapping."""
    tmp_path = torf_tsv + ".tmp"
    total = 0
    missing = []
    with open(torf_tsv, "r", encoding="utf-8") as f_in, open(tmp_path, "w", encoding="utf-8", newline="") as f_out:
        r = csv.DictReader(f_in, delimiter="\t")
        fieldnames = list(r.fieldnames or [])
        if "torf_id" not in fieldnames:
            raise RuntimeError(f"Missing torf_id column in {torf_tsv}")
        if "gorf_id" not in fieldnames:
            fieldnames.append("gorf_id")
        w = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in r:
            total += 1
            tid = str(row.get("torf_id", "")).strip()
            gid = torf2gorf.get(tid, "")
            if not gid:
                if len(missing) < 10:
                    missing.append(tid)
            row["gorf_id"] = gid
            w.writerow(row)
    if missing:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise RuntimeError(
            f"Incomplete tORF→gORF mapping while updating {torf_tsv}: "
            f"examples={missing}. Use --min-aa for filtering; do not filter only during collapse."
        )
    os.replace(tmp_path, torf_tsv)
    sys.stderr.write(f"[ok] updated {torf_tsv} with gorf_id for {total:,} tORFs\n")


def merge_tsv_files_with_identical_header(input_paths: List[str], output_path: str) -> None:
    """Merge per-chrom TSV shards after gorf_id backfill."""
    header = None
    total = 0
    with open(output_path, "w", encoding="utf-8", newline="") as f_out:
        for path in input_paths:
            with open(path, "r", encoding="utf-8") as f_in:
                r = csv.reader(f_in, delimiter="\t")
                try:
                    current_header = next(r)
                except StopIteration:
                    continue
                if header is None:
                    header = current_header
                    f_out.write("\t".join(header) + "\n")
                elif current_header != header:
                    raise RuntimeError(f"TSV header mismatch while merging: {path}")
                for fields in r:
                    f_out.write("\t".join(fields) + "\n")
                    total += 1
    if header is None:
        raise RuntimeError("No per-chrom tORF TSV rows available for merge")
    sys.stderr.write(f"[all] [ok] merged {len(input_paths)} tORF shards into {output_path} ({total:,} rows)\n")


def concatenate_text_files(input_paths: List[str], output_path: str) -> None:
    with open(output_path, "w", encoding="utf-8") as f_out:
        for path in input_paths:
            with open(path, "r", encoding="utf-8") as f_in:
                for line in f_in:
                    f_out.write(line)
    sys.stderr.write(f"[all] [ok] merged {len(input_paths)} FASTA shards into {output_path}\n")


def validate_gorf_outputs(torf_tsv: str, gorf_tsv: str, members_tsv: str, report_path: str) -> None:
    """Require exact agreement among tORF, gORF and membership files."""
    gorf_ids = set()
    with open(gorf_tsv, "r", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            gid = str(row.get("gorf_id", "")).strip()
            if not gid:
                raise RuntimeError(f"Blank gorf_id in {gorf_tsv}")
            if gid in gorf_ids:
                raise RuntimeError(f"Duplicate gorf_id in {gorf_tsv}: {gid}")
            gorf_ids.add(gid)

    member_map: Dict[str, str] = {}
    with open(members_tsv, "r", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            gid = str(row.get("gorf_id", "")).strip()
            tid = str(row.get("torf_id", "")).strip()
            if not gid or not tid:
                raise RuntimeError(f"Blank gorf_id or torf_id in {members_tsv}")
            if gid not in gorf_ids:
                raise RuntimeError(f"Membership references absent gORF: {gid}")
            previous = member_map.get(tid)
            if previous is not None and previous != gid:
                raise RuntimeError(f"tORF {tid} maps to multiple gORFs in members TSV")
            member_map[tid] = gid

    torf_ids = set()
    missing_member = []
    mismatch = []
    with open(torf_tsv, "r", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        fields = set(r.fieldnames or [])
        if "gorf_id" not in fields:
            raise RuntimeError(f"Missing gorf_id column in {torf_tsv}")
        for row in r:
            tid = str(row.get("torf_id", "")).strip()
            gid = str(row.get("gorf_id", "")).strip()
            if not tid or not gid:
                if len(missing_member) < 10:
                    missing_member.append((tid, gid))
                continue
            if tid in torf_ids:
                raise RuntimeError(f"Duplicate torf_id in {torf_tsv}: {tid}")
            torf_ids.add(tid)
            if member_map.get(tid) != gid and len(mismatch) < 10:
                mismatch.append((tid, gid, member_map.get(tid, "")))
    if missing_member:
        raise RuntimeError(f"Blank tORF/gORF mapping in {torf_tsv}; examples={missing_member}")
    if mismatch:
        raise RuntimeError(f"tORF TSV and members TSV disagree; examples={mismatch}")
    extra_members = sorted(set(member_map) - torf_ids)
    missing_members = sorted(torf_ids - set(member_map))
    if extra_members or missing_members:
        raise RuntimeError(
            f"Membership set mismatch: extra_members={len(extra_members):,}, missing_members={len(missing_members):,}"
        )

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"program\t{__program__}\n")
        f.write(f"version\t{__version__}\n")
        f.write(f"author\t{__author__}\n")
        f.write(f"contact\t{__contact__}\n")
        f.write(f"torf_rows\t{len(torf_ids)}\n")
        f.write(f"gorf_rows\t{len(gorf_ids)}\n")
        f.write(f"membership_rows\t{len(member_map)}\n")
        f.write("all_torf_have_gorf_id\t1\n")
        f.write("torf_membership_sets_identical\t1\n")
        f.write("validation_passed\t1\n")
    sys.stderr.write(f"[validate] [ok] exact tORF↔gORF consistency verified; report: {report_path}\n")

def main():
    a = parse_args()

    if a.min_aa < 1:
        sys.exit("ERROR: --min-aa must be >= 1")
    sys.stderr.write(f"[init] single ORF-length cutoff: --min-aa {a.min_aa} aa (applied before tORF output and shared by gORF collapse)\n")

    start_set = {x.upper() for x in a.start}
    if a.near_cognate:
        start_set |= {"CTG","TTG","GTG"}


    sys.stderr.write("[init] building chromosome whitelist from FASTA index...\n")
    chrom_whitelist = build_chrom_whitelist(a.fa, a.primary_only)
    sys.stderr.write(f"[ok] contigs kept: {len(chrom_whitelist)}\n")

    sys.stderr.write("[init] parsing GTF...\n")
    gene_dict, tx_map = build_indices(a.gtf, chrom_whitelist)
    sys.stderr.write(f"[ok] transcripts with exons: {len(tx_map)}\n")
    annotated_stops_by_gene, all_annotated_stops, invalid_annotated_stops = build_annotated_stop_indexes(tx_map)
    sys.stderr.write(
        f"[ok] indexed explicit annotated CDS stops: {len(all_annotated_stops)} "
        f"across {len(annotated_stops_by_gene)} genes\n"
    )
    if invalid_annotated_stops:
        sys.stderr.write(
            f"[warn] skipped {invalid_annotated_stops} transcript stop_codon annotations "
            "that did not resolve to exactly 3 unique genomic positions\n"
        )

    # Common header for tORF TSV
    tsv_header = [
        "torf_id",
        "gene_id","gene_name","gene_biotype",
        "transcript_id","transcript_name","transcript_biotype","tx_tags","has_MANE_Select_tag",
        "chrom","strand",
        "tx_model","tx_complete_cds","tx_has_start_annot","tx_has_stop_annot","tx_stop_recovered","tx_cds_synthesized",
        "t_start","t_end","frame","aa_len","start_codon",
        "kozak_seq","kozak_class",
        "orf_biotype","overlap_cds_bp","overlap_utr5_bp","overlap_utr3_bp",
        "cds_first_t","cds_last_t",
        "is_annotated_cds","annotated_cds_id",
        "orf_has_upstream_start_annot","orf_shares_annotated_stop",
        "intorf_shares_annotated_stop_same_gene","intorf_shares_annotated_stop_any_gene",
        "intorf_annotated_stop_confounded",
        "chromStart1","chromEnd1","blockCount","blockSizes","genomic_block_starts1","genomic_block_ends1",
        "flank5_genomic_starts1","flank5_genomic_ends1",
        "flank3_genomic_starts1","flank3_genomic_ends1",
        "flank5_tx_len_requested","flank5_tx_len_actual","flank5_complete","flank5_piece_count",
        "flank3_tx_len_requested","flank3_tx_len_actual","flank3_complete","flank3_piece_count",
        "peptide_len","peptide"
    ]

    # ---------- Per-chromosome mode ----------
    if a.by_chrom:
        # group transcripts by chrom
        chrom2tids = {}
        for tid, tx in tx_map.items():
            if tx.chrom not in chrom_whitelist:
                continue
            chrom2tids.setdefault(tx.chrom, []).append(tid)

        # optional subset
        chrom_list = sorted(chrom2tids.keys())
        if a.chroms:
            want = set(a.chroms)
            chrom_list = [c for c in chrom_list if c in want]
            if not chrom_list:
                sys.exit("ERROR: --chroms filtered out all chromosomes present in GTF/FAI.")

        # prepare per-chrom outdir (optional)
        def per_chrom_prefix(chrom: str) -> str:
            import os
            if a.perchrom_outdir:
                os.makedirs(a.perchrom_outdir, exist_ok=True)
                base = os.path.basename(a.out_prefix)
                return os.path.join(a.perchrom_outdir, f"{base}.{chrom}")
            return f"{a.out_prefix}.{chrom}"

        # reset global merged outputs before appending / rebuilding
        for p in [f"{a.out_prefix}.gorf.tsv", f"{a.out_prefix}.gorf_members.tsv", f"{a.out_prefix}.gorf.faa",
                  f"{a.out_prefix}.torf.tsv", f"{a.out_prefix}.torf.faa", f"{a.out_prefix}.gorf_validation.tsv"]:
            if os.path.exists(p):
                os.remove(p)

        grand_tx = 0
        grand_orf = 0
        torf_shards: List[str] = []
        faa_shards: List[str] = []
        remaining_transcripts = a.max_transcripts if (a.max_transcripts and a.max_transcripts > 0) else None
        # Persistent registry for early cross-chromosome collision detection.
        # Final validation still runs as an independent second line of defense.
        seen_global_gorf_signatures: Dict[str, Tuple[str, str]] = {}

        for chrom in chrom_list:
            tids_this_chrom = chrom2tids[chrom]
            if remaining_transcripts is not None:
                if remaining_transcripts <= 0:
                    break
                tids_this_chrom = tids_this_chrom[:remaining_transcripts]
            payloads = [TxPayload(tid, tx_map[tid]) for tid in tids_this_chrom]
            if not payloads:
                continue
            if remaining_transcripts is not None:
                remaining_transcripts -= len(payloads)

            out_prefix_chrom = per_chrom_prefix(chrom)
            out_tsv = f"{out_prefix_chrom}.torf.tsv"
            out_faa = f"{out_prefix_chrom}.torf.faa"
            if a.skip_peptide_fasta and os.path.exists(out_faa):
                os.remove(out_faa)

            sys.stderr.write(f"[{chrom}] scanning transcripts: {len(payloads)}\n")

            n_tx = 0
            n_orf = 0
            w_faa = None
            with open(out_tsv, "w", encoding="utf-8") as w_tsv:
                if not a.skip_peptide_fasta:
                    w_faa = open(out_faa, "w", encoding="utf-8")
                w_tsv.write("\t".join(tsv_header) + "\n")
                try:
                    if a.workers <= 1:
                        for p in payloads:
                            n_tx += 1
                            rows, faa = process_tx((a.fa, p, start_set, a.min_aa, a.flank_nt))
                            annotate_intorf_altframe_stop_matches(rows, annotated_stops_by_gene, all_annotated_stops)
                            for row in rows:
                                w_tsv.write("\t".join(str(row[k]) for k in tsv_header) + "\n")
                            if w_faa is not None:
                                for tid, pep in faa:
                                    w_faa.write(f">{tid}\n{pep}\n")
                            n_orf += len(rows)
                    else:
                        with Pool(
                            processes=a.workers,
                            initializer=init_fasta_worker,
                            initargs=(a.fa,),
                        ) as pool:
                            it = pool.imap(
                                process_tx,
                                ((a.fa, p, start_set, a.min_aa, a.flank_nt) for p in payloads),
                                chunksize=max(1, a.mp_chunksize)
                            )
                            for rows, faa in it:
                                n_tx += 1
                                annotate_intorf_altframe_stop_matches(rows, annotated_stops_by_gene, all_annotated_stops)
                                for row in rows:
                                    w_tsv.write("\t".join(str(row[k]) for k in tsv_header) + "\n")
                                if w_faa is not None:
                                    for tid, pep in faa:
                                        w_faa.write(f">{tid}\n{pep}\n")
                                n_orf += len(rows)
                finally:
                    if w_faa is not None:
                        w_faa.close()

            grand_tx += n_tx
            grand_orf += n_orf
            sys.stderr.write(f"[{chrom}] [ok] transcripts processed: {n_tx}\n")
            sys.stderr.write(f"[{chrom}] [ok] ORFs reported: {n_orf}\n")
            sys.stderr.write(f"[{chrom}] [ok] wrote {out_tsv}\n")
            if not a.skip_peptide_fasta:
                sys.stderr.write(f"[{chrom}] [ok] wrote {out_faa}\n")

            # Track shards for the final merged global tORF/FAA outputs.
            torf_shards.append(out_tsv)
            if not a.skip_peptide_fasta:
                faa_shards.append(out_faa)

            if a.no_collapse:
                sys.stderr.write(f"[{chrom}] [skip] collapse disabled by --no-collapse\n")
                continue

            # -------- collapse on this chrom tORF, append gORFs globally, and update this chrom's tORF with gorf_id
            sys.stderr.write(f"[{chrom}] [collapse] reading tORFs...\n")
            rows = []
            with open(out_tsv, "r", encoding="utf-8") as f:
                r = csv.DictReader(f, delimiter="\t")
                hdr = r.fieldnames or []
                need = ["torf_id","chrom","strand","orf_biotype","chromStart1","chromEnd1","blockSizes","genomic_block_starts1","peptide"]
                for k in need:
                    if k not in hdr:
                        sys.exit(f"ERROR: required column '{k}' not found in {out_tsv}")
                for row in r:
                    if "peptide_len" not in row or row["peptide_len"] in ("", None):
                        row["peptide_len"] = str(len(row.get("peptide","")))
                    rows.append(row)

            g_list, members_map, faa_list, torf2gorf = collapse_torf_rows_igv(
                rows, a.collapse_merge_by
            )

            sys.stderr.write(f"[{chrom}] [collapse] appending gORF outputs to {a.out_prefix}.* ...\n")
            write_gorf_outputs_append(
                a.out_prefix,
                g_list,
                members_map,
                faa_list,
                seen_global_gorf_signatures=seen_global_gorf_signatures,
                source_label=chrom,
                write_faa=not a.skip_peptide_fasta,
            )

            # update per-chrom tORF with gorf_id; fail on any missing mapping
            update_torf_tsv_with_gorf_id(out_tsv, torf2gorf)

        sys.stderr.write(f"[all] [ok] transcripts processed (sum): {grand_tx}\n")
        sys.stderr.write(f"[all] [ok] ORFs reported (sum): {grand_orf}\n")

        # Build one canonical global tORF TSV from the already-backfilled shards.
        # Downstream tools should use this file, never a hand-concatenated intermediate.
        merged_torf = f"{a.out_prefix}.torf.tsv"
        merged_faa = f"{a.out_prefix}.torf.faa"
        merge_tsv_files_with_identical_header(torf_shards, merged_torf)
        if not a.skip_peptide_fasta:
            concatenate_text_files(faa_shards, merged_faa)
        if a.no_collapse:
            sys.stderr.write("[all] [skip] gORF validation disabled by --no-collapse\n")
            return
        validate_gorf_outputs(
            merged_torf,
            f"{a.out_prefix}.gorf.tsv",
            f"{a.out_prefix}.gorf_members.tsv",
            f"{a.out_prefix}.gorf_validation.tsv",
        )
        return

    # ---------- Original (single global run) ----------
    # Invalidate any success marker from an earlier run before replacing tORF
    # outputs.  A no-collapse run must also not leave stale gORF products under
    # the same prefix.
    stale_outputs = [f"{a.out_prefix}.gorf_validation.tsv"]
    if a.skip_peptide_fasta:
        stale_outputs.extend([
            f"{a.out_prefix}.torf.faa",
            f"{a.out_prefix}.gorf.faa",
        ])
    if a.no_collapse:
        stale_outputs.extend([
            f"{a.out_prefix}.gorf.tsv",
            f"{a.out_prefix}.gorf_members.tsv",
            f"{a.out_prefix}.gorf.faa",
        ])
    for path in stale_outputs:
        if os.path.exists(path):
            os.remove(path)
            sys.stderr.write(f"[init] removed stale output: {path}\n")

    tids = list(tx_map.keys())
    if a.max_transcripts and a.max_transcripts > 0:
        tids = tids[:a.max_transcripts]
    payloads = [TxPayload(tid, tx_map[tid]) for tid in tids]

    out_tsv = f"{a.out_prefix}.torf.tsv"
    out_faa = f"{a.out_prefix}.torf.faa"

    n_tx = 0; n_orf = 0
    w_faa = None
    with open(out_tsv, "w", encoding="utf-8") as w_tsv:
        if not a.skip_peptide_fasta:
            w_faa = open(out_faa, "w", encoding="utf-8")
        w_tsv.write("\t".join(tsv_header) + "\n")
        try:
            if a.workers <= 1:
                for p in payloads:
                    n_tx += 1
                    rows, faa = process_tx((a.fa, p, start_set, a.min_aa, a.flank_nt))
                    annotate_intorf_altframe_stop_matches(rows, annotated_stops_by_gene, all_annotated_stops)
                    for row in rows:
                        w_tsv.write("\t".join(str(row[k]) for k in tsv_header) + "\n")
                    if w_faa is not None:
                        for tid, pep in faa:
                            w_faa.write(f">{tid}\n{pep}\n")
                    n_orf += len(rows)
            else:
                with Pool(
                    processes=a.workers,
                    initializer=init_fasta_worker,
                    initargs=(a.fa,),
                ) as pool:
                    it = pool.imap(
                        process_tx,
                        ((a.fa, p, start_set, a.min_aa, a.flank_nt) for p in payloads),
                        chunksize=max(1, a.mp_chunksize)
                    )
                    for rows, faa in it:
                        n_tx += 1
                        annotate_intorf_altframe_stop_matches(rows, annotated_stops_by_gene, all_annotated_stops)
                        for row in rows:
                            w_tsv.write("\t".join(str(row[k]) for k in tsv_header) + "\n")
                        if w_faa is not None:
                            for tid, pep in faa:
                                w_faa.write(f">{tid}\n{pep}\n")
                        n_orf += len(rows)
        finally:
            if w_faa is not None:
                w_faa.close()

    sys.stderr.write(f"[ok] transcripts processed: {n_tx}\n")
    sys.stderr.write(f"[ok] ORFs reported: {n_orf}\n")
    sys.stderr.write(f"[ok] wrote {out_tsv}\n")
    if not a.skip_peptide_fasta:
        sys.stderr.write(f"[ok] wrote {out_faa}\n")

    if a.no_collapse:
        sys.stderr.write("[skip] collapse disabled by --no-collapse\n")
        return

    sys.stderr.write("[collapse] reading tORFs...\n")
    rows = []
    with open(out_tsv, "r", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        hdr = r.fieldnames or []
        need = ["torf_id","chrom","strand","orf_biotype","chromStart1","chromEnd1","blockSizes","genomic_block_starts1","peptide"]
        for k in need:
            if k not in hdr:
                sys.exit(f"ERROR: required column '{k}' not found in {out_tsv}")
        for row in r:
            if "peptide_len" not in row or row["peptide_len"] in ("", None):
                row["peptide_len"] = str(len(row.get("peptide","")))
            rows.append(row)

    g_list, members_map, faa_list, torf2gorf = collapse_torf_rows_igv(rows, a.collapse_merge_by)

    sys.stderr.write("[collapse] writing gORF outputs...\n")
    write_gorf_outputs(a.out_prefix, g_list, members_map, faa_list, write_faa=not a.skip_peptide_fasta)

    update_torf_tsv_with_gorf_id(out_tsv, torf2gorf)
    validate_gorf_outputs(
        out_tsv,
        f"{a.out_prefix}.gorf.tsv",
        f"{a.out_prefix}.gorf_members.tsv",
        f"{a.out_prefix}.gorf_validation.tsv",
    )
if __name__ == "__main__":
    main()
