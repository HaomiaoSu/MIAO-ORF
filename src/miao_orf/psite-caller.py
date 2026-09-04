#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#suhaomiao@gmail.com

"""
Multiprocessing P-site caller (Ribo-seq, MD required)

- Ribo-TISH quality/offset import with corrected frame-0 screening
- 5' anchor offset (Ribo-TISH definition) + fast CIGAR walk
- Per-contig multiprocessing + BGZF I/O threads
- Emits 1-nt P-site BAM; BED is derived from BAM for exact coordinate match
- Optional strand-split bedGraph from BAM

CLI example (one line):
python psite-caller.py --bam sample.mapped.bam --offsets sample.mapped.para.py --ribotish-quality sample.mapped_qual.txt --out-prefix sample --workers 8 --merge
python psite-caller.py --bam sample.mapped.bam --length-offsets 28:12 29:12 30:12 --out-prefix sample --workers 8 --merge
"""

import argparse
import json
import os
import re
import sys
import shutil
import subprocess
import time
from collections import defaultdict
from multiprocessing import Pool, get_context

import pysam

try:
    from .ribotish_offsets import (
        load_ribotish_offset_dicts,
        parse_length_offset_specs,
        select_explicit_length_offsets,
        select_ribotish_offsets,
        write_selection_tsv,
    )
except ImportError:
    from ribotish_offsets import (
        load_ribotish_offset_dicts,
        parse_length_offset_specs,
        select_explicit_length_offsets,
        select_ribotish_offsets,
        write_selection_tsv,
    )

PROGRAM = "miao-orf-psite"
VERSION = "1.0.0"


# ------------------------ CLI ------------------------
def parse_args():
    ap = argparse.ArgumentParser(
        prog=PROGRAM,
        description="MIAO P-site caller (MD required), multiprocessing, BED-from-BAM, bedGraph.",
    )
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    ap.add_argument("--bam", required=True, help="Input BAM (aligned; all aligned reads must carry MD tag)")
    ap.add_argument("--offsets",
                    help="Automatic mode: Ribo-TISH mapped.para.py containing offdict")
    ap.add_argument(
        "--length-offsets", nargs="+", metavar="LENGTH:OFFSET",
        help=(
            "Traditional mode: explicit read-length/P-site-offset pairs, for example "
            "--length-offsets 28:12 29:12 30:12"
        ),
    )
    ap.add_argument(
        "--ribotish-quality", "--quality", dest="ribotish_quality",
        help=(
            "Ribo-TISH mapped_qual.txt. By default, select lengths from the 5'-matched "
            "group whose offset-corrected frame-0 proportion is strictly above 2/3, "
            "then retain the read-richest contiguous passing-length block. "
            "Required unless --keep-lengths is supplied as an explicit legacy override."
        ),
    )
    ap.add_argument("--out-prefix", required=True, help="Output prefix; parts go into <prefix>.parts/")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument(
        "--min-frame0-proportion", type=float, default=2 / 3,
        help="Automatic mode: require corrected frame-0 proportion above this value",
    )
    ap.add_argument(
        "--length-selection-policy",
        choices=("dominant_contiguous", "all_passing"),
        default="dominant_contiguous",
        help="Automatic mode: keep the read-richest contiguous block or every passing length",
    )
    ap.add_argument(
        "--keep-lengths", nargs="+", type=int,
        help=(
            "Optional read-length whitelist applied after quality screening. If "
            "--ribotish-quality is omitted, this becomes an explicit legacy override."
        ),
    )
    ap.add_argument("--require-unique", action=argparse.BooleanOptionalAction, default=True,
                    help="Require NH==1 if present (default: True; use --no-require-unique to allow multimappers)")
    ap.add_argument("--workers", type=int, default=4, help="CPU processes")
    ap.add_argument("--bgzf-threads", type=int, default=2, help="htslib BGZF I/O threads per file handle")
    ap.add_argument("--no-bed", action="store_true", help="Do not emit BED (default: emit BED derived from BAM)")
    ap.add_argument("--bed-flush", type=int, default=200000, help="Flush BED every N records; -1 = only at end")
    ap.add_argument("--merge", action=argparse.BooleanOptionalAction, default=True,
                    help="Merge per-contig BAMs (samtools merge+sort+index), then derive BED/bedGraph (default: True; use --no-merge to disable)")
    ap.add_argument("--no-bedgraph", action="store_true", help="Disable bedGraph output")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


# --------------------- Utilities ---------------------
def log(msg: str):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def sanitize(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._+-]", "_", name)


def load_offsets_from_parapy(path: str) -> dict:
    """Backward-compatible helper returning all main Ribo-TISH offsets."""
    try:
        selected, _mismatch = load_ribotish_offset_dicts(path)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    return selected


def find_q5_query_index(r: pysam.AlignedSegment):
    c = r.cigartuples or []
    if not r.is_reverse:
        q = 0
        for op, ln in c:
            if op == 5:  # H
                continue
            if op in (2, 3):  # D/N
                continue
            if op in (4, 1):  # S/I
                q += ln
                continue
            if op in (0, 7, 8):  # M/= /X
                return q
        return None
    else:
        q = r.query_length
        for op, ln in reversed(c):
            if op == 5:
                continue
            if op in (2, 3):
                continue
            if op in (4, 1):
                q -= ln
                continue
            if op in (0, 7, 8):
                return q - 1
        return None


def refpos_from_qidx_fast(r: pysam.AlignedSegment, qidx: int):
    if qidx is None or qidx < 0 or qidx >= r.query_length:
        return None
    q = 0
    ref = r.reference_start
    for op, ln in (r.cigartuples or []):
        if op == 5:  # H
            continue
        if op == 4:  # S
            q += ln
            continue
        if op in (0, 7, 8):  # M/= /X
            if qidx < q + ln:
                return ref + (qidx - q)
            q += ln
            ref += ln
            continue
        if op == 1:  # I
            if qidx < q + ln:
                return None
            q += ln
            continue
        if op in (2, 3):  # D/N
            ref += ln
            continue
    return None


def fiveprime_is_aligned(r: pysam.AlignedSegment) -> bool:
    c = r.cigartuples or []
    if not r.is_reverse:
        for op, _ in c:
            if op == 5:
                continue
            if op in (2, 3):
                continue
            if op in (4, 1):
                return False
            if op in (0, 7, 8):
                return True
        return False
    else:
        for op, _ in reversed(c):
            if op == 5:
                continue
            if op in (2, 3):
                continue
            if op in (4, 1):
                return False
            if op in (0, 7, 8):
                return True
        return False


_MD_TOKEN = re.compile(r"(\d+|\^[A-Z]+|[A-Z])")


def fiveprime_match_by_eqx(r: pysam.AlignedSegment):
    c = r.cigartuples or []
    if not r.is_reverse:
        for op, _ in c:
            if op in (5,):
                continue
            if op in (2, 3):
                continue
            if op in (4, 1):
                return False
            if op == 7:
                return True
            if op == 8:
                return False
            if op == 0:
                return None
        return None
    else:
        for op, _ in reversed(c):
            if op in (5,):
                continue
            if op in (2, 3):
                continue
            if op in (4, 1):
                return False
            if op == 7:
                return True
            if op == 8:
                return False
            if op == 0:
                return None
        return None


def fiveprime_match_by_md_or_die(r: pysam.AlignedSegment) -> bool:
    try:
        md = r.get_tag("MD")
    except KeyError:
        raise RuntimeError("Missing MD tag. Fix with: samtools calmd -b <in.bam> <ref.fa> > out.bam")
    toks = _MD_TOKEN.findall(md)
    if not toks:
        return False
    if not r.is_reverse:
        i = 0
        while i < len(toks):
            t = toks[i]
            if t.startswith("^"):
                i += 1
                continue
            if t.isdigit():
                if int(t) >= 1:
                    return True
                i += 1
                continue
            return False
        return False
    else:
        i = len(toks) - 1
        while i >= 0:
            t = toks[i]
            if t.startswith("^"):
                i -= 1
                continue
            if t.isdigit():
                if int(t) >= 1:
                    return True
                i -= 1
                continue
            return False
        return False


def base_at_psite(r: pysam.AlignedSegment, qidx: int) -> str:
    if qidx is None or qidx < 0 or qidx >= r.query_length:
        return "N"
    b = r.query_sequence[qidx]
    return (b or "N").upper()


# --------------------- Worker ---------------------
def worker_task(args):
    (bam_path, contig, off_map, allowed_lengths, out_dir, out_prefix,
     min_mapq, require_unique, bgzf_threads, verbose) = args

    tag = sanitize(contig)
    out_stem = os.path.basename(os.path.normpath(out_prefix))
    if not out_stem:
        raise ValueError(f"Invalid output prefix: {out_prefix!r}")
    bam_out_path = os.path.join(out_dir, f"{out_stem}.psite.{tag}.bam")

    in_bam = pysam.AlignmentFile(bam_path, "rb", threads=max(1, bgzf_threads))
    header = in_bam.header
    out_bam = pysam.AlignmentFile(bam_out_path, "wb", header=header, threads=max(1, bgzf_threads))

    stats = dict(kept=0, skipped_len=0, skipped_mapq=0, skipped_unique=0,
                 skipped_5p_not_aligned=0, skipped_5p_mismatch=0, skipped_qmap=0)

    start = time.time()
    log(f"[start] {contig}")
    try:
        for r in in_bam.fetch(contig):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if r.mapping_quality < min_mapq:
                stats["skipped_mapq"] += 1
                continue
            if require_unique:
                try:
                    if r.get_tag("NH") != 1:
                        stats["skipped_unique"] += 1
                        continue
                except KeyError:
                    pass

            L = r.query_length
            if L not in allowed_lengths:
                stats["skipped_len"] += 1
                continue
            if not fiveprime_is_aligned(r):
                stats["skipped_5p_not_aligned"] += 1
                continue

            ex = fiveprime_match_by_eqx(r)
            ok = fiveprime_match_by_md_or_die(r) if ex is None else ex
            if not ok:
                stats["skipped_5p_mismatch"] += 1
                continue

            off = off_map.get(L)
            if off is None:
                stats["skipped_len"] += 1
                continue

            q5 = find_q5_query_index(r)
            if q5 is None:
                stats["skipped_5p_not_aligned"] += 1
                continue
            qidx = q5 + off if not r.is_reverse else q5 - off

            ref_pos = refpos_from_qidx_fast(r, qidx)
            if ref_pos is None:
                stats["skipped_qmap"] += 1
                continue

            new = pysam.AlignedSegment(out_bam.header)
            new.query_name = r.query_name
            new.flag = 16 if r.is_reverse else 0
            new.reference_id = r.reference_id
            new.reference_start = ref_pos
            new.mapping_quality = r.mapping_quality
            new.cigar = ((0, 1),)  # 1M
            new.next_reference_id = -1
            new.next_reference_start = -1
            new.template_length = 0
            new.query_sequence = base_at_psite(r, qidx)
            new.query_qualities = pysam.qualitystring_to_array("I")
            new.set_tag("OL", L)
            new.set_tag("AL", r.query_alignment_length)
            new.set_tag("PO", off)
            new.set_tag("RL", L)
            try:
                new.set_tag("NH", r.get_tag("NH"))
            except KeyError:
                pass

            out_bam.write(new)
            stats["kept"] += 1

    except RuntimeError as e:
        out_bam.close()
        in_bam.close()
        return {"contig": contig, "error": str(e)}

    out_bam.close()
    in_bam.close()

    dur = time.time() - start
    if verbose:
        log(f"[done] {contig} kept={stats['kept']} len={stats['skipped_len']} mapq={stats['skipped_mapq']} unique={stats['skipped_unique']} 5pNotAligned={stats['skipped_5p_not_aligned']} 5pMismatch={stats['skipped_5p_mismatch']} qmapFail={stats['skipped_qmap']} ({dur:.1f}s)")
    else:
        log(f"[done] {contig} kept={stats['kept']} ({dur:.1f}s)")

    return {"contig": contig, "stats": stats, "bam": bam_out_path}


# --------------- Derive BED / bedGraph ---------------
def write_bed_from_bam(bam_path: str, bed_path: str, buffer_n: int = 200000):
    bam = pysam.AlignmentFile(bam_path, "rb")
    with open(bed_path, "w") as w:
        buf = []
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped:
                continue
            chrom = bam.get_reference_name(r.reference_id)
            start = r.reference_start
            end = start + 1
            strand = "-" if r.is_reverse else "+"
            buf.append(f"{chrom}\t{start}\t{end}\t.\t1\t{strand}\n")
            if buffer_n != -1 and len(buf) >= buffer_n:
                w.writelines(buf)
                buf.clear()
        if buf:
            w.writelines(buf)
    bam.close()


def make_bedgraph_from_bam(bam_path: str, out_prefix: str):
    plus = f"{out_prefix}.plus.bedGraph"
    minus = f"{out_prefix}.minus.bedGraph"
    try:
        with open(plus, "w") as f1:
            subprocess.run(["bedtools", "genomecov", "-bg", "-ibam", bam_path, "-strand", "+"],
                           check=True, stdout=f1)
        with open(minus, "w") as f2:
            subprocess.run(["bedtools", "genomecov", "-bg", "-ibam", bam_path, "-strand", "-"],
                           check=True, stdout=f2)
        log(f"[ok] bedGraph: {plus} , {minus}")
    except FileNotFoundError:
        log("[warn] bedtools not found; skip bedGraph")
    except subprocess.CalledProcessError as e:
        log(f"[warn] bedGraph failed: {e}")


# ------------------------ Main ------------------------
def main():
    a = parse_args()

    try:
        if a.length_offsets is not None:
            if a.offsets or a.ribotish_quality or a.keep_lengths:
                raise ValueError(
                    "--length-offsets cannot be combined with --offsets, "
                    "--ribotish-quality or --keep-lengths"
                )
            off_map, selection_rows = select_explicit_length_offsets(
                parse_length_offset_specs(a.length_offsets)
            )
            mode = "explicit length:offset"
        else:
            if not a.offsets:
                raise ValueError(
                    "choose traditional --length-offsets mode or provide --offsets "
                    "for Ribo-TISH mode"
                )
            off_map, selection_rows = select_ribotish_offsets(
                a.offsets,
                quality_path=a.ribotish_quality,
                min_frame0_proportion=a.min_frame0_proportion,
                keep_lengths=a.keep_lengths,
                length_selection_policy=a.length_selection_policy,
            )
            mode = "Ribo-TISH automatic" if a.ribotish_quality else "legacy para.py + lengths"
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    if not off_map:
        sys.exit("ERROR: No read length passed the Ribo-TISH offset quality contract.")
    allowed = set(off_map)
    selection_path = f"{a.out_prefix}.offset_selection.tsv"
    write_selection_tsv(selection_path, selection_rows)
    selected_text = ", ".join(f"{length}:{off_map[length]}" for length in sorted(off_map))
    log(f"[offsets] mode={mode} selected={selected_text}")
    log(f"[offsets] audit={selection_path}")

    try:
        h = pysam.AlignmentFile(a.bam, "rb")
    except ValueError as e:
        sys.exit(f"ERROR: cannot open BAM: {e}")
    contigs = [rname for rname, rlen in zip(h.references, h.lengths) if rlen and rlen > 0]
    h.close()
    if not contigs:
        sys.exit("ERROR: No contigs in header.")

    parts_dir = f"{a.out_prefix}.parts"
    os.makedirs(parts_dir, exist_ok=True)
    log(f"[init] contigs={len(contigs)} workers={a.workers} bgzf-threads={a.bgzf_threads} out={parts_dir}")

    tasks = []
    for ctg in contigs:
        tasks.append((a.bam, ctg, off_map, allowed, parts_dir, a.out_prefix,
                      a.min_mapq, a.require_unique, max(1, a.bgzf_threads), a.verbose))

    agg = defaultdict(int)
    part_bams = []
    ctx = get_context("spawn")
    with ctx.Pool(processes=max(1, a.workers)) as pool:
        for res in pool.imap_unordered(worker_task, tasks):
            if "error" in res:
                pool.terminate()
                pool.join()
                sys.exit(f"ERROR (contig {res.get('contig')}): {res['error']}")
            for k, v in res["stats"].items():
                agg[k] += v
            part_bams.append(res["bam"])
            log(f"[progress] {len(part_bams)}/{len(tasks)} parts done")

    log(f"[total] kept={agg['kept']} skip_len={agg['skipped_len']} skip_mapq={agg['skipped_mapq']} skip_unique={agg['skipped_unique']} skip_5pNotAligned={agg['skipped_5p_not_aligned']} skip_5pMismatch={agg['skipped_5p_mismatch']} qidx_map_fail={agg['skipped_qmap']}")

    if a.merge:
        merged_unsorted = f"{a.out_prefix}.psite.unsorted.bam"
        merged_bam = f"{a.out_prefix}.psite.bam"
        merge_cmd = ["samtools", "merge", "-f", "-@", str(max(1, a.workers)), merged_unsorted] + part_bams
        sort_cmd = ["samtools", "sort", "-@", str(max(1, a.workers)), "-o", merged_bam, merged_unsorted]
        index_cmd = ["samtools", "index", merged_bam]
        try:
            subprocess.run(merge_cmd, check=True)
            subprocess.run(sort_cmd, check=True)
            subprocess.run(index_cmd, check=True)
        except FileNotFoundError:
            sys.exit("ERROR: 'samtools' not found. Install or run without --merge.")
        except subprocess.CalledProcessError as e:
            sys.exit(f"ERROR during samtools merge/sort/index: {e}")
        try:
            os.remove(merged_unsorted)
        except OSError:
            pass
        log(f"[ok] merged BAM: {merged_bam} (+.bai)")

        if not a.no_bed:
            merged_bed = f"{a.out_prefix}.psite.bed"
            write_bed_from_bam(merged_bam, merged_bed, buffer_n=a.bed_flush)
            try:
                subprocess.run(["bash", "-lc", f"sort -k1,1 -k2,2n {merged_bed} -o {merged_bed}"], check=True)
            except Exception as e:
                log(f"[warn] BED sort failed: {e}")
            log(f"[ok] merged BED: {merged_bed}")

        if not a.no_bedgraph:
            make_bedgraph_from_bam(merged_bam, f"{a.out_prefix}.psite")

        try:
            shutil.rmtree(parts_dir)
            log(f"[ok] removed parts dir: {parts_dir}")
        except Exception as e:
            log(f"[warn] could not remove parts dir {parts_dir}: {e}")
    else:
        if not a.no_bed:
            for pb in part_bams:
                bedp = pb.replace(".bam", ".bed")
                write_bed_from_bam(pb, bedp, buffer_n=a.bed_flush)
            log("[ok] per-contig BEDs generated from per-contig BAMs")
        if not a.no_bedgraph:
            for pb in part_bams:
                prefix = pb[:-4]
                make_bedgraph_from_bam(pb, prefix)
        log(f"[ok] per-contig outputs in: {parts_dir}")


if __name__ == "__main__":
    main()
