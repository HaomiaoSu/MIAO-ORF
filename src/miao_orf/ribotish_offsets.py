#!/usr/bin/env python3
"""Resolve explicit or Ribo-TISH-derived read lengths and P-site offsets.

Ribo-TISH ``quality`` writes five Python-literal dictionaries for reads whose
5' nucleotide matches the reference, followed (when enabled) by five
dictionaries for reads with a mismatch at that nucleotide.  The fourth
dictionary in each group stores the three frame counts.  Those counts use
Ribo-TISH's fixed QC offset of 12 nt, so they must be rotated by the
length-specific offset in ``mapped.para.py`` before testing corrected frame 0.
"""

from __future__ import annotations

import argparse
import ast
import csv
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple


PROGRAM = "miao-orf-offsets"
VERSION = "1.0.0"
RIBOTISH_QC_OFFSET = 12


@dataclass(frozen=True)
class OffsetSelection:
    selection_mode: str
    quality_group: str
    ribotish_qc_offset: Optional[int]
    min_frame0_proportion: Optional[float]
    length_selection_policy: str
    length: int
    offset: Optional[int]
    matched_reads: Optional[int]
    raw_frame0: Optional[int]
    raw_frame1: Optional[int]
    raw_frame2: Optional[int]
    corrected_frame0: Optional[int]
    corrected_frame1: Optional[int]
    corrected_frame2: Optional[int]
    corrected_frame0_proportion: Optional[float]
    selected: bool
    reason: str


def _literal_dict(text: str, label: str) -> dict:
    try:
        value = ast.literal_eval(text)
    except (SyntaxError, ValueError) as exc:
        raise ValueError(f"{label} is not a valid Python literal dictionary") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a dictionary")
    return value


def _int_dict(value: Mapping, label: str) -> Dict[int, int]:
    result: Dict[int, int] = {}
    for key, item in value.items():
        if not isinstance(key, int) or isinstance(key, bool):
            continue
        if not isinstance(item, int) or isinstance(item, bool) or item < 0:
            raise ValueError(f"{label}[{key!r}] must be a non-negative integer")
        result[int(key)] = int(item)
    if not result:
        raise ValueError(f"{label} contains no integer-keyed values")
    return result


def _frame_dict(value: Mapping, label: str) -> Dict[int, Tuple[int, int, int]]:
    result: Dict[int, Tuple[int, int, int]] = {}
    for key, item in value.items():
        if not isinstance(key, int) or isinstance(key, bool):
            continue
        if not isinstance(item, (list, tuple)) or len(item) != 3:
            raise ValueError(f"{label}[{key!r}] must contain exactly three frame counts")
        if any(not isinstance(count, int) or isinstance(count, bool) or count < 0 for count in item):
            raise ValueError(f"{label}[{key!r}] frame counts must be non-negative integers")
        result[int(key)] = (int(item[0]), int(item[1]), int(item[2]))
    if not result:
        raise ValueError(f"{label} contains no integer-keyed frame counts")
    return result


def load_ribotish_offset_dicts(path: str | Path) -> Tuple[Dict[int, int], Dict[int, int]]:
    """Return the 5'-match and 5'-mismatch offset dictionaries from para.py."""
    source = Path(path).read_text(encoding="utf-8")
    try:
        tree = ast.parse(source, filename=str(path), mode="exec")
    except SyntaxError as exc:
        raise ValueError(f"invalid Ribo-TISH parameter file: {path}") from exc

    offdict = None
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == "offdict" for target in node.targets):
            try:
                offdict = ast.literal_eval(node.value)
            except (SyntaxError, ValueError) as exc:
                raise ValueError("offdict must be a literal dictionary") from exc
            break
    if not isinstance(offdict, dict):
        raise ValueError(f"could not find a literal offdict assignment in {path}")

    main = _int_dict(offdict, "offdict")
    m0_value = offdict.get("m0", {})
    if m0_value and not isinstance(m0_value, dict):
        raise ValueError("offdict['m0'] must be a dictionary")
    mismatch = _int_dict(m0_value, "offdict['m0']") if m0_value else {}
    return main, mismatch


def load_ribotish_quality(
    path: str | Path,
) -> Tuple[Dict[int, int], Dict[int, Tuple[int, int, int]], Dict[int, int], Dict[int, Tuple[int, int, int]]]:
    """Read Ribo-TISH quality output without executing its Python literals."""
    lines = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) not in (5, 10):
        raise ValueError(
            f"Ribo-TISH quality file must contain 5 or 10 non-empty dictionary lines; found {len(lines)}"
        )

    values = [_literal_dict(line, f"quality line {index + 1}") for index, line in enumerate(lines)]
    matched_lengths = _int_dict(values[0], "matched length counts")
    matched_frames = _frame_dict(values[3], "matched frame counts")
    if len(values) == 10:
        mismatch_lengths = _int_dict(values[5], "5'-mismatch length counts")
        mismatch_frames = _frame_dict(values[8], "5'-mismatch frame counts")
    else:
        mismatch_lengths, mismatch_frames = {}, {}
    return matched_lengths, matched_frames, mismatch_lengths, mismatch_frames


def corrected_frame_counts(
    raw_counts: Sequence[int], offset: int, qc_offset: int = RIBOTISH_QC_OFFSET
) -> Tuple[int, int, int]:
    """Rotate QC frames from the fixed 12-nt anchor to a selected P-site offset."""
    if len(raw_counts) != 3:
        raise ValueError("raw_counts must contain exactly three values")
    corrected = [0, 0, 0]
    for raw_frame, count in enumerate(raw_counts):
        corrected_frame = (raw_frame + offset - qc_offset) % 3
        corrected[corrected_frame] += int(count)
    return corrected[0], corrected[1], corrected[2]


def dominant_contiguous_block(
    lengths: Iterable[int], weights: Mapping[int, int]
) -> set[int]:
    """Return the consecutive passing-length block with the greatest total weight."""
    ordered = sorted(set(int(length) for length in lengths))
    if not ordered:
        return set()
    blocks: list[list[int]] = [[ordered[0]]]
    for length in ordered[1:]:
        if length == blocks[-1][-1] + 1:
            blocks[-1].append(length)
        else:
            blocks.append([length])
    best = max(
        blocks,
        key=lambda block: (
            sum(int(weights.get(length, 0)) for length in block),
            len(block),
            -block[0],
        ),
    )
    return set(best)


def parse_length_offset_specs(values: Iterable[str]) -> Dict[int, int]:
    """Parse one or more comma/space-separated ``LENGTH:OFFSET`` specifications."""
    result: Dict[int, int] = {}
    for value in values:
        for item in str(value).split(","):
            token = item.strip()
            if not token:
                continue
            parts = token.split(":")
            if len(parts) != 2:
                raise ValueError(
                    f"invalid length:offset value {token!r}; expected, for example, 28:12"
                )
            try:
                length, offset = (int(part.strip()) for part in parts)
            except ValueError as exc:
                raise ValueError(
                    f"invalid length:offset value {token!r}; both fields must be integers"
                ) from exc
            if length <= 0:
                raise ValueError(f"read length must be positive: {length}")
            if not 0 <= offset < length:
                raise ValueError(
                    f"offset for read length {length} must satisfy 0 <= offset < {length}"
                )
            if length in result:
                raise ValueError(f"duplicate read length in --length-offsets: {length}")
            result[length] = offset
    if not result:
        raise ValueError("--length-offsets must contain at least one LENGTH:OFFSET pair")
    return result


def select_explicit_length_offsets(
    length_offsets: Mapping[int, int],
) -> Tuple[Dict[int, int], list[OffsetSelection]]:
    """Validate and audit a traditional user-supplied length/offset mapping."""
    normalized = parse_length_offset_specs(
        f"{length}:{offset}" for length, offset in length_offsets.items()
    )
    rows = [
        OffsetSelection(
            selection_mode="explicit_length_offsets",
            quality_group="",
            ribotish_qc_offset=None,
            min_frame0_proportion=None,
            length_selection_policy="explicit_length_offsets",
            length=length,
            offset=normalized[length],
            matched_reads=None,
            raw_frame0=None,
            raw_frame1=None,
            raw_frame2=None,
            corrected_frame0=None,
            corrected_frame1=None,
            corrected_frame2=None,
            corrected_frame0_proportion=None,
            selected=True,
            reason="selected_explicit_length_offset",
        )
        for length in sorted(normalized)
    ]
    return dict(sorted(normalized.items())), rows


def select_ribotish_offsets(
    offsets_path: str | Path,
    quality_path: Optional[str | Path] = None,
    min_frame0_proportion: float = 2 / 3,
    keep_lengths: Optional[Iterable[int]] = None,
    length_selection_policy: str = "dominant_contiguous",
) -> Tuple[Dict[int, int], list[OffsetSelection]]:
    """Select main-group offsets using corrected frame-0 periodicity.

    ``quality_path`` is the default evidence-based mode.  ``keep_lengths`` can
    additionally restrict that selection.  Supplying only ``keep_lengths`` is
    retained as an explicit legacy override for frozen analyses.
    """
    if not 0 <= min_frame0_proportion < 1:
        raise ValueError("min_frame0_proportion must be in [0, 1)")
    if length_selection_policy not in {"dominant_contiguous", "all_passing"}:
        raise ValueError(
            "length_selection_policy must be 'dominant_contiguous' or 'all_passing'"
        )
    main_offsets, _mismatch_offsets = load_ribotish_offset_dicts(offsets_path)
    whitelist = None if keep_lengths is None else {int(length) for length in keep_lengths}
    if whitelist is not None and (not whitelist or any(length <= 0 for length in whitelist)):
        raise ValueError("keep_lengths must contain positive read lengths")

    if quality_path is None:
        if whitelist is None:
            raise ValueError("provide --ribotish-quality or an explicit --keep-lengths override")
        rows = []
        selected: Dict[int, int] = {}
        for length in sorted(set(main_offsets) | whitelist):
            offset = main_offsets.get(length)
            use = length in whitelist and offset is not None and 0 <= offset < length
            if length not in whitelist:
                reason = "not_in_keep_lengths"
            elif offset is None:
                reason = "missing_main_offset"
            elif not 0 <= offset < length:
                reason = "invalid_offset"
            else:
                reason = "selected_explicit_override"
                selected[length] = offset
            rows.append(
                OffsetSelection(
                    selection_mode="explicit_keep_lengths",
                    quality_group="",
                    ribotish_qc_offset=None,
                    min_frame0_proportion=None,
                    length_selection_policy="explicit_keep_lengths",
                    length=length,
                    offset=offset,
                    matched_reads=None,
                    raw_frame0=None,
                    raw_frame1=None,
                    raw_frame2=None,
                    corrected_frame0=None,
                    corrected_frame1=None,
                    corrected_frame2=None,
                    corrected_frame0_proportion=None,
                    selected=use,
                    reason=reason,
                )
            )
        return selected, rows

    matched_lengths, matched_frames, _mismatch_lengths, _mismatch_frames = load_ribotish_quality(quality_path)
    preliminary = []
    passing_lengths = set()
    for length in sorted(set(main_offsets) | set(matched_frames)):
        offset = main_offsets.get(length)
        raw = matched_frames.get(length)
        corrected = None
        proportion = None
        if offset is not None and raw is not None and 0 <= offset < length:
            corrected = corrected_frame_counts(raw, offset)
            total = sum(corrected)
            proportion = corrected[0] / total if total else None

        if offset is None:
            reason = "missing_main_offset"
        elif not 0 <= offset < length:
            reason = "invalid_offset"
        elif raw is None:
            reason = "missing_matched_frame_counts"
        elif proportion is None:
            reason = "no_matched_frame_reads"
        elif proportion <= min_frame0_proportion:
            reason = "frame0_not_above_threshold"
        elif whitelist is not None and length not in whitelist:
            reason = "not_in_keep_lengths"
        else:
            reason = "passes_frame0_threshold"
            passing_lengths.add(length)

        preliminary.append((length, offset, raw, corrected, proportion, reason))

    if length_selection_policy == "dominant_contiguous":
        retained_lengths = dominant_contiguous_block(
            passing_lengths,
            {
                length: matched_lengths.get(length, sum(matched_frames.get(length, ())))
                for length in passing_lengths
            },
        )
    else:
        retained_lengths = set(passing_lengths)

    selected = {}
    rows = []
    for length, offset, raw, corrected, proportion, preliminary_reason in preliminary:
        if preliminary_reason != "passes_frame0_threshold":
            reason = preliminary_reason
        elif length not in retained_lengths:
            reason = "outside_dominant_contiguous_block"
        else:
            reason = "selected_quality"
            assert offset is not None
            selected[length] = offset

        rows.append(
            OffsetSelection(
                selection_mode="ribotish_quality",
                quality_group="5prime_match",
                ribotish_qc_offset=RIBOTISH_QC_OFFSET,
                min_frame0_proportion=min_frame0_proportion,
                length_selection_policy=length_selection_policy,
                length=length,
                offset=offset,
                matched_reads=matched_lengths.get(length),
                raw_frame0=raw[0] if raw else None,
                raw_frame1=raw[1] if raw else None,
                raw_frame2=raw[2] if raw else None,
                corrected_frame0=corrected[0] if corrected else None,
                corrected_frame1=corrected[1] if corrected else None,
                corrected_frame2=corrected[2] if corrected else None,
                corrected_frame0_proportion=proportion,
                selected=reason == "selected_quality",
                reason=reason,
            )
        )
    return selected, rows


def write_selection_tsv(path: str | Path, rows: Sequence[OffsetSelection]) -> None:
    fieldnames = list(asdict(rows[0]).keys()) if rows else list(OffsetSelection.__dataclass_fields__)
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            values = asdict(row)
            if values["corrected_frame0_proportion"] is not None:
                values["corrected_frame0_proportion"] = f"{values['corrected_frame0_proportion']:.10g}"
            values["selected"] = int(values["selected"])
            writer.writerow({key: "" if value is None else value for key, value in values.items()})


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog=PROGRAM,
        description=(
            "Resolve traditional explicit length:offset pairs or automatically "
            "select Ribo-TISH offsets using corrected frame-0 periodicity."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--offsets", help="Ribo-TISH mapped.para.py for automatic mode")
    parser.add_argument("--ribotish-quality", help="Ribo-TISH mapped_qual.txt for automatic mode")
    parser.add_argument(
        "--length-offsets", nargs="+", metavar="LENGTH:OFFSET",
        help="Traditional mode: explicit pairs such as 28:12 29:12 30:12",
    )
    parser.add_argument("--min-frame0-proportion", type=float, default=2 / 3)
    parser.add_argument(
        "--length-selection-policy",
        choices=("dominant_contiguous", "all_passing"),
        default="dominant_contiguous",
    )
    parser.add_argument("--keep-lengths", nargs="+", type=int)
    parser.add_argument("--out", help="Output audit TSV; default: print selected length:offset pairs")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    try:
        if args.length_offsets is not None:
            if args.offsets or args.ribotish_quality or args.keep_lengths:
                raise ValueError(
                    "--length-offsets cannot be combined with --offsets, "
                    "--ribotish-quality or --keep-lengths"
                )
            selected, rows = select_explicit_length_offsets(
                parse_length_offset_specs(args.length_offsets)
            )
        else:
            if not args.offsets:
                raise ValueError(
                    "choose traditional --length-offsets mode or provide --offsets "
                    "for Ribo-TISH mode"
                )
            selected, rows = select_ribotish_offsets(
                args.offsets,
                args.ribotish_quality,
                args.min_frame0_proportion,
                args.keep_lengths,
                args.length_selection_policy,
            )
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    if not selected:
        raise SystemExit("ERROR: no read length passed the Ribo-TISH offset quality contract")
    if args.out:
        write_selection_tsv(args.out, rows)
    print(",".join(f"{length}:{selected[length]}" for length in sorted(selected)))


if __name__ == "__main__":
    main()
