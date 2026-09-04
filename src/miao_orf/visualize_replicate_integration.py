#!/usr/bin/env python3
"""Visualize evidence-preserving MIAO replicate integration outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Mapping, NamedTuple, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
from scipy.stats import spearmanr


PROGRAM = "miao-orf-visualize-replicates"
VERSION = "1.0.0"
PFPKM_REQUIRED_COLUMNS = {
    "gorf_id", "overlap_type", "intorf_pFPKM", "abundance_status",
}

COLORS = {
    "total": "#4C78A8",
    "significant": "#F58518",
    "primary": "#54A24B",
    "neutral": "#9D9D9D",
    "absent": "#E6E6E6",
    "unavailable": "#B8B8B8",
}


def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        return fields, [dict(row) for row in reader]


def parse_assignment(value: str, label: str) -> Tuple[str, Path]:
    if "=" not in value:
        raise ValueError(f"{label} must use REPLICATE_ID=PATH")
    replicate, assigned = value.split("=", 1)
    replicate = replicate.strip()
    assigned = assigned.strip()
    if not replicate or not assigned:
        raise ValueError(f"{label} must use non-empty REPLICATE_ID=PATH")
    return replicate, Path(assigned)


def read_abundance_tables(
    assignments: Sequence[str], replicates: Sequence[str]
) -> Tuple[Dict[str, Dict[str, Dict[str, str]]], Dict[str, Path]]:
    """Read optional pFPKM tables without changing DM integration semantics."""
    values: Dict[str, Dict[str, Dict[str, str]]] = {}
    paths: Dict[str, Path] = {}
    expected = set(replicates)
    for assignment in assignments:
        replicate, path = parse_assignment(assignment, "--abundance")
        if replicate not in expected:
            raise ValueError(f"abundance supplied for unexpected replicate: {replicate}")
        if replicate in values:
            raise ValueError(f"duplicate abundance table for replicate: {replicate}")
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"abundance table is missing or empty for {replicate}: {path}")
        fields, rows = read_tsv(path)
        missing = sorted(PFPKM_REQUIRED_COLUMNS - set(fields))
        if missing:
            raise ValueError(
                f"abundance table for {replicate} is missing column(s): {', '.join(missing)}"
            )
        by_candidate: Dict[str, Dict[str, str]] = {}
        for line_number, row in enumerate(rows, start=2):
            gorf_id = str(row.get("gorf_id", "")).strip()
            overlap_type = str(row.get("overlap_type", "")).strip()
            if not gorf_id or not overlap_type:
                raise ValueError(
                    f"{path}:{line_number} has empty gorf_id or overlap_type"
                )
            key = f"{gorf_id}|{overlap_type}"
            if key in by_candidate:
                raise ValueError(f"{path}:{line_number} duplicates candidate {key}")
            value = safe_float(row.get("intorf_pFPKM"))
            status = str(row.get("abundance_status", "")).strip()
            if math.isfinite(value) and value < 0:
                raise ValueError(f"{path}:{line_number} has negative intorf_pFPKM")
            if status == "quantified" and not math.isfinite(value):
                raise ValueError(
                    f"{path}:{line_number} is quantified but intorf_pFPKM is not finite"
                )
            by_candidate[key] = dict(row)
        values[replicate] = by_candidate
        paths[replicate] = path
    return values, paths


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


def replicate_columns(fields: Sequence[str]) -> List[str]:
    prefix = "replicate::"
    suffix = "::present"
    return [field[len(prefix):-len(suffix)] for field in fields if field.startswith(prefix) and field.endswith(suffix)]


def read_sample_metadata(path: Path) -> Dict[str, Dict[str, str]]:
    """Read optional sample-level sequencing/QC metrics used as heatmap annotations."""
    fields, rows = read_tsv(path)
    canonical = {"replicate_id", "psite_alignments", "frame0_prop", "A0"}
    batch_summary = {
        "replicate_id", "psite_alignments", "qc_frame0_prop", "qc_A0",
    }
    if not canonical.issubset(fields) and not batch_summary.issubset(fields):
        missing = sorted(canonical - set(fields))
        raise ValueError(
            "sample metadata must contain canonical columns or batch-summary QC "
            f"columns; missing canonical column(s): {', '.join(missing)}"
        )
    result: Dict[str, Dict[str, str]] = {}
    for row in rows:
        row = dict(row)
        if "frame0_prop" not in row:
            row["frame0_prop"] = row.get("qc_frame0_prop", "")
        if "A0" not in row:
            row["A0"] = row.get("qc_A0", "")
        replicate = str(row.get("replicate_id", "")).strip()
        if not replicate:
            raise ValueError("sample metadata contains an empty replicate_id")
        if replicate in result:
            raise ValueError(f"sample metadata contains duplicate replicate_id: {replicate}")
        result[replicate] = dict(row)
    return result


def compact_count(value: object) -> str:
    count = safe_float(value)
    if not math.isfinite(count):
        return "NA"
    if count >= 1_000_000:
        return f"{count / 1_000_000:.1f}M"
    if count >= 1_000:
        return f"{count / 1_000:.1f}k"
    return f"{count:.0f}"


def replicate_quality_labels(
    replicates: Sequence[str], sample_metadata: Mapping[str, Mapping[str, str]]
) -> List[str]:
    """Return compact sample labels with P-site depth, frame-0 fraction, and A0."""
    labels = []
    for replicate in replicates:
        metadata = sample_metadata.get(replicate)
        if not metadata:
            labels.append(replicate)
            continue
        frame0 = safe_float(metadata.get("frame0_prop"))
        a0 = safe_float(metadata.get("A0"))
        frame0_text = f"{100 * frame0:.1f}%" if math.isfinite(frame0) else "NA"
        a0_text = f"{a0:.2f}" if math.isfinite(a0) else "NA"
        lengths = [
            token.split(":", 1)[0].strip()
            for token in str(metadata.get("selected_length_offsets", "")).split(",")
            if token.strip()
        ]
        depth_text = f"{compact_count(metadata.get('psite_alignments'))} P-sites"
        if lengths:
            depth_text += f" · L {'/'.join(lengths)}"
        labels.append(
            f"{replicate}\n{depth_text}\n"
            f"F0 {frame0_text} · A0 {a0_text}"
        )
    return labels


def save_figure(fig: plt.Figure, prefix: Path, formats: Sequence[str], dpi: int) -> List[Path]:
    outputs = []
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in formats:
        destination = Path(f"{prefix}.{extension}")
        temporary = Path(f"{destination}.tmp.{extension}")
        fig.savefig(temporary, dpi=dpi, bbox_inches="tight")
        temporary.replace(destination)
        outputs.append(destination)
    plt.close(fig)
    return outputs


def call_counts(
    long_rows: Sequence[Mapping[str, str]], replicates: Sequence[str], fdr_threshold: float
) -> Dict[str, Dict[str, int]]:
    result = {
        replicate: {"total": 0, "significant": 0, "primary": 0}
        for replicate in replicates
    }
    for row in long_rows:
        replicate = str(row.get("replicate_id", ""))
        if replicate not in result:
            continue
        result[replicate]["total"] += 1
        q_value = safe_float(row.get("q_BH"))
        result[replicate]["significant"] += int(math.isfinite(q_value) and q_value <= fdr_threshold)
        result[replicate]["primary"] += int(safe_int(row.get("primary_credible_call")) == 1)
    return result


def plot_call_counts(
    counts: Mapping[str, Mapping[str, int]], replicates: Sequence[str]
) -> plt.Figure:
    height = max(4.2, 0.38 * len(replicates) + 2.2)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, height), sharey=True)
    y = np.arange(len(replicates))
    panels = (
        ("total", "Candidates observed", COLORS["total"]),
        ("significant", "BH significant", COLORS["significant"]),
        ("primary", "Primary credible", COLORS["primary"]),
    )
    for index, (metric, title, color) in enumerate(panels):
        values = [counts[replicate][metric] for replicate in replicates]
        axes[index].barh(y, values, color=color, edgecolor="none")
        axes[index].set_title(title)
        axes[index].set_xlabel("Candidate count")
        axes[index].grid(axis="x", color="#DDDDDD", linewidth=0.6)
        axes[index].set_axisbelow(True)
        for yi, value in zip(y, values):
            axes[index].text(value, yi, f" {value:,}", va="center", fontsize=8)
    axes[0].set_yticks(y, replicates)
    axes[0].invert_yaxis()
    fig.suptitle("Replicate-level intORF call counts", fontweight="normal")
    fig.tight_layout()
    return fig


def support_combinations(
    consensus_rows: Sequence[Mapping[str, str]], replicates: Sequence[str]
) -> Counter[Tuple[str, ...]]:
    counts: Counter[Tuple[str, ...]] = Counter()
    for row in consensus_rows:
        support = tuple(
            replicate for replicate in replicates
            if safe_int(row.get(f"replicate::{replicate}::primary_credible_call")) == 1
        )
        if support:
            counts[support] += 1
    return counts


def plot_support_combinations(
    combination_counts: Mapping[Tuple[str, ...], int],
    replicates: Sequence[str],
    max_combinations: int,
) -> plt.Figure:
    ordered = sorted(
        combination_counts.items(), key=lambda item: (-item[1], -len(item[0]), item[0])
    )[:max_combinations]
    if not ordered:
        fig, ax = plt.subplots(figsize=(8.5, 3.5))
        ax.text(0.5, 0.5, "No primary credible calls in available replicates", ha="center", va="center")
        ax.axis("off")
        fig.suptitle("Primary credible support combinations")
        return fig

    width = max(8.5, 0.55 * len(ordered) + 4.0)
    height = max(5.0, 0.35 * len(replicates) + 3.5)
    fig = plt.figure(figsize=(width, height))
    grid = fig.add_gridspec(2, 1, height_ratios=(2.1, max(1.2, 0.28 * len(replicates))), hspace=0.05)
    bar_ax = fig.add_subplot(grid[0])
    matrix_ax = fig.add_subplot(grid[1], sharex=bar_ax)
    x = np.arange(len(ordered))
    values = [count for _support, count in ordered]
    bar_ax.bar(x, values, color=COLORS["primary"], edgecolor="none")
    bar_ax.set_ylabel("Candidate count")
    bar_ax.set_title("Primary credible support combinations")
    bar_ax.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    bar_ax.set_axisbelow(True)
    bar_ax.tick_params(axis="x", labelbottom=False)
    for xi, value in zip(x, values):
        bar_ax.text(xi, value, f"{value:,}", ha="center", va="bottom", fontsize=8)

    for yi, replicate in enumerate(replicates):
        for xi, (support, _count) in enumerate(ordered):
            active = replicate in support
            matrix_ax.scatter(
                xi, yi, s=45 if active else 24,
                color=COLORS["primary"] if active else COLORS["absent"],
                edgecolor="none", zorder=2,
            )
        matrix_ax.axhline(yi, color="#EEEEEE", linewidth=0.5, zorder=0)
    for xi, (support, _count) in enumerate(ordered):
        active_y = [replicates.index(item) for item in support]
        if len(active_y) > 1:
            matrix_ax.plot([xi, xi], [min(active_y), max(active_y)], color="#555555", linewidth=1.0, zorder=1)
    matrix_ax.set_yticks(np.arange(len(replicates)), replicates)
    matrix_ax.set_xticks(x, [str(index + 1) for index in x])
    matrix_ax.set_xlabel("Support combination (ranked by candidate count)")
    matrix_ax.set_ylim(len(replicates) - 0.5, -0.5)
    for spine in ("top", "right"):
        matrix_ax.spines[spine].set_visible(False)
        bar_ax.spines[spine].set_visible(False)
    return fig


class LambdaPair(NamedTuple):
    left: str
    right: str
    x: np.ndarray
    y: np.ndarray
    left_primary: np.ndarray
    right_primary: np.ndarray


def paired_lambda_values(
    long_rows: Sequence[Mapping[str, str]], replicates: Sequence[str]
) -> List[LambdaPair]:
    values: Dict[str, Dict[str, Tuple[float, bool]]] = {replicate: {} for replicate in replicates}
    for row in long_rows:
        replicate = str(row.get("replicate_id", ""))
        value = safe_float(row.get("lambda_hat"))
        key = str(row.get("candidate_key", ""))
        if replicate in values and key and math.isfinite(value):
            values[replicate][key] = (
                value,
                safe_int(row.get("primary_credible_call")) == 1,
            )
    pairs = []
    for left, right in combinations(replicates, 2):
        shared = sorted(set(values[left]) & set(values[right]))
        if not shared:
            continue
        x = np.asarray([values[left][key][0] for key in shared], dtype=float)
        y = np.asarray([values[right][key][0] for key in shared], dtype=float)
        left_primary = np.asarray([values[left][key][1] for key in shared], dtype=bool)
        right_primary = np.asarray([values[right][key][1] for key in shared], dtype=bool)
        pairs.append(LambdaPair(left, right, x, y, left_primary, right_primary))
    return sorted(pairs, key=lambda item: (-len(item.x), item.left, item.right))


def lambda_values_by_replicate(
    long_rows: Sequence[Mapping[str, str]],
    replicates: Sequence[str],
    *,
    primary_union_only: bool,
) -> Dict[str, np.ndarray]:
    """Return finite lambda values for diagonal distribution panels."""
    primary_union = {
        str(row.get("candidate_key", ""))
        for row in long_rows
        if safe_int(row.get("primary_credible_call")) == 1
    }
    values: Dict[str, List[float]] = {replicate: [] for replicate in replicates}
    for row in long_rows:
        replicate = str(row.get("replicate_id", ""))
        key = str(row.get("candidate_key", ""))
        value = safe_float(row.get("lambda_hat"))
        if (
            replicate in values
            and math.isfinite(value)
            and (not primary_union_only or key in primary_union)
        ):
            values[replicate].append(value)
    return {
        replicate: np.asarray(replicate_values, dtype=float)
        for replicate, replicate_values in values.items()
    }


def _plot_lambda_concordance_matrix(
    pairs: Sequence[LambdaPair],
    long_rows: Sequence[Mapping[str, str]],
    replicates: Sequence[str],
    *,
    primary_union_only: bool,
) -> plt.Figure:
    title = (
        "Primary-union λ concordance (lower-triangle scatter matrix)"
        if primary_union_only
        else "All-eligible λ concordance (lower-triangle scatter matrix)"
    )
    if not replicates:
        fig, ax = plt.subplots(figsize=(8.5, 3.5))
        ax.text(0.5, 0.5, "No replicate is available", ha="center", va="center")
        ax.axis("off")
        fig.suptitle(title)
        return fig

    size = len(replicates)
    panel_size = max(2.1, min(3.35, 19.0 / size))
    fig, axes = plt.subplots(
        size, size,
        figsize=(panel_size * size, 0.96 * panel_size * size),
        squeeze=False,
    )
    pair_by_names = {(pair.left, pair.right): pair for pair in pairs}
    diagonal_values = lambda_values_by_replicate(
        long_rows, replicates, primary_union_only=primary_union_only
    )
    for row_index, row_name in enumerate(replicates):
        for column_index, column_name in enumerate(replicates):
            ax = axes[row_index, column_index]
            if column_index > row_index:
                ax.axis("off")
                continue
            if row_index == column_index:
                values = diagonal_values[row_name]
                if len(values):
                    bins = min(35, max(10, int(math.sqrt(len(values)))))
                    ax.hist(
                        values, bins=bins, range=(0, 1),
                        color=COLORS["primary"] if primary_union_only else COLORS["total"],
                        alpha=0.82,
                    )
                    ax.text(
                        0.04, 0.94,
                        f"{row_name}\nn={len(values):,}\nmedian={np.median(values):.3f}",
                        transform=ax.transAxes, va="top", fontsize=8,
                    )
                else:
                    qualifier = "primary-union " if primary_union_only else ""
                    ax.text(
                        0.5, 0.5, f"{row_name}\nno {qualifier}finite λ",
                        ha="center", va="center",
                    )
                ax.set_xlim(-0.02, 1.02)
                ax.set_ylabel("Candidates")
                ax.grid(axis="y", color="#E5E5E5", linewidth=0.5)
            else:
                pair = pair_by_names.get((column_name, row_name))
                if pair is None:
                    ax.text(0.5, 0.5, "No shared candidates with finite λ", ha="center", va="center")
                    ax.set_xticks([])
                    ax.set_yticks([])
                    continue
                if primary_union_only:
                    union = pair.left_primary | pair.right_primary
                    both = pair.left_primary & pair.right_primary
                    x = pair.x[union]
                    y = pair.y[union]
                    both_union = both[union]
                    one_sided_union = ~both_union
                    ax.scatter(
                        x[one_sided_union], y[one_sided_union],
                        s=10, alpha=0.32, color=COLORS["significant"],
                        edgecolor="none", rasterized=True,
                    )
                    ax.scatter(
                        x[both_union], y[both_union],
                        s=11, alpha=0.48, color=COLORS["primary"],
                        edgecolor="none", rasterized=True,
                    )
                else:
                    x, y = pair.x, pair.y
                    both_union = np.asarray([], dtype=bool)
                    ax.scatter(
                        x, y, s=7, alpha=0.22, color=COLORS["total"],
                        edgecolor="none", rasterized=True,
                    )
                if not len(x):
                    ax.text(
                        0.5, 0.5,
                        "No shared candidates primary credible in either replicate",
                        ha="center", va="center",
                    )
                    ax.set_xticks([])
                    ax.set_yticks([])
                    continue
                ax.plot([0, 1], [0, 1], color="#666666", linewidth=0.8, linestyle="--")
                pearson = (
                    float(np.corrcoef(x, y)[0, 1])
                    if len(x) > 1 and np.std(x) > 0 and np.std(y) > 0
                    else float("nan")
                )
                statistics = [
                    f"union n={len(x):,}" if primary_union_only else f"n={len(x):,}"
                ]
                if primary_union_only:
                    statistics.append(f"both={int(np.sum(both_union)):,}")
                    spearman = (
                        float(spearmanr(x, y).statistic)
                        if len(x) > 1 and np.std(x) > 0 and np.std(y) > 0
                        else float("nan")
                    )
                    median_abs_difference = float(np.median(np.abs(x - y)))
                    if math.isfinite(spearman):
                        statistics.append(f"rho={spearman:.3f}")
                    statistics.append(f"median |delta lambda|={median_abs_difference:.3f}")
                if math.isfinite(pearson):
                    statistics.insert(2 if primary_union_only else 1, f"r={pearson:.3f}")
                ax.text(
                    0.04, 0.96, "\n".join(statistics),
                    transform=ax.transAxes, va="top", fontsize=8,
                )
                ax.set_xlim(-0.02, 1.02)
                ax.set_ylim(-0.02, 1.02)
                ax.set_aspect("equal", adjustable="box")
                ax.grid(color="#E5E5E5", linewidth=0.5)
            if row_index == size - 1:
                ax.set_xlabel(f"λ\n{column_name}")
            elif row_index != column_index:
                ax.set_xticklabels([])
            if column_index == 0 and row_index > 0:
                ax.set_ylabel(f"{row_name}\nλ")
            elif column_index > 0 and row_index != column_index:
                ax.set_yticklabels([])
    if primary_union_only:
        legend_items = [
            plt.Line2D(
                [0], [0], marker="o", color="none",
                markerfacecolor=COLORS["primary"], markeredgecolor="none",
                markersize=6, label="Primary in both",
            ),
            plt.Line2D(
                [0], [0], marker="o", color="none",
                markerfacecolor=COLORS["significant"], markeredgecolor="none",
                markersize=6, label="Primary in one replicate",
            ),
        ]
        fig.legend(
            handles=legend_items, loc="upper center", bbox_to_anchor=(0.5, 0.975),
            ncol=2, frameon=False,
        )
        fig.suptitle(title, y=0.998)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
    else:
        fig.suptitle(title)
        fig.tight_layout(rect=(0, 0, 1, 0.975))
    return fig


def plot_all_eligible_lambda_concordance(
    pairs: Sequence[LambdaPair],
    long_rows: Sequence[Mapping[str, str]],
    replicates: Sequence[str],
) -> plt.Figure:
    return _plot_lambda_concordance_matrix(
        pairs, long_rows, replicates, primary_union_only=False
    )


def plot_lambda_concordance(
    pairs: Sequence[LambdaPair],
    long_rows: Sequence[Mapping[str, str]],
    replicates: Sequence[str],
) -> plt.Figure:
    return _plot_lambda_concordance_matrix(
        pairs, long_rows, replicates, primary_union_only=True
    )


def reproducibility_matrix(
    consensus_rows: Sequence[Mapping[str, str]],
    replicates: Sequence[str],
    max_candidates: int,
    *,
    group_by_pattern: bool = False,
) -> Tuple[np.ndarray, List[str], List[Mapping[str, str]], set[str]]:
    unavailable = {
        replicate for replicate in replicates
        if not any(safe_int(row.get(f"replicate::{replicate}::present")) == 1 for row in consensus_rows)
    }
    eligible = [row for row in consensus_rows if safe_int(row.get("n_primary_credible")) > 0]
    def row_state(row: Mapping[str, str], replicate: str) -> int:
        if replicate in unavailable:
            return -1
        if safe_int(row.get(f"replicate::{replicate}::present")) != 1:
            return 0
        return 2 if safe_int(row.get(f"replicate::{replicate}::primary_credible_call")) == 1 else 1

    if group_by_pattern:
        eligible.sort(
            key=lambda row: (
                -safe_int(row.get("n_primary_credible")),
                tuple(-row_state(row, replicate) for replicate in replicates),
                -safe_int(row.get("max_core_reads")),
                str(row.get("candidate_key", "")),
            )
        )
    else:
        eligible.sort(
            key=lambda row: (
                -safe_int(row.get("n_primary_credible")),
                -safe_int(row.get("max_core_reads")),
                str(row.get("candidate_key", "")),
            )
        )
    shown = eligible[:max_candidates]
    matrix = np.zeros((len(shown), len(replicates)), dtype=int)
    for yi, row in enumerate(shown):
        for xi, replicate in enumerate(replicates):
            matrix[yi, xi] = row_state(row, replicate)
    labels = candidate_display_labels(shown)
    return matrix, labels, shown, unavailable


def candidate_display_labels(rows: Sequence[Mapping[str, str]]) -> List[str]:
    """Return readable, unambiguous labels while preserving candidate identity in tables."""
    bases = []
    for row in rows:
        gene_name = str(row.get("gene_name", "")).strip()
        gene_id = str(row.get("gene_id", "")).strip()
        candidate_key = str(row.get("candidate_key", "")).strip()
        bases.append(gene_name or gene_id or candidate_key or "unknown candidate")

    base_counts = Counter(bases)
    labels = []
    for base, row in zip(bases, rows):
        if base_counts[base] == 1:
            labels.append(base)
            continue

        details = []
        overlap_type = str(row.get("overlap_type", "")).strip()
        if overlap_type:
            details.append(overlap_type)

        gorf_id = str(row.get("gorf_id", "")).strip()
        hash_token = next(
            (part for part in gorf_id.split("|") if part.startswith("h") and len(part) > 1),
            "",
        )
        if hash_token:
            details.append(f"h...{hash_token[-6:]}")

        candidate_key = str(row.get("candidate_key", "")).strip()
        labels.append(f"{base} ({'; '.join(details)})" if details else candidate_key or base)
    return labels


class PFPKMPair(NamedTuple):
    left: str
    right: str
    x: np.ndarray
    y: np.ndarray
    pearson: float
    spearman: float


def abundance_value(
    abundance: Mapping[str, Mapping[str, Mapping[str, str]]],
    replicate: str,
    candidate_key: str,
) -> float:
    row = abundance.get(replicate, {}).get(candidate_key)
    if not row or str(row.get("abundance_status", "")) != "quantified":
        return float("nan")
    value = safe_float(row.get("intorf_pFPKM"))
    return value if math.isfinite(value) and value >= 0 else float("nan")


def pfpkm_candidate_rows(
    consensus_rows: Sequence[Mapping[str, str]],
    replicates: Sequence[str],
    abundance: Mapping[str, Mapping[str, Mapping[str, str]]],
) -> Tuple[List[Mapping[str, str]], np.ndarray]:
    """Return candidates primary in any replicate and their raw pFPKM matrix."""
    eligible = [
        row for row in consensus_rows
        if safe_int(row.get("n_primary_credible")) > 0
    ]
    eligible.sort(
        key=lambda row: (
            -safe_int(row.get("n_primary_credible")),
            str(row.get("candidate_key", "")),
        )
    )
    matrix = np.full((len(eligible), len(replicates)), np.nan, dtype=float)
    for row_index, row in enumerate(eligible):
        key = str(row.get("candidate_key", ""))
        for column_index, replicate in enumerate(replicates):
            matrix[row_index, column_index] = abundance_value(abundance, replicate, key)
    keep = np.sum(np.isfinite(matrix), axis=1) > 0 if len(eligible) else np.asarray([], dtype=bool)
    kept_rows = [row for row, selected in zip(eligible, keep) if selected]
    kept_matrix = matrix[keep] if matrix.size else np.empty((0, len(replicates)))
    if len(kept_rows):
        medians = np.asarray([
            float(np.nanmedian(row)) if np.any(np.isfinite(row)) else float("-inf")
            for row in kept_matrix
        ])
        support = np.asarray([safe_int(row.get("n_primary_credible")) for row in kept_rows])
        order = np.lexsort((
            np.asarray([str(row.get("candidate_key", "")) for row in kept_rows]),
            -medians,
            -support,
        ))
        kept_rows = [kept_rows[index] for index in order]
        kept_matrix = kept_matrix[order]
    return kept_rows, kept_matrix


def paired_pfpkm_values(
    rows: Sequence[Mapping[str, str]],
    matrix: np.ndarray,
    replicates: Sequence[str],
) -> List[PFPKMPair]:
    del rows  # row selection is already encoded by the aligned matrix
    pairs: List[PFPKMPair] = []
    for left_index, right_index in combinations(range(len(replicates)), 2):
        complete = np.isfinite(matrix[:, left_index]) & np.isfinite(matrix[:, right_index])
        x = np.log2(1.0 + matrix[complete, left_index])
        y = np.log2(1.0 + matrix[complete, right_index])
        pearson = (
            float(np.corrcoef(x, y)[0, 1])
            if len(x) > 1 and np.std(x) > 0 and np.std(y) > 0
            else float("nan")
        )
        spearman = (
            float(spearmanr(x, y).statistic)
            if len(x) > 1 and np.std(x) > 0 and np.std(y) > 0
            else float("nan")
        )
        pairs.append(
            PFPKMPair(
                replicates[left_index], replicates[right_index], x, y, pearson, spearman
            )
        )
    return pairs


def plot_pfpkm_correlation(
    pairs: Sequence[PFPKMPair],
    matrix: np.ndarray,
    replicates: Sequence[str],
) -> plt.Figure:
    if not replicates:
        fig, ax = plt.subplots(figsize=(8.5, 3.5))
        ax.text(0.5, 0.5, "No abundance-enabled replicate is available", ha="center", va="center")
        ax.axis("off")
        fig.suptitle("Pairwise intORF pFPKM concordance")
        return fig

    size = len(replicates)
    # Keep small batches spacious while capping the physical canvas for larger
    # batches.  The number and arrangement of cells always follow the input.
    panel_size = max(2.1, min(3.35, 19.0 / size))
    fig, axes = plt.subplots(
        size, size,
        figsize=(panel_size * size, 0.96 * panel_size * size),
        squeeze=False,
    )
    pair_by_names = {(pair.left, pair.right): pair for pair in pairs}
    for row_index, row_name in enumerate(replicates):
        for column_index, column_name in enumerate(replicates):
            ax = axes[row_index, column_index]
            if column_index > row_index:
                ax.axis("off")
                continue
            if row_index == column_index:
                values = matrix[:, row_index]
                values = np.log2(1.0 + values[np.isfinite(values)])
                if len(values):
                    bins = min(35, max(10, int(math.sqrt(len(values)))))
                    ax.hist(values, bins=bins, color=COLORS["primary"], alpha=0.82)
                    ax.text(
                        0.04, 0.94,
                        f"{row_name}\nn={len(values):,}\nmedian={np.median(values):.2f}",
                        transform=ax.transAxes, va="top", fontsize=8,
                    )
                else:
                    ax.text(0.5, 0.5, f"{row_name}\nno quantified pFPKM", ha="center", va="center")
                ax.set_ylabel("Candidates")
                ax.grid(axis="y", color="#E5E5E5", linewidth=0.5)
            else:
                pair = pair_by_names.get((column_name, row_name))
                if pair is None or not len(pair.x):
                    ax.text(0.5, 0.5, "No shared quantified candidates", ha="center", va="center")
                    ax.set_xticks([])
                    ax.set_yticks([])
                    continue
                combined = np.concatenate((pair.x, pair.y))
                low = min(0.0, float(np.min(combined)))
                high = max(1.0, float(np.max(combined)))
                padding = 0.04 * max(high - low, 1.0)
                ax.scatter(
                    pair.x, pair.y, s=8, alpha=0.28, color=COLORS["primary"],
                    edgecolor="none", rasterized=True,
                )
                ax.plot(
                    [low, high], [low, high], color="#666666",
                    linewidth=0.8, linestyle="--",
                )
                statistics_text = [f"n={len(pair.x):,}"]
                if math.isfinite(pair.pearson):
                    statistics_text.append(f"r={pair.pearson:.3f}")
                if math.isfinite(pair.spearman):
                    statistics_text.append(f"rho={pair.spearman:.3f}")
                ax.text(
                    0.04, 0.96, "\n".join(statistics_text),
                    transform=ax.transAxes, va="top", fontsize=8,
                )
                ax.set_xlim(low - padding, high + padding)
                ax.set_ylim(low - padding, high + padding)
                ax.set_aspect("equal", adjustable="box")
                ax.grid(color="#E5E5E5", linewidth=0.5)
            if row_index == size - 1:
                ax.set_xlabel(f"log2(1 + pFPKM)\n{column_name}")
            elif row_index != column_index:
                ax.set_xticklabels([])
            if column_index == 0 and row_index > 0:
                ax.set_ylabel(f"{row_name}\nlog2(1 + pFPKM)")
            elif column_index > 0 and row_index != column_index:
                ax.set_yticklabels([])
    fig.suptitle("intORF pFPKM concordance (lower-triangle scatter matrix)")
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    return fig


def pfpkm_correlation_matrices(
    pairs: Sequence[PFPKMPair], replicates: Sequence[str]
) -> Tuple[np.ndarray, np.ndarray]:
    correlation = np.full((len(replicates), len(replicates)), np.nan, dtype=float)
    counts = np.zeros((len(replicates), len(replicates)), dtype=int)
    for index in range(len(replicates)):
        correlation[index, index] = 1.0
    index_by_name = {name: index for index, name in enumerate(replicates)}
    for pair in pairs:
        left = index_by_name[pair.left]
        right = index_by_name[pair.right]
        correlation[left, right] = correlation[right, left] = pair.spearman
        counts[left, right] = counts[right, left] = len(pair.x)
    return correlation, counts


def plot_pfpkm_correlation_matrix(
    correlation: np.ndarray, counts: np.ndarray, replicates: Sequence[str]
) -> plt.Figure:
    size = max(5.2, 0.9 * len(replicates) + 2.2)
    fig, ax = plt.subplots(figsize=(size, size))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(COLORS["unavailable"])
    image = ax.imshow(np.ma.masked_invalid(correlation), vmin=-1, vmax=1, cmap=cmap)
    ax.set_xticks(np.arange(len(replicates)), replicates, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(replicates)), replicates)
    for row in range(len(replicates)):
        for column in range(len(replicates)):
            value = correlation[row, column]
            if not math.isfinite(value):
                label = "NA"
            elif row == column:
                label = "1.000"
            else:
                label = f"{value:.3f}\nn={counts[row, column]:,}"
            ax.text(column, row, label, ha="center", va="center", fontsize=8)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Spearman rho of log2(1 + pFPKM)")
    ax.set_title("intORF pFPKM replicate correlation")
    fig.tight_layout()
    return fig


def lambda_candidate_rows(
    consensus_rows: Sequence[Mapping[str, str]],
    long_rows: Sequence[Mapping[str, str]],
    replicates: Sequence[str],
) -> Tuple[List[Mapping[str, str]], np.ndarray]:
    """Align lambda_hat for candidates primary credible in at least one replicate."""
    eligible = [
        row for row in consensus_rows
        if safe_int(row.get("n_primary_credible")) > 0
    ]
    eligible.sort(
        key=lambda row: (
            -safe_int(row.get("n_primary_credible")),
            -safe_float(row.get("max_core_reads")),
            str(row.get("candidate_key", "")),
        )
    )
    replicate_index = {replicate: index for index, replicate in enumerate(replicates)}
    row_index = {
        str(row.get("candidate_key", "")): index for index, row in enumerate(eligible)
    }
    matrix = np.full((len(eligible), len(replicates)), np.nan, dtype=float)
    for row in long_rows:
        key = str(row.get("candidate_key", ""))
        replicate = str(row.get("replicate_id", ""))
        if key not in row_index or replicate not in replicate_index:
            continue
        value = safe_float(row.get("lambda_hat"))
        if math.isfinite(value) and 0 <= value <= 1:
            matrix[row_index[key], replicate_index[replicate]] = value
    keep = np.sum(np.isfinite(matrix), axis=1) > 0 if len(eligible) else np.asarray([], dtype=bool)
    rows = [row for row, selected in zip(eligible, keep) if selected]
    return rows, matrix[keep] if matrix.size else np.empty((0, len(replicates)))


def plot_lambda_heatmap(
    rows: Sequence[Mapping[str, str]],
    matrix: np.ndarray,
    replicates: Sequence[str],
    total_candidates: int,
    sample_metadata: Mapping[str, Mapping[str, str]] | None = None,
) -> plt.Figure:
    height = max(5.0, min(15.0, 0.065 * max(len(rows), 1) + 2.8))
    width = max(7.5, 1.3 * len(replicates) + 5.5)
    fig, ax = plt.subplots(figsize=(width, height))
    if matrix.size == 0:
        ax.text(0.5, 0.5, "No primary-union candidate has finite lambda_hat", ha="center", va="center")
        ax.axis("off")
        fig.suptitle("intORF mixture-fraction heatmap")
        return fig
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(COLORS["unavailable"])
    image = ax.imshow(
        np.ma.masked_invalid(matrix), aspect="auto", interpolation="nearest",
        cmap=cmap, vmin=0, vmax=1,
    )
    quality_labels = replicate_quality_labels(replicates, sample_metadata or {})
    labels = candidate_display_labels(rows)
    ax.set_xticks(np.arange(len(replicates)), quality_labels, fontsize=8)
    if len(labels) <= 60:
        ax.set_yticks(np.arange(len(labels)), labels, fontsize=7)
    else:
        tick_count = min(15, len(labels))
        ticks = np.unique(np.linspace(0, len(labels) - 1, tick_count, dtype=int))
        ax.set_yticks(ticks, [labels[index] for index in ticks], fontsize=7)
    ax.set_xlabel("Replicate")
    ax.set_ylabel("Gene / candidate (primary support, then core reads)")
    shown = f"top {len(rows):,} of {total_candidates:,}" if len(rows) < total_candidates else f"all {len(rows):,}"
    ax.set_title(f"Model-estimated intORF mixture fraction ({shown} primary-union candidates)")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("lambda_hat (0 = annotated CDS; 1 = intORF)")
    fig.tight_layout()
    return fig


def write_lambda_matrix(
    path: Path,
    rows: Sequence[Mapping[str, str]],
    matrix: np.ndarray,
    replicates: Sequence[str],
) -> None:
    columns = [
        "candidate_key", "gorf_id", "overlap_type", "gene_id", "gene_name",
        "n_primary_credible", "max_core_reads",
        *[f"replicate::{replicate}::lambda_hat" for replicate in replicates],
    ]
    records: List[Dict[str, object]] = []
    for row_index, row in enumerate(rows):
        record: Dict[str, object] = {
            column: row.get(column, "") for column in columns[:7]
        }
        for column_index, replicate in enumerate(replicates):
            value = matrix[row_index, column_index]
            record[f"replicate::{replicate}::lambda_hat"] = (
                f"{value:.12g}" if math.isfinite(value) else ""
            )
        records.append(record)
    write_table(path, columns, records)


def plot_pfpkm_heatmap(
    rows: Sequence[Mapping[str, str]],
    matrix: np.ndarray,
    replicates: Sequence[str],
    total_candidates: int,
    sample_metadata: Mapping[str, Mapping[str, str]] | None = None,
) -> plt.Figure:
    height = max(5.0, min(15.0, 0.085 * max(len(rows), 1) + 2.8))
    width = max(8.5, 1.4 * len(replicates) + 5.5)
    fig, ax = plt.subplots(figsize=(width, height))
    if matrix.size == 0:
        ax.text(0.5, 0.5, "No primary-union candidate has quantified pFPKM", ha="center", va="center")
        ax.axis("off")
        fig.suptitle("intORF pFPKM heatmap")
        return fig

    log_matrix = np.log2(1.0 + matrix)
    abundance_cmap = plt.get_cmap("viridis").copy()
    abundance_cmap.set_bad(COLORS["unavailable"])
    finite_log = log_matrix[np.isfinite(log_matrix)]
    upper = float(np.quantile(finite_log, 0.98)) if len(finite_log) else 1.0
    upper = max(upper, 1.0)
    image_abundance = ax.imshow(
        np.ma.masked_invalid(log_matrix), aspect="auto", interpolation="nearest",
        cmap=abundance_cmap, vmin=0, vmax=upper,
    )
    quality_labels = replicate_quality_labels(replicates, sample_metadata or {})
    labels = candidate_display_labels(rows)
    ax.set_xticks(np.arange(len(replicates)), quality_labels, fontsize=8)
    ax.set_xlabel("Replicate")
    if len(labels) <= 60:
        ax.set_yticks(np.arange(len(labels)), labels, fontsize=7)
    else:
        tick_count = min(15, len(labels))
        ticks = np.unique(np.linspace(0, len(labels) - 1, tick_count, dtype=int))
        ax.set_yticks(ticks, [labels[index] for index in ticks], fontsize=7)
    ax.set_ylabel("Gene / candidate (primary support, then median pFPKM)")
    ax.set_title("Absolute abundance")
    colorbar_absolute = fig.colorbar(image_abundance, ax=ax, fraction=0.046, pad=0.04)
    colorbar_absolute.set_label("log2(1 + intORF pFPKM)")
    shown = f"top {len(rows):,} of {total_candidates:,}" if len(rows) < total_candidates else f"all {len(rows):,}"
    fig.suptitle(f"intORF pFPKM heatmap ({shown} candidates primary in any replicate)")
    fig.tight_layout()
    return fig


def write_pfpkm_matrix(
    path: Path,
    rows: Sequence[Mapping[str, str]],
    matrix: np.ndarray,
    replicates: Sequence[str],
    abundance: Mapping[str, Mapping[str, Mapping[str, str]]],
) -> None:
    columns = [
        "candidate_key", "gorf_id", "overlap_type", "gene_id", "gene_name",
        "n_primary_credible",
    ]
    for replicate in replicates:
        columns.extend([
            f"replicate::{replicate}::intorf_pFPKM",
            f"replicate::{replicate}::abundance_status",
        ])
    records: List[Dict[str, object]] = []
    for row_index, row in enumerate(rows):
        key = str(row.get("candidate_key", ""))
        record: Dict[str, object] = {column: row.get(column, "") for column in columns[:6]}
        for column_index, replicate in enumerate(replicates):
            value = matrix[row_index, column_index]
            abundance_row = abundance.get(replicate, {}).get(key, {})
            record[f"replicate::{replicate}::intorf_pFPKM"] = (
                f"{value:.12g}" if math.isfinite(value) else ""
            )
            record[f"replicate::{replicate}::abundance_status"] = abundance_row.get(
                "abundance_status", "unavailable"
            )
        records.append(record)
    write_table(path, columns, records)


def write_pfpkm_correlations(path: Path, pairs: Sequence[PFPKMPair]) -> None:
    columns = [
        "replicate_left", "replicate_right", "shared_candidates",
        "pearson_log2p1_pFPKM", "spearman_log2p1_pFPKM",
    ]
    records = [
        {
            "replicate_left": pair.left,
            "replicate_right": pair.right,
            "shared_candidates": len(pair.x),
            "pearson_log2p1_pFPKM": f"{pair.pearson:.12g}" if math.isfinite(pair.pearson) else "",
            "spearman_log2p1_pFPKM": f"{pair.spearman:.12g}" if math.isfinite(pair.spearman) else "",
        }
        for pair in pairs
    ]
    write_table(path, columns, records)


def plot_reproducibility_heatmap(
    matrix: np.ndarray,
    labels: Sequence[str],
    replicates: Sequence[str],
    total_eligible: int,
    sample_metadata: Mapping[str, Mapping[str, str]] | None = None,
) -> plt.Figure:
    height = max(4.5, min(14.0, 0.075 * max(len(labels), 1) + 2.6))
    width = max(
        7.5,
        (1.7 * len(replicates) + 5.8) if sample_metadata else (0.55 * len(replicates) + 4.8),
    )
    fig, ax = plt.subplots(figsize=(width, height))
    if matrix.size == 0:
        ax.text(0.5, 0.5, "No candidate was primary credible in any replicate", ha="center", va="center")
        ax.axis("off")
        fig.suptitle("Primary credible reproducibility")
        return fig
    cmap = ListedColormap([COLORS["unavailable"], COLORS["absent"], COLORS["neutral"], COLORS["primary"]])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5, 2.5], cmap.N)
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
    quality_labels = replicate_quality_labels(replicates, sample_metadata or {})
    ax.set_xticks(np.arange(len(replicates)), quality_labels, rotation=0, ha="center", fontsize=8)
    if len(labels) <= 40:
        ax.set_yticks(np.arange(len(labels)), labels, fontsize=7)
    else:
        tick_count = min(12, len(labels))
        ticks = np.unique(np.linspace(0, len(labels) - 1, tick_count, dtype=int))
        ax.set_yticks(ticks, [labels[index] for index in ticks], fontsize=7)
    ax.set_xlabel("Replicate")
    ax.set_ylabel("Gene / candidate (ranked by replicate support, then core reads)")
    shown_note = f"top {len(labels):,} of {total_eligible:,}" if len(labels) < total_eligible else f"all {len(labels):,}"
    ax.set_title(f"Primary credible reproducibility ({shown_note} candidates with any primary call)")
    legend_items = [
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=COLORS["primary"], markeredgecolor="none", markersize=8, label="Primary credible"),
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=COLORS["neutral"], markeredgecolor="none", markersize=8, label="Observed, not primary"),
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=COLORS["absent"], markeredgecolor="none", markersize=8, label="Candidate absent"),
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=COLORS["unavailable"], markeredgecolor="none", markersize=8, label="Replicate unavailable"),
    ]
    ax.legend(handles=legend_items, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    fig.tight_layout()
    return fig


def plot_global_reproducibility_heatmap(
    matrix: np.ndarray,
    rows: Sequence[Mapping[str, str]],
    replicates: Sequence[str],
    height: float,
    sample_metadata: Mapping[str, Mapping[str, str]] | None = None,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(max(8.2, 0.55 * len(replicates) + 5.2), height))
    if matrix.size == 0:
        ax.text(0.5, 0.5, "No candidate was primary credible in any replicate", ha="center", va="center")
        ax.axis("off")
        fig.suptitle("Global primary credible reproducibility")
        return fig

    cmap = ListedColormap([COLORS["unavailable"], COLORS["absent"], COLORS["neutral"], COLORS["primary"]])
    norm = BoundaryNorm([-1.5, -0.5, 0.5, 1.5, 2.5], cmap.N)
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
    quality_labels = replicate_quality_labels(replicates, sample_metadata or {})
    ax.set_xticks(np.arange(len(replicates)), quality_labels, rotation=0, ha="center", fontsize=8)
    ax.set_yticks([])
    ax.set_xlabel("Replicate")
    ax.set_ylabel("All candidates with any primary call (grouped by support and call pattern)")
    ax.set_title(f"Global primary credible reproducibility (all {len(rows):,} candidates)")

    support = [safe_int(row.get("n_primary_credible")) for row in rows]
    for index in range(1, len(support)):
        if support[index] != support[index - 1]:
            ax.axhline(index - 0.5, color="white", linewidth=0.7)

    legend_items = [
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=COLORS["primary"], markeredgecolor="none", markersize=8, label="Primary credible"),
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=COLORS["neutral"], markeredgecolor="none", markersize=8, label="Observed, not primary"),
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=COLORS["absent"], markeredgecolor="none", markersize=8, label="Candidate absent"),
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=COLORS["unavailable"], markeredgecolor="none", markersize=8, label="Replicate unavailable"),
    ]
    ax.legend(handles=legend_items, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    fig.tight_layout()
    return fig


def write_summary(path: Path, rows: Sequence[Tuple[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{path}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(rows)
    temporary.replace(path)


def write_table(
    path: Path, columns: Sequence[str], rows: Sequence[Mapping[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{path}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def filter_gene_level_pure_intorfs(
    long_rows: Sequence[Mapping[str, str]],
    consensus_rows: Sequence[Mapping[str, str]],
) -> Tuple[List[Mapping[str, str]], List[Mapping[str, str]], int]:
    """Filter integrated rows using the post-DM same-gene CDS-context decision."""
    filtered_consensus = [
        row for row in consensus_rows
        if safe_int(row.get("gene_level_pure_intorf_eligible")) == 1
    ]
    allowed_keys = {
        str(row.get("candidate_key", "")) for row in filtered_consensus
    }
    filtered_long = [
        row for row in long_rows
        if str(row.get("candidate_key", "")) in allowed_keys
    ]
    return filtered_long, filtered_consensus, len(consensus_rows) - len(filtered_consensus)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog=PROGRAM, description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--replicate-long", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument(
        "--gene-level-pure-intorf-only",
        action="store_true",
        help=(
            "Require gene-level context columns on the consensus table and exclude "
            "annotated-CDS-derived candidates from final replicate figures and counts"
        ),
    )
    parser.add_argument(
        "--sample-metadata",
        help=(
            "Optional TSV with replicate_id, psite_alignments, frame0_prop, and A0; "
            "selected_length_offsets is also shown when present. Batch-summary "
            "qc_frame0_prop/qc_A0 columns are accepted directly"
        ),
    )
    parser.add_argument(
        "--abundance",
        action="append",
        default=[],
        metavar="REPLICATE_ID=INTORF_ABUNDANCE.tsv",
        help="Optional model-allocated intORF pFPKM table; repeat per replicate",
    )
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--max-combinations", type=int, default=20)
    parser.add_argument(
        "--max-lambda-pairs", type=int, default=12,
        help=(
            "Deprecated compatibility option; lambda scatter matrices always "
            "show every replicate pair"
        ),
    )
    parser.add_argument("--max-heatmap-candidates", type=int, default=200)
    parser.add_argument("--max-pfpkm-heatmap-candidates", type=int, default=100)
    parser.add_argument(
        "--global-heatmap-height",
        type=float,
        default=7.2,
        help="Height in inches of the compact all-primary-candidate heatmap (default: 7.2)",
    )
    parser.add_argument("--formats", default="png,pdf")
    parser.add_argument("--dpi", type=int, default=220)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not 0 < args.fdr_threshold < 1:
        raise SystemExit("ERROR: --fdr-threshold must be in (0, 1)")
    if min(
        args.max_combinations, args.max_lambda_pairs, args.max_heatmap_candidates,
        args.max_pfpkm_heatmap_candidates, args.dpi,
    ) < 1:
        raise SystemExit("ERROR: plot limits and --dpi must be positive")
    if not math.isfinite(args.global_heatmap_height) or args.global_heatmap_height <= 0:
        raise SystemExit("ERROR: --global-heatmap-height must be a finite positive number")
    formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    if not formats or any(item not in {"png", "pdf"} for item in formats):
        raise SystemExit("ERROR: --formats must contain png and/or pdf")
    try:
        long_fields, long_rows = read_tsv(Path(args.replicate_long))
        consensus_fields, consensus_rows = read_tsv(Path(args.consensus))
        sample_metadata = (
            read_sample_metadata(Path(args.sample_metadata)) if args.sample_metadata else {}
        )
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    required_long = {"candidate_key", "replicate_id", "q_BH", "lambda_hat", "core_reads", "primary_credible_call"}
    missing_long = sorted(required_long - set(long_fields))
    if missing_long:
        raise SystemExit(f"ERROR: replicate-long missing column(s): {', '.join(missing_long)}")
    required_consensus = {"candidate_key", "n_primary_credible", "max_core_reads"}
    if args.gene_level_pure_intorf_only:
        required_consensus.update(
            {"gene_level_orf_class", "gene_level_pure_intorf_eligible"}
        )
    missing_consensus = sorted(required_consensus - set(consensus_fields))
    if missing_consensus:
        raise SystemExit(f"ERROR: consensus missing column(s): {', '.join(missing_consensus)}")
    replicates = replicate_columns(consensus_fields)
    if not replicates:
        raise SystemExit("ERROR: consensus contains no replicate::<ID>::present columns")
    raw_consensus_rows = consensus_rows
    raw_consensus_keys = {
        str(row.get("candidate_key", "")) for row in raw_consensus_rows
    }
    gene_context_excluded = 0
    if args.gene_level_pure_intorf_only:
        long_rows, consensus_rows, gene_context_excluded = filter_gene_level_pure_intorfs(
            long_rows, raw_consensus_rows
        )
    try:
        abundance, abundance_paths = read_abundance_tables(args.abundance, replicates)
        for replicate, rows in abundance.items():
            unexpected = sorted(set(rows) - raw_consensus_keys)
            if unexpected:
                raise ValueError(
                    f"abundance table for {replicate} contains candidate(s) absent from "
                    f"the integration consensus; first: {unexpected[0]}"
                )
    except (OSError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc

    prefix = Path(args.out_prefix)
    outputs: List[Path] = []
    counts = call_counts(long_rows, replicates, args.fdr_threshold)
    outputs.extend(save_figure(plot_call_counts(counts, replicates), Path(f"{prefix}.call_counts"), formats, args.dpi))
    combinations_count = support_combinations(consensus_rows, replicates)
    outputs.extend(save_figure(
        plot_support_combinations(combinations_count, replicates, args.max_combinations),
        Path(f"{prefix}.support_combinations"), formats, args.dpi,
    ))
    pairs = paired_lambda_values(long_rows, replicates)
    outputs.extend(save_figure(
        plot_lambda_concordance(pairs, long_rows, replicates),
        Path(f"{prefix}.lambda_concordance"), formats, args.dpi,
    ))
    outputs.extend(save_figure(
        plot_all_eligible_lambda_concordance(pairs, long_rows, replicates),
        Path(f"{prefix}.lambda_concordance_all_eligible"), formats, args.dpi,
    ))
    any_primary_total = sum(safe_int(row.get("n_primary_credible")) > 0 for row in consensus_rows)
    matrix, labels, shown, unavailable = reproducibility_matrix(
        consensus_rows, replicates, args.max_heatmap_candidates
    )
    outputs.extend(save_figure(
        plot_reproducibility_heatmap(
            matrix, labels, replicates, any_primary_total, sample_metadata
        ),
        Path(f"{prefix}.primary_reproducibility"), formats, args.dpi,
    ))
    global_matrix, _global_labels, global_rows, _global_unavailable = reproducibility_matrix(
        consensus_rows, replicates, any_primary_total, group_by_pattern=True
    )
    outputs.extend(save_figure(
        plot_global_reproducibility_heatmap(
            global_matrix, global_rows, replicates, args.global_heatmap_height,
            sample_metadata,
        ),
        Path(f"{prefix}.primary_reproducibility_all"), formats, args.dpi,
    ))

    lambda_rows, lambda_matrix = lambda_candidate_rows(
        consensus_rows, long_rows, replicates
    )
    lambda_heatmap_count = min(len(lambda_rows), args.max_heatmap_candidates)
    outputs.extend(save_figure(
        plot_lambda_heatmap(
            lambda_rows[:lambda_heatmap_count], lambda_matrix[:lambda_heatmap_count],
            replicates, len(lambda_rows), sample_metadata,
        ),
        Path(f"{prefix}.lambda_heatmap"), formats, args.dpi,
    ))
    lambda_matrix_path = Path(f"{prefix}.lambda_matrix.tsv")
    write_lambda_matrix(lambda_matrix_path, lambda_rows, lambda_matrix, replicates)
    outputs.append(lambda_matrix_path)

    pfpkm_rows: List[Mapping[str, str]] = []
    pfpkm_matrix = np.empty((0, len(replicates)), dtype=float)
    pfpkm_pairs: List[PFPKMPair] = []
    if abundance:
        pfpkm_rows, pfpkm_matrix = pfpkm_candidate_rows(
            consensus_rows, replicates, abundance
        )
        pfpkm_pairs = paired_pfpkm_values(pfpkm_rows, pfpkm_matrix, replicates)
        outputs.extend(save_figure(
            plot_pfpkm_correlation(pfpkm_pairs, pfpkm_matrix, replicates),
            Path(f"{prefix}.pfpkm_correlation"), formats, args.dpi,
        ))
        correlation_matrix, correlation_counts = pfpkm_correlation_matrices(
            pfpkm_pairs, replicates
        )
        outputs.extend(save_figure(
            plot_pfpkm_correlation_matrix(
                correlation_matrix, correlation_counts, replicates
            ),
            Path(f"{prefix}.pfpkm_correlation_matrix"), formats, args.dpi,
        ))
        heatmap_count = min(len(pfpkm_rows), args.max_pfpkm_heatmap_candidates)
        outputs.extend(save_figure(
            plot_pfpkm_heatmap(
                pfpkm_rows[:heatmap_count], pfpkm_matrix[:heatmap_count], replicates,
                len(pfpkm_rows), sample_metadata,
            ),
            Path(f"{prefix}.pfpkm_heatmap"), formats, args.dpi,
        ))
        pfpkm_matrix_path = Path(f"{prefix}.pfpkm_matrix.tsv")
        pfpkm_correlations_path = Path(f"{prefix}.pfpkm_correlations.tsv")
        write_pfpkm_matrix(
            pfpkm_matrix_path, pfpkm_rows, pfpkm_matrix, replicates, abundance
        )
        write_pfpkm_correlations(pfpkm_correlations_path, pfpkm_pairs)
        outputs.extend([pfpkm_matrix_path, pfpkm_correlations_path])

    summary_path = Path(f"{prefix}.plot_summary.tsv")
    manifest_path = Path(f"{prefix}.plot_manifest.json")
    summary_rows: List[Tuple[str, object]] = [
        ("program", PROGRAM), ("version", VERSION),
        ("replicates", ",".join(replicates)),
        ("unavailable_replicates", ",".join(sorted(unavailable))),
        ("candidate_union_total", len(consensus_rows)),
        ("candidates_with_any_primary_call", any_primary_total),
        ("support_combinations_total", len(combinations_count)),
        ("lambda_pairs_available", len(pairs)),
        ("lambda_pairs_plotted", len(pairs)),
        ("lambda_candidates_primary_union", len(lambda_rows)),
        ("lambda_heatmap_candidates_plotted", lambda_heatmap_count),
        ("lambda_candidate_set", "primary_credible_in_any_available_replicate"),
        ("heatmap_candidates_plotted", len(shown)),
        ("global_heatmap_candidates_plotted", len(global_rows)),
        ("sample_metadata_provided", int(bool(args.sample_metadata))),
        ("gene_level_pure_intorf_only", int(args.gene_level_pure_intorf_only)),
        ("candidate_union_total_raw", len(raw_consensus_rows)),
        ("gene_context_excluded_candidates", gene_context_excluded),
        ("pfpkm_abundance_replicates", ",".join(abundance)),
        ("pfpkm_candidates_primary_union", len(pfpkm_rows)),
        ("pfpkm_pairs_available", len(pfpkm_pairs)),
        (
            "pfpkm_heatmap_candidates_plotted",
            min(len(pfpkm_rows), args.max_pfpkm_heatmap_candidates),
        ),
        ("pfpkm_transform", "log2(1 + intorf_pFPKM)"),
        ("pfpkm_candidate_set", "primary_credible_in_any_available_replicate"),
    ]
    for replicate in replicates:
        for metric in ("total", "significant", "primary"):
            summary_rows.append((f"replicate::{replicate}::{metric}", counts[replicate][metric]))
        metadata = sample_metadata.get(replicate, {})
        for metric in ("psite_alignments", "frame0_prop", "A0"):
            summary_rows.append((f"replicate::{replicate}::{metric}", metadata.get(metric, "")))
    write_summary(summary_path, summary_rows)
    manifest = {
        "program": PROGRAM,
        "version": VERSION,
        "inputs": {
            "replicate_long": args.replicate_long,
            "consensus": args.consensus,
            "sample_metadata": args.sample_metadata or "",
            "abundance": {key: str(value) for key, value in abundance_paths.items()},
        },
        "parameters": {
            "fdr_threshold": args.fdr_threshold,
            "max_combinations": args.max_combinations,
            "max_lambda_pairs": args.max_lambda_pairs,
            "max_heatmap_candidates": args.max_heatmap_candidates,
            "max_pfpkm_heatmap_candidates": args.max_pfpkm_heatmap_candidates,
            "global_heatmap_height": args.global_heatmap_height,
            "formats": formats,
            "dpi": args.dpi,
            "gene_level_pure_intorf_only": args.gene_level_pure_intorf_only,
        },
        "outputs": [str(path) for path in [*outputs, summary_path]],
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{manifest_path}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    temporary.replace(manifest_path)
    print(f"[ok] wrote {len(outputs)} plot/data files")
    print(f"[ok] summary={summary_path}")
    print(f"[ok] manifest={manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
