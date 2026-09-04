#!/usr/bin/env python3
"""Create standalone, auditable visual summaries from an intORF DM result table.

The module never changes caller results.  It rotates F1/F2 phase compositions
into one annotated-CDS coordinate system and visualizes the fitted H->T model,
post-BH gates, and their margins.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


SQRT3 = math.sqrt(3.0)
EXCLUDED_PRE_FDR = {
    "invalid_overlap_type",
    "too_short",
    "insufficient_data",
    "template_separation_too_small",
}
CLASS_COLORS = {
    "credible_extra_ORF_like_signal": "#0072B2",
    "ultra_short_exploratory": "#56B4E9",
    "atypical_pattern": "#D55E00",
    "localized_core_signal": "#CC79A7",
    "low_target_frame_residual_breadth": "#E69F00",
    "host_only_supported": "#8A8A8A",
}
CLASS_LABELS = {
    "credible_extra_ORF_like_signal": "credible",
    "ultra_short_exploratory": "ultra-short exploratory",
    "atypical_pattern": "atypical geometry",
    "localized_core_signal": "localized core signal",
    "low_target_frame_residual_breadth": "low residual breadth",
    "host_only_supported": "host only",
}
PROGRAM = "miao-orf-visualize"
VERSION = "1.0.0"


def parse_fraction(value: str) -> float:
    text = str(value).strip()
    if "/" in text:
        numerator, denominator = text.split("/", 1)
        result = float(numerator) / float(denominator)
    else:
        result = float(text)
    if not 0.0 <= result <= 1.0:
        raise argparse.ArgumentTypeError("fraction must be in [0,1]")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog=PROGRAM, description=__doc__)
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    parser.add_argument("--dm-results", required=True)
    parser.add_argument(
        "--gene-level-pure-intorf-only",
        action="store_true",
        help=(
            "Require gene-level context columns and exclude candidates that reuse "
            "an annotated CDS N terminus before making final figures and count tables"
        ),
    )
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--lambda-min", type=float, default=0.05)
    parser.add_argument("--lambda-abs-diff-max", type=float, default=0.10)
    parser.add_argument("--lambda-rel-diff-max", type=float, default=0.30)
    parser.add_argument("--distance-to-segment-max", type=float, default=0.10)
    parser.add_argument("--min-active-core-codons", type=int, default=5)
    parser.add_argument("--min-active-core-frac", type=float, default=0.15)
    parser.add_argument("--min-target-residual-frac", type=parse_fraction, default=1.0 / 3.0)
    parser.add_argument("--highlight-gene", action="append", default=[])
    parser.add_argument("--max-background-points", type=int, default=30000)
    parser.add_argument("--heatmap-rows", type=int, default=36)
    parser.add_argument(
        "--include-diagnostics",
        action="store_true",
        help=(
            "Also write developer-oriented branch projection, gain/drop, and "
            "candidate gate-margin diagnostics. These are omitted from the "
            "default end-user result summary."
        ),
    )
    parser.add_argument("--formats", default="png,pdf")
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--seed", type=int, default=20260803)
    return parser.parse_args()


def ffloat(value: object, default: float = float("nan")) -> float:
    try:
        result = float(str(value))
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def fint(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def rotate_to_cds(values: Sequence[float], overlap_type: str) -> np.ndarray:
    """Rotate candidate-frame P0/P1/P2 into CDS P0/+1/+2 coordinates."""
    p0, p1, p2 = map(float, values)
    if overlap_type == "F1":
        return np.asarray([p2, p0, p1], dtype=float)
    if overlap_type == "F2":
        return np.asarray([p1, p2, p0], dtype=float)
    raise ValueError(f"Unsupported overlap type: {overlap_type}")


def ternary_xy(comp: Sequence[float]) -> Tuple[float, float]:
    p0, p1, p2 = map(float, comp)
    return (0.5 * p0 + p2, (SQRT3 / 2.0) * p0)


def segment_distance(point: np.ndarray, start: np.ndarray, end: np.ndarray) -> float:
    direction = end - start
    denom = float(direction @ direction)
    if denom <= 0:
        return float(np.linalg.norm(point - start))
    lam = float(np.clip(((point - start) @ direction) / denom, 0.0, 1.0))
    return float(np.linalg.norm(point - (start + lam * direction)))


@dataclass
class Point:
    torf_id: str
    gorf_id: str
    gene_name: str
    gene_id: str
    overlap_type: str
    classification: str
    classification_pre_fdr: str
    fdr_family: str
    peptide_len: int
    pi_trans: np.ndarray
    pi_obs_candidate: np.ndarray
    pi_obs_cds: np.ndarray
    x: float
    y: float
    q_bh: float
    p_final: float
    core_reads: int
    n_active: int
    active_frac: float
    residual_frac: float
    lambda_hat: float
    lambda_gain: float
    lambda_drop: float
    lambda_abs_diff: float
    lambda_rel_diff: float
    lambda_proj: float
    distance: float
    alt_gain: float
    host_drop: float
    geometry_consistent: int
    gene_level_orf_class: str
    gene_level_pure_intorf_eligible: int


def require_columns(columns: Sequence[str]) -> None:
    required = {
        "torf_id", "gorf_id", "gene_name", "gene_id", "overlap_type",
        "classification", "classification_pre_fdr", "fdr_family", "peptide_len",
        "pi_trans_0", "pi_trans_1", "pi_trans_2", "pi_obs_0", "pi_obs_1", "pi_obs_2",
        "q_BH", "p_final", "core_reads", "n_active_core_codons", "active_core_codon_frac",
        "target_vs_unused_z05_frac_active", "lambda_hat", "lambda_gain_raw",
        "lambda_drop_raw", "lambda_abs_diff", "lambda_rel_diff", "lambda_proj_raw",
        "distance_to_mixture_segment", "alt_frame_gain", "host_frame_drop",
        "mixture_geometry_consistent",
    }
    missing = sorted(required.difference(columns))
    if missing:
        raise SystemExit(f"ERROR: DM table lacks required columns: {missing}")


def load_points(
    path: str, gene_level_pure_intorf_only: bool = False
) -> Tuple[List[Point], Dict[str, int]]:
    stats: Dict[str, int] = {
        "candidate_total_raw": 0,
        "candidate_after_gene_context": 0,
        "gene_context_excluded_rows": 0,
        "eligible_geometry_rows": 0,
    }
    points: List[Point] = []
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"ERROR: input has no header: {path}")
        require_columns(reader.fieldnames)
        if gene_level_pure_intorf_only:
            required_context = {
                "gene_level_orf_class", "gene_level_pure_intorf_eligible",
            }
            missing_context = sorted(required_context.difference(reader.fieldnames))
            if missing_context:
                raise SystemExit(
                    "ERROR: gene-level pure-intORF filtering requested but input lacks "
                    f"context column(s): {missing_context}"
                )
        for row in reader:
            stats["candidate_total_raw"] += 1
            pure_eligible = fint(row.get("gene_level_pure_intorf_eligible"), 1)
            if gene_level_pure_intorf_only and pure_eligible != 1:
                stats["gene_context_excluded_rows"] += 1
                continue
            stats["candidate_after_gene_context"] += 1
            overlap = str(row.get("overlap_type", ""))
            pi_obs = np.asarray([
                ffloat(row.get("pi_obs_0")),
                ffloat(row.get("pi_obs_1")),
                ffloat(row.get("pi_obs_2")),
            ])
            pi_trans = np.asarray([
                ffloat(row.get("pi_trans_0")),
                ffloat(row.get("pi_trans_1")),
                ffloat(row.get("pi_trans_2")),
            ])
            if overlap not in {"F1", "F2"} or not np.all(np.isfinite(pi_obs)) or not np.all(np.isfinite(pi_trans)):
                continue
            cds = rotate_to_cds(pi_obs, overlap)
            x, y = ternary_xy(cds)
            points.append(Point(
                torf_id=str(row.get("torf_id", "")),
                gorf_id=str(row.get("gorf_id", "")),
                gene_name=str(row.get("gene_name", "")),
                gene_id=str(row.get("gene_id", "")),
                overlap_type=overlap,
                classification=str(row.get("classification", "")),
                classification_pre_fdr=str(row.get("classification_pre_fdr", "")),
                fdr_family=str(row.get("fdr_family", "")),
                peptide_len=fint(row.get("peptide_len")),
                pi_trans=pi_trans,
                pi_obs_candidate=pi_obs,
                pi_obs_cds=cds,
                x=x,
                y=y,
                q_bh=ffloat(row.get("q_BH"), 1.0),
                p_final=ffloat(row.get("p_final"), 1.0),
                core_reads=fint(row.get("core_reads")),
                n_active=fint(row.get("n_active_core_codons")),
                active_frac=ffloat(row.get("active_core_codon_frac"), 0.0),
                residual_frac=ffloat(row.get("target_vs_unused_z05_frac_active"), 0.0),
                lambda_hat=ffloat(row.get("lambda_hat")),
                lambda_gain=ffloat(row.get("lambda_gain_raw")),
                lambda_drop=ffloat(row.get("lambda_drop_raw")),
                lambda_abs_diff=ffloat(row.get("lambda_abs_diff")),
                lambda_rel_diff=ffloat(row.get("lambda_rel_diff")),
                lambda_proj=ffloat(row.get("lambda_proj_raw")),
                distance=ffloat(row.get("distance_to_mixture_segment")),
                alt_gain=ffloat(row.get("alt_frame_gain")),
                host_drop=ffloat(row.get("host_frame_drop")),
                geometry_consistent=fint(row.get("mixture_geometry_consistent")),
                gene_level_orf_class=str(row.get("gene_level_orf_class", "")),
                gene_level_pure_intorf_eligible=pure_eligible,
            ))
    stats["eligible_geometry_rows"] = len(points)
    return points, stats


def model_templates(point: Point) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    host = point.pi_trans.copy()
    if point.overlap_type == "F1":
        target = np.asarray([point.pi_trans[2], point.pi_trans[0], point.pi_trans[1]])
    else:
        target = np.asarray([point.pi_trans[1], point.pi_trans[2], point.pi_trans[0]])
    return host, target, point.pi_obs_cds


def gate_passes(point: Point, args: argparse.Namespace) -> Dict[str, bool]:
    return {
        "BH": point.q_bh < args.fdr_threshold,
        "direction": (
            math.isfinite(point.lambda_hat) and point.lambda_hat > args.lambda_min
            and math.isfinite(point.alt_gain) and point.alt_gain > 0
            and math.isfinite(point.host_drop) and point.host_drop < 0
        ),
        "geometry": (
            point.geometry_consistent == 1
            and math.isfinite(point.distance)
            and point.distance < args.distance_to_segment_max
        ),
        "active breadth": (
            point.n_active >= args.min_active_core_codons
            and point.active_frac >= args.min_active_core_frac
        ),
        "residual breadth": point.residual_frac > args.min_target_residual_frac,
    }


def gate_margins(point: Point, args: argparse.Namespace) -> List[float]:
    def finite(value: float, fallback: float = -2.0) -> float:
        return value if math.isfinite(value) else fallback

    return [
        finite((args.fdr_threshold - point.q_bh) / args.fdr_threshold),
        finite((point.lambda_hat - args.lambda_min) / args.lambda_min),
        finite((args.distance_to_segment_max - point.distance) / args.distance_to_segment_max),
        finite((args.lambda_abs_diff_max - point.lambda_abs_diff) / args.lambda_abs_diff_max),
        finite((args.lambda_rel_diff_max - point.lambda_rel_diff) / args.lambda_rel_diff_max),
        (point.n_active - args.min_active_core_codons) / max(args.min_active_core_codons, 1),
        finite((point.active_frac - args.min_active_core_frac) / max(args.min_active_core_frac, 1e-12)),
        finite((point.residual_frac - args.min_target_residual_frac) / max(args.min_target_residual_frac, 1e-12)),
    ]


def save_figure(fig: plt.Figure, prefix: str, formats: Sequence[str], dpi: int) -> None:
    for extension in formats:
        fig.savefig(f"{prefix}.{extension}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def draw_ternary_frame(ax: plt.Axes) -> None:
    vertices = np.asarray([[0.5, SQRT3 / 2.0], [0.0, 0.0], [1.0, 0.0], [0.5, SQRT3 / 2.0]])
    ax.plot(vertices[:, 0], vertices[:, 1], color="#333333", lw=1.2)
    for value in (0.2, 0.4, 0.6, 0.8):
        # Constant CDS P0, +1, and +2 grid lines.
        segments = [
            ([value, 1 - value, 0], [value, 0, 1 - value]),
            ([1 - value, value, 0], [0, value, 1 - value]),
            ([1 - value, 0, value], [0, 1 - value, value]),
        ]
        for a, b in segments:
            xa, ya = ternary_xy(a)
            xb, yb = ternary_xy(b)
            ax.plot([xa, xb], [ya, yb], color="#D7D7D7", lw=0.55, zorder=0)
    ax.text(0.5, SQRT3 / 2.0 + 0.035, "CDS P0", ha="center", va="bottom", fontsize=10)
    ax.text(-0.025, -0.025, "CDS +1 (F1 target)", ha="left", va="top", fontsize=10)
    ax.text(1.025, -0.025, "CDS +2 (F2 target)", ha="right", va="top", fontsize=10)
    ax.set_xlim(-0.06, 1.06)
    ax.set_ylim(-0.06, SQRT3 / 2.0 + 0.08)
    ax.set_aspect("equal")
    ax.axis("off")


def add_distance_corridors(ax: plt.Axes, template: np.ndarray, distance_max: float) -> None:
    comps: List[Tuple[float, float, float]] = []
    step = 0.02
    n = int(round(1.0 / step))
    for i in range(n + 1):
        for j in range(n + 1 - i):
            p0 = i * step
            p1 = j * step
            comps.append((p0, p1, 1.0 - p0 - p1))
    arr = np.asarray(comps)
    f1 = np.asarray([template[2], template[0], template[1]])
    f2 = np.asarray([template[1], template[2], template[0]])
    values = np.asarray([
        min(segment_distance(comp, template, f1), segment_distance(comp, template, f2))
        for comp in arr
    ])
    xy = np.asarray([ternary_xy(comp) for comp in arr])
    tri = mtri.Triangulation(xy[:, 0], xy[:, 1])
    ax.tricontourf(
        tri, values, levels=[0.0, distance_max, float(values.max()) + 1.0],
        colors=["#CFE8F3", "#FFFFFF"], alpha=0.34, zorder=0,
    )


def add_single_distance_corridor(
    ax: plt.Axes,
    template: np.ndarray,
    endpoint: np.ndarray,
    distance_max: float,
) -> None:
    comps: List[Tuple[float, float, float]] = []
    step = 0.02
    n = int(round(1.0 / step))
    for i in range(n + 1):
        for j in range(n + 1 - i):
            p0 = i * step
            p1 = j * step
            comps.append((p0, p1, 1.0 - p0 - p1))
    arr = np.asarray(comps)
    values = np.asarray([segment_distance(comp, template, endpoint) for comp in arr])
    xy = np.asarray([ternary_xy(comp) for comp in arr])
    tri = mtri.Triangulation(xy[:, 0], xy[:, 1])
    ax.tricontourf(
        tri, values, levels=[0.0, distance_max, float(values.max()) + 1.0],
        colors=["#CFE8F3", "#FFFFFF"], alpha=0.34, zorder=0,
    )


def plot_ternary(points: Sequence[Point], out_prefix: str, args: argparse.Namespace, formats: Sequence[str]) -> None:
    fig, ax = plt.subplots(figsize=(9.2, 8.0))
    draw_ternary_frame(ax)
    template = points[0].pi_trans
    add_distance_corridors(ax, template, args.distance_to_segment_max)

    rng = random.Random(args.seed)
    background = list(points)
    if len(background) > args.max_background_points:
        background = rng.sample(background, args.max_background_points)
    ax.scatter([p.x for p in background], [p.y for p in background], s=4, c="#777777", alpha=0.09,
               linewidths=0, rasterized=True, label=f"Eligible candidates (n={len(points):,})")

    significant = [p for p in points if p.q_bh < args.fdr_threshold]
    order = ["atypical_pattern", "localized_core_signal", "low_target_frame_residual_breadth",
             "ultra_short_exploratory", "credible_extra_ORF_like_signal"]
    for classification in order:
        subset = [p for p in significant if p.classification == classification]
        if not subset:
            continue
        marker = "o" if classification != "credible_extra_ORF_like_signal" else "D"
        ax.scatter([p.x for p in subset], [p.y for p in subset], s=18 if marker == "o" else 24,
                   c=CLASS_COLORS[classification], marker=marker, alpha=0.72,
                   edgecolors="none" if marker == "o" else "#1F1F1F", linewidths=0.25,
                   rasterized=True, label=f"{classification} (n={len(subset):,})")

    host = template
    endpoints = {
        "F1": np.asarray([template[2], template[0], template[1]]),
        "F2": np.asarray([template[1], template[2], template[0]]),
    }
    hx, hy = ternary_xy(host)
    ax.scatter([hx], [hy], s=85, c="#000000", marker="X", zorder=7)
    ax.text(hx, hy + 0.045, "H: annotated CDS template", ha="center", fontsize=9)
    for label, endpoint in endpoints.items():
        ex, ey = ternary_xy(endpoint)
        ax.plot([hx, ex], [hy, ey], lw=2.2, color="#222222", zorder=4)
        for lam in (0.25, 0.5, 0.75, 1.0):
            comp = (1.0 - lam) * host + lam * endpoint
            tx, ty = ternary_xy(comp)
            ax.scatter([tx], [ty], s=20, facecolor="#FFFFFF", edgecolor="#222222", linewidth=0.7, zorder=6)
            if lam in (0.5, 1.0):
                ax.text(tx, ty - 0.022, f"{lam:g}", fontsize=7.5, ha="center", va="top")
        ax.text(ex, ey + (0.02 if label == "F1" else -0.025), f"T({label})", fontsize=9, ha="center")

    highlights = [p for p in points if p.gene_name in set(args.highlight_gene)]
    if highlights:
        names = ",".join(sorted({p.gene_name for p in highlights}))
        ax.scatter([p.x for p in highlights], [p.y for p in highlights], s=92, marker="*",
                   c="#F0E442", edgecolor="#111111", linewidth=0.7, zorder=9,
                   label=f"Highlighted: {names} (n={len(highlights)})")

    ax.set_title("intORF DM geometry in annotated-CDS phase coordinates", fontsize=13)
    ax.text(0.5, -0.075,
            f"Shaded corridor: distance < {args.distance_to_segment_max:g}; λ ticks are reference mixtures, not λ-hat gates",
            ha="center", va="top", fontsize=9)
    ax.legend(loc="upper right", bbox_to_anchor=(1.30, 1.0), frameon=False, fontsize=8)
    save_figure(fig, f"{out_prefix}.model_geometry_ternary", formats, args.dpi)


def plot_ternary_by_frame(
    points: Sequence[Point],
    out_prefix: str,
    args: argparse.Namespace,
    formats: Sequence[str],
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 6.5))
    template = points[0].pi_trans
    endpoints = {
        "F1": np.asarray([template[2], template[0], template[1]]),
        "F2": np.asarray([template[1], template[2], template[0]]),
    }
    rng = random.Random(args.seed + 11)
    legend_handles = []
    legend_labels = []
    for ax, overlap in zip(axes, ("F1", "F2")):
        draw_ternary_frame(ax)
        endpoint = endpoints[overlap]
        add_single_distance_corridor(ax, template, endpoint, args.distance_to_segment_max)
        subset = [p for p in points if p.overlap_type == overlap]
        max_points = max(1000, args.max_background_points // 2)
        background = subset if len(subset) <= max_points else rng.sample(subset, max_points)
        background_handle = ax.scatter(
            [p.x for p in background], [p.y for p in background], s=4, c="#777777",
            alpha=0.09, linewidths=0, rasterized=True,
        )
        for classification in (
            "atypical_pattern", "localized_core_signal",
            "low_target_frame_residual_breadth", "ultra_short_exploratory",
            "credible_extra_ORF_like_signal",
        ):
            selected = [
                p for p in subset
                if p.q_bh < args.fdr_threshold and p.classification == classification
            ]
            if not selected:
                continue
            marker = "D" if classification == "credible_extra_ORF_like_signal" else "o"
            handle = ax.scatter(
                [p.x for p in selected], [p.y for p in selected],
                s=20 if marker == "D" else 15, c=CLASS_COLORS[classification],
                marker=marker, alpha=0.72, linewidths=0, rasterized=True,
            )
            if overlap == "F1":
                legend_handles.append(handle)
                legend_labels.append(classification)
        hx, hy = ternary_xy(template)
        ex, ey = ternary_xy(endpoint)
        ax.plot([hx, ex], [hy, ey], color="#222222", lw=2.0)
        ax.scatter([hx], [hy], s=60, c="#000000", marker="X", zorder=6)
        for lam in (0.25, 0.5, 0.75, 1.0):
            comp = (1.0 - lam) * template + lam * endpoint
            tx, ty = ternary_xy(comp)
            ax.scatter([tx], [ty], s=17, facecolor="#FFFFFF", edgecolor="#222222", linewidth=0.6, zorder=6)
        ax.set_title(f"{overlap}: eligible n={len(subset):,}")
    legend_handles.insert(0, background_handle)
    legend_labels.insert(0, "eligible background")
    fig.legend(legend_handles, legend_labels, loc="center right", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
    fig.suptitle("Frame-specific DM geometry in annotated-CDS coordinates")
    fig.subplots_adjust(right=0.84, wspace=0.12)
    save_figure(fig, f"{out_prefix}.model_geometry_ternary_by_frame", formats, args.dpi)


def signed_branch_distance(point: Point) -> float:
    host, target, observed = model_templates(point)
    hx, hy = ternary_xy(host)
    tx, ty = ternary_xy(target)
    ox, oy = ternary_xy(observed)
    cross = (tx - hx) * (oy - hy) - (ty - hy) * (ox - hx)
    sign = 1.0 if cross >= 0 else -1.0
    return sign * point.distance


def plot_branch_geometry(points: Sequence[Point], out_prefix: str, args: argparse.Namespace, formats: Sequence[str]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.0), sharex=True, sharey=True)
    rng = random.Random(args.seed + 1)
    for ax, overlap in zip(axes, ("F1", "F2")):
        subset = [p for p in points if p.overlap_type == overlap and math.isfinite(p.lambda_hat) and math.isfinite(p.distance)]
        background = subset if len(subset) <= args.max_background_points // 2 else rng.sample(subset, args.max_background_points // 2)
        ax.scatter([p.lambda_hat for p in background], [signed_branch_distance(p) for p in background],
                   s=5, c="#777777", alpha=0.10, linewidths=0, rasterized=True)
        for classification in ("atypical_pattern", "localized_core_signal", "low_target_frame_residual_breadth",
                               "ultra_short_exploratory", "credible_extra_ORF_like_signal"):
            selected = [p for p in subset if p.q_bh < args.fdr_threshold and p.classification == classification]
            if selected:
                ax.scatter([p.lambda_hat for p in selected], [signed_branch_distance(p) for p in selected],
                           s=14, c=CLASS_COLORS[classification], alpha=0.70, linewidths=0, rasterized=True)
        ax.axhspan(-args.distance_to_segment_max, args.distance_to_segment_max, color="#56B4E9", alpha=0.10)
        ax.axhline(args.distance_to_segment_max, color="#555555", lw=0.8, ls="--")
        ax.axhline(-args.distance_to_segment_max, color="#555555", lw=0.8, ls="--")
        ax.axvline(args.lambda_min, color="#222222", lw=1.0, ls=":")
        ax.set_title(f"{overlap} candidates")
        ax.set_xlabel("Codon-level DM fitted λ-hat")
        ax.grid(alpha=0.18, lw=0.5)
        ax.set_xlim(-0.25, 1.25)
        ax.set_ylim(-0.35, 0.35)
    axes[0].set_ylabel("Signed distance to assigned branch")
    fig.suptitle("Branch-aligned mixture geometry")
    save_figure(fig, f"{out_prefix}.branch_aligned_geometry", formats, args.dpi)


def plot_lambda_consistency(points: Sequence[Point], out_prefix: str, args: argparse.Namespace, formats: Sequence[str]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharex=True, sharey=True)
    rng = random.Random(args.seed + 2)
    grid = np.linspace(-0.25, 1.25, 240)
    gx, gy = np.meshgrid(grid, grid)
    absdiff = np.abs(gx - gy)
    reldiff = absdiff / (np.maximum(np.abs(gx), np.abs(gy)) + 0.01)
    accepted = (absdiff < args.lambda_abs_diff_max) & (reldiff < args.lambda_rel_diff_max)
    for ax, overlap in zip(axes, ("F1", "F2")):
        ax.contourf(gx, gy, accepted.astype(float), levels=[0.5, 1.5], colors=["#CFE8F3"], alpha=0.45)
        subset = [p for p in points if p.overlap_type == overlap and math.isfinite(p.lambda_gain) and math.isfinite(p.lambda_drop)]
        background = subset if len(subset) <= args.max_background_points // 2 else rng.sample(subset, args.max_background_points // 2)
        ax.scatter([p.lambda_gain for p in background], [p.lambda_drop for p in background],
                   s=5, c="#777777", alpha=0.10, linewidths=0, rasterized=True)
        credible = [p for p in subset if p.classification == "credible_extra_ORF_like_signal"]
        atypical = [p for p in subset if p.q_bh < args.fdr_threshold and p.classification == "atypical_pattern"]
        if atypical:
            ax.scatter([p.lambda_gain for p in atypical], [p.lambda_drop for p in atypical], s=14,
                       c=CLASS_COLORS["atypical_pattern"], alpha=0.65, linewidths=0, rasterized=True)
        if credible:
            ax.scatter([p.lambda_gain for p in credible], [p.lambda_drop for p in credible], s=16,
                       c=CLASS_COLORS["credible_extra_ORF_like_signal"], alpha=0.70, linewidths=0, rasterized=True)
        ax.plot([-0.25, 1.25], [-0.25, 1.25], color="#222222", lw=1.0)
        ax.axvline(args.lambda_min, color="#555555", lw=0.8, ls=":")
        ax.axhline(args.lambda_min, color="#555555", lw=0.8, ls=":")
        ax.set_title(overlap)
        ax.set_xlabel("λ from target-frame gain")
        ax.grid(alpha=0.15, lw=0.5)
    axes[0].set_ylabel("λ from host-frame drop")
    axes[0].set_xlim(-0.25, 1.25)
    axes[0].set_ylim(-0.25, 1.25)
    fig.suptitle("Gain/drop consistency gate (shaded region passes)")
    save_figure(fig, f"{out_prefix}.lambda_gain_vs_drop", formats, args.dpi)


def waterfall_counts(points: Sequence[Point], total: int, args: argparse.Namespace) -> List[Tuple[str, int]]:
    eligible = list(points)
    significant = [p for p in eligible if p.q_bh < args.fdr_threshold]
    direction_geometry = [p for p in significant if gate_passes(p, args)["direction"] and gate_passes(p, args)["geometry"]]
    active = [p for p in direction_geometry if gate_passes(p, args)["active breadth"]]
    residual = [p for p in active if gate_passes(p, args)["residual breadth"]]
    primary = [p for p in residual if p.classification == "credible_extra_ORF_like_signal"]
    return [
        ("All intORF_altframe", total),
        ("Coverage eligible", len(eligible)),
        (f"BH q < {args.fdr_threshold:g}", len(significant)),
        ("Direction + geometry", len(direction_geometry)),
        ("Active-codon breadth", len(active)),
        ("Target-residual breadth", len(residual)),
        ("Primary credible", len(primary)),
    ]


def plot_waterfall(points: Sequence[Point], total: int, out_prefix: str, args: argparse.Namespace, formats: Sequence[str]) -> List[Tuple[str, int]]:
    counts = waterfall_counts(points, total, args)
    labels = [x[0] for x in counts]
    values = [x[1] for x in counts]
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    y = np.arange(len(labels))
    colors = ["#777777", "#777777", "#56B4E9", "#009E73", "#009E73", "#009E73", "#0072B2"]
    bars = ax.barh(y, values, color=colors, alpha=0.85)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Candidates (log scale)")
    ax.set_title("Sequential DM evidence gates")
    ax.grid(axis="x", alpha=0.20, lw=0.5)
    for bar, value in zip(bars, values):
        ax.text(value * 1.06, bar.get_y() + bar.get_height() / 2, f"{value:,}", va="center", fontsize=9)
    save_figure(fig, f"{out_prefix}.gate_waterfall", formats, args.dpi)
    return counts


def heatmap_selection(points: Sequence[Point], args: argparse.Namespace) -> List[Point]:
    highlight_names = set(args.highlight_gene)
    highlighted = [p for p in points if p.gene_name in highlight_names]
    near = sorted(
        [p for p in points if math.isfinite(p.q_bh) and p.q_bh < 0.15 and p.gene_name not in highlight_names],
        key=lambda p: abs(p.q_bh - args.fdr_threshold),
    )
    selected = highlighted + near[: max(0, args.heatmap_rows - len(highlighted))]
    return selected[: args.heatmap_rows]


def plot_gate_heatmap(points: Sequence[Point], out_prefix: str, args: argparse.Namespace, formats: Sequence[str]) -> List[Point]:
    selected = heatmap_selection(points, args)
    if not selected:
        return []
    data = np.asarray([gate_margins(p, args) for p in selected], dtype=float)
    data = np.clip(data, -1.5, 1.5)
    columns = ["BH q", "λ", "distance", "λ abs", "λ rel", "active n", "active frac", "residual frac"]
    labels = [f"{p.gene_name or p.gene_id} | {p.overlap_type} | {p.classification}" for p in selected]
    fig_height = max(6.0, 0.27 * len(selected) + 1.8)
    fig, ax = plt.subplots(figsize=(12.5, fig_height))
    image = ax.imshow(data, aspect="auto", cmap="RdBu", vmin=-1.5, vmax=1.5)
    ax.set_xticks(np.arange(len(columns)), columns, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(labels)), labels, fontsize=7.5)
    ax.set_title("Normalized gate margins: positive passes, negative fails")
    cbar = fig.colorbar(image, ax=ax, shrink=0.65)
    cbar.set_label("Normalized margin")
    save_figure(fig, f"{out_prefix}.gate_margin_heatmap", formats, args.dpi)
    return selected


def plot_depth_effects(points: Sequence[Point], out_prefix: str, args: argparse.Namespace, formats: Sequence[str]) -> None:
    usable = [p for p in points if p.core_reads > 0 and math.isfinite(p.lambda_hat) and math.isfinite(p.q_bh)]
    rng = random.Random(args.seed + 3)
    background = usable if len(usable) <= args.max_background_points else rng.sample(usable, args.max_background_points)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.0))
    axes[0].scatter([p.core_reads for p in background], [p.lambda_hat for p in background],
                    s=5, c="#777777", alpha=0.10, linewidths=0, rasterized=True)
    axes[1].scatter([p.core_reads for p in background], [-math.log10(max(p.q_bh, 1e-300)) for p in background],
                    s=5, c="#777777", alpha=0.10, linewidths=0, rasterized=True)
    for classification in ("atypical_pattern", "localized_core_signal", "low_target_frame_residual_breadth",
                           "credible_extra_ORF_like_signal"):
        subset = [p for p in usable if p.q_bh < args.fdr_threshold and p.classification == classification]
        color = CLASS_COLORS[classification]
        axes[0].scatter([p.core_reads for p in subset], [p.lambda_hat for p in subset], s=13, c=color,
                        alpha=0.65, linewidths=0, rasterized=True)
        axes[1].scatter([p.core_reads for p in subset], [-math.log10(max(p.q_bh, 1e-300)) for p in subset],
                        s=13, c=color, alpha=0.65, linewidths=0, rasterized=True)
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel("Core P-site reads")
        ax.grid(alpha=0.18, lw=0.5)
    axes[0].axhline(args.lambda_min, color="#222222", lw=0.9, ls="--")
    axes[0].set_ylabel("λ estimate")
    axes[0].set_title("Effect size versus depth")
    axes[1].axhline(-math.log10(args.fdr_threshold), color="#222222", lw=0.9, ls="--")
    axes[1].set_ylabel("−log10(BH q)")
    axes[1].set_title("Significance versus depth")
    save_figure(fig, f"{out_prefix}.depth_effects", formats, args.dpi)


def plot_breadth_plane(
    points: Sequence[Point],
    out_prefix: str,
    args: argparse.Namespace,
    formats: Sequence[str],
) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 6.8))
    rng = random.Random(args.seed + 12)
    background = list(points)
    if len(background) > args.max_background_points:
        background = rng.sample(background, args.max_background_points)
    ax.scatter(
        [p.active_frac for p in background], [p.residual_frac for p in background],
        s=5, c="#777777", alpha=0.08, linewidths=0, rasterized=True,
    )
    for classification in (
        "atypical_pattern", "localized_core_signal",
        "low_target_frame_residual_breadth", "ultra_short_exploratory",
        "credible_extra_ORF_like_signal",
    ):
        subset = [
            p for p in points
            if p.q_bh < args.fdr_threshold and p.classification == classification
        ]
        if not subset:
            continue
        sizes = [8.0 + 5.0 * math.log10(max(p.core_reads, 1)) for p in subset]
        ax.scatter(
            [p.active_frac for p in subset], [p.residual_frac for p in subset],
            s=sizes, c=CLASS_COLORS[classification], alpha=0.68,
            linewidths=0, rasterized=True,
            label=f"{CLASS_LABELS.get(classification, classification)} (n={len(subset):,})",
        )
    ax.axvline(args.min_active_core_frac, color="#222222", lw=1.0, ls="--")
    ax.axhline(args.min_target_residual_frac, color="#222222", lw=1.0, ls="--")
    ax.text(args.min_active_core_frac + 0.01, 0.985, f"active fraction ≥ {args.min_active_core_frac:g}", va="top", fontsize=8)
    ax.text(0.985, args.min_target_residual_frac + 0.015,
            f"residual fraction > {args.min_target_residual_frac:.3g}", ha="right", fontsize=8)
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.set_xlabel("Active core codon fraction")
    ax.set_ylabel("Target-residual supported fraction among active codons")
    ax.set_title("Codon-breadth evidence gates")
    ax.grid(alpha=0.16, lw=0.5)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
    fig.subplots_adjust(right=0.78)
    save_figure(fig, f"{out_prefix}.breadth_gate_plane", formats, args.dpi)


def q_score(point: Point) -> float:
    return -math.log10(max(point.q_bh, 1e-300))


def plot_effect_significance(
    points: Sequence[Point],
    out_prefix: str,
    args: argparse.Namespace,
    formats: Sequence[str],
) -> None:
    usable = [p for p in points if math.isfinite(p.lambda_hat) and math.isfinite(p.q_bh)]
    all_scores = np.asarray([q_score(p) for p in usable])
    ymax = max(5.0, min(70.0, float(np.max(all_scores)) * 1.03))
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.2), sharex=True, sharey=True)
    rng = random.Random(args.seed + 13)
    for ax, overlap in zip(axes, ("F1", "F2")):
        subset = [p for p in usable if p.overlap_type == overlap]
        max_points = max(1000, args.max_background_points // 2)
        background = subset if len(subset) <= max_points else rng.sample(subset, max_points)
        ax.scatter(
            [p.lambda_hat for p in background], [min(q_score(p), ymax) for p in background],
            s=5, c="#777777", alpha=0.09, linewidths=0, rasterized=True,
        )
        for classification in (
            "atypical_pattern", "localized_core_signal",
            "low_target_frame_residual_breadth", "ultra_short_exploratory",
            "credible_extra_ORF_like_signal",
        ):
            selected = [
                p for p in subset
                if p.q_bh < args.fdr_threshold and p.classification == classification
            ]
            if selected:
                ax.scatter(
                    [p.lambda_hat for p in selected], [min(q_score(p), ymax) for p in selected],
                    s=14, c=CLASS_COLORS[classification], alpha=0.70,
                    linewidths=0, rasterized=True,
                )
        ax.axvline(args.lambda_min, color="#222222", lw=0.9, ls="--")
        ax.axhline(-math.log10(args.fdr_threshold), color="#222222", lw=0.9, ls="--")
        ax.set_title(overlap)
        ax.set_xlabel("Codon-level DM fitted λ-hat")
        ax.grid(alpha=0.16, lw=0.5)
    axes[0].set_ylabel("−log10(BH q)")
    axes[0].set_xlim(-0.02, 1.02)
    axes[0].set_ylim(-0.1, ymax)
    fig.suptitle("DM effect size and statistical significance")
    save_figure(fig, f"{out_prefix}.effect_significance", formats, args.dpi)


def gate_combination_rows(
    points: Sequence[Point],
    args: argparse.Namespace,
) -> List[Dict[str, object]]:
    gate_names = ["direction", "geometry", "active breadth", "residual breadth"]
    counts: Dict[Tuple[bool, ...], int] = {}
    for point in points:
        if point.q_bh >= args.fdr_threshold:
            continue
        passes = gate_passes(point, args)
        key = tuple(bool(passes[name]) for name in gate_names)
        counts[key] = counts.get(key, 0) + 1
    ordered = sorted(counts.items(), key=lambda item: (-item[1], tuple(not x for x in item[0])))
    rows: List[Dict[str, object]] = []
    for key, count in ordered:
        row: Dict[str, object] = {"count": count}
        for name, passed in zip(gate_names, key):
            row[name.replace(" ", "_")] = int(passed)
        row["combination"] = "+".join(name for name, passed in zip(gate_names, key) if passed) or "none"
        rows.append(row)
    return rows


def plot_gate_combinations(
    points: Sequence[Point],
    out_prefix: str,
    args: argparse.Namespace,
    formats: Sequence[str],
) -> List[Dict[str, object]]:
    rows = gate_combination_rows(points, args)
    shown = rows[:12]
    gate_names = ["direction", "geometry", "active_breadth", "residual_breadth"]
    display_names = ["direction", "geometry", "active breadth", "residual breadth"]
    x = np.arange(len(shown))
    fig = plt.figure(figsize=(11.2, 6.5))
    grid = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.35], hspace=0.05)
    ax_bar = fig.add_subplot(grid[0])
    ax_matrix = fig.add_subplot(grid[1], sharex=ax_bar)
    all_pass = [all(int(row[name]) == 1 for name in gate_names) for row in shown]
    colors = ["#0072B2" if passed else "#E69F00" for passed in all_pass]
    bars = ax_bar.bar(x, [int(row["count"]) for row in shown], color=colors, alpha=0.88)
    for bar, row in zip(bars, shown):
        ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height() * 1.015,
                    f"{int(row['count']):,}", ha="center", va="bottom", fontsize=8)
    ax_bar.set_ylabel("BH-significant candidates")
    ax_bar.set_title("Post-BH gate combinations")
    ax_bar.grid(axis="y", alpha=0.18, lw=0.5)
    ax_bar.tick_params(axis="x", labelbottom=False)
    for col, row in enumerate(shown):
        passed_rows = []
        for y_index, gate in enumerate(gate_names):
            passed = int(row[gate]) == 1
            ax_matrix.scatter(
                [col], [y_index], s=52 if passed else 30,
                c="#0072B2" if passed else "#D0D0D0", edgecolors="none", zorder=3,
            )
            if passed:
                passed_rows.append(y_index)
        if len(passed_rows) >= 2:
            ax_matrix.plot([col, col], [min(passed_rows), max(passed_rows)], color="#0072B2", lw=1.4, zorder=2)
    ax_matrix.set_yticks(np.arange(len(display_names)), display_names)
    ax_matrix.set_xticks(x, [str(i + 1) for i in x])
    ax_matrix.set_xlabel("Gate combination rank")
    ax_matrix.set_ylim(-0.6, len(display_names) - 0.4)
    ax_matrix.invert_yaxis()
    ax_matrix.grid(axis="x", alpha=0.12, lw=0.5)
    save_figure(fig, f"{out_prefix}.gate_combinations", formats, args.dpi)
    return rows


def length_bin_label(peptide_len: int) -> str:
    if peptide_len <= 9:
        return "6–9"
    if peptide_len <= 14:
        return "10–14"
    if peptide_len <= 29:
        return "15–29"
    if peptide_len <= 59:
        return "30–59"
    return "≥60"


def length_stratification_rows(
    points: Sequence[Point],
    args: argparse.Namespace,
) -> List[Dict[str, object]]:
    labels = ["6–9", "10–14", "15–29", "30–59", "≥60"]
    rows: List[Dict[str, object]] = []
    for label in labels:
        subset = [p for p in points if length_bin_label(p.peptide_len) == label]
        significant = [p for p in subset if p.q_bh < args.fdr_threshold]
        primary = [p for p in subset if p.classification == "credible_extra_ORF_like_signal"]
        exploratory = [p for p in subset if p.classification == "ultra_short_exploratory"]
        rows.append({
            "length_bin_aa": label,
            "eligible": len(subset),
            "bh_significant": len(significant),
            "primary_credible": len(primary),
            "ultra_short_exploratory": len(exploratory),
            "bh_significant_frac": len(significant) / len(subset) if subset else 0.0,
            "primary_credible_frac": len(primary) / len(subset) if subset else 0.0,
        })
    return rows


def plot_length_stratification(
    points: Sequence[Point],
    out_prefix: str,
    args: argparse.Namespace,
    formats: Sequence[str],
) -> List[Dict[str, object]]:
    rows = length_stratification_rows(points, args)
    plotted_rows = [row for row in rows if int(row["eligible"]) > 0]
    labels = [str(row["length_bin_aa"]) for row in plotted_rows]
    x = np.arange(len(plotted_rows))
    width = 0.24
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.2))
    axes[0].bar(x - width, [int(row["eligible"]) for row in plotted_rows], width, color="#8A8A8A", label="eligible")
    axes[0].bar(x, [int(row["bh_significant"]) for row in plotted_rows], width, color="#56B4E9", label="BH significant")
    axes[0].bar(x + width, [int(row["primary_credible"]) for row in plotted_rows], width, color="#0072B2", label="primary credible")
    axes[0].set_yscale("log")
    axes[0].set_xticks(x, labels)
    axes[0].set_xlabel("Peptide length (aa)")
    axes[0].set_ylabel("Candidates (log scale)")
    axes[0].set_title("Counts by length")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(axis="y", alpha=0.18, lw=0.5)
    axes[1].plot(x, [100.0 * float(row["bh_significant_frac"]) for row in plotted_rows], marker="o", color="#D55E00", label="BH significant / eligible")
    axes[1].plot(x, [100.0 * float(row["primary_credible_frac"]) for row in plotted_rows], marker="D", color="#0072B2", label="primary credible / eligible")
    axes[1].set_xticks(x, labels)
    axes[1].set_xlabel("Peptide length (aa)")
    axes[1].set_ylabel("Candidates (%)")
    axes[1].set_title("Call rates by length")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(alpha=0.18, lw=0.5)
    fig.suptitle("Length-stratified DM results")
    save_figure(fig, f"{out_prefix}.length_stratification", formats, args.dpi)
    return rows


def write_tsv(path: str, rows: Iterable[Mapping[str, object]], columns: Sequence[str]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in columns})


def plot_data_rows(points: Sequence[Point], args: argparse.Namespace) -> Iterable[Dict[str, object]]:
    for point in points:
        passes = gate_passes(point, args)
        yield {
            "torf_id": point.torf_id,
            "gorf_id": point.gorf_id,
            "gene_id": point.gene_id,
            "gene_name": point.gene_name,
            "overlap_type": point.overlap_type,
            "classification": point.classification,
            "gene_level_orf_class": point.gene_level_orf_class,
            "gene_level_pure_intorf_eligible": point.gene_level_pure_intorf_eligible,
            "cds_p0": f"{point.pi_obs_cds[0]:.10g}",
            "cds_plus1": f"{point.pi_obs_cds[1]:.10g}",
            "cds_plus2": f"{point.pi_obs_cds[2]:.10g}",
            "ternary_x": f"{point.x:.10g}",
            "ternary_y": f"{point.y:.10g}",
            "lambda_hat": f"{point.lambda_hat:.10g}",
            "lambda_proj": f"{point.lambda_proj:.10g}",
            "signed_branch_distance": f"{signed_branch_distance(point):.10g}",
            "q_BH": f"{point.q_bh:.10g}",
            "core_reads": point.core_reads,
            "n_active_core_codons": point.n_active,
            "active_core_codon_frac": f"{point.active_frac:.10g}",
            "target_residual_supported_frac": f"{point.residual_frac:.10g}",
            "gate_bh": int(passes["BH"]),
            "gate_direction": int(passes["direction"]),
            "gate_geometry": int(passes["geometry"]),
            "gate_active_breadth": int(passes["active breadth"]),
            "gate_residual_breadth": int(passes["residual breadth"]),
        }


def main() -> None:
    args = parse_args()
    input_path = os.path.abspath(args.dm_results)
    if not os.path.isfile(input_path) or os.path.getsize(input_path) == 0:
        raise SystemExit(f"ERROR: missing or empty DM result: {input_path}")
    formats = [x.strip().lower() for x in args.formats.split(",") if x.strip()]
    unsupported = sorted(set(formats).difference({"png", "pdf", "svg"}))
    if unsupported:
        raise SystemExit(f"ERROR: unsupported formats: {unsupported}")
    out_prefix = os.path.abspath(args.out_prefix)
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    points, stats = load_points(input_path, args.gene_level_pure_intorf_only)
    if not points:
        raise SystemExit("ERROR: no eligible F1/F2 rows with observed phase compositions")

    # Default figures are limited to summaries that an end user needs to
    # interpret one completed caller run. Developer/model diagnostics are
    # opt-in and never affect result classification.
    plot_ternary(points, out_prefix, args, formats)
    plot_ternary_by_frame(points, out_prefix, args, formats)
    candidate_total = stats["candidate_after_gene_context"]
    waterfall = plot_waterfall(points, candidate_total, out_prefix, args, formats)
    plot_breadth_plane(points, out_prefix, args, formats)
    plot_effect_significance(points, out_prefix, args, formats)
    combinations = plot_gate_combinations(points, out_prefix, args, formats)
    length_rows = plot_length_stratification(points, out_prefix, args, formats)
    plot_depth_effects(points, out_prefix, args, formats)
    selected: List[Point] = []
    if args.include_diagnostics:
        plot_branch_geometry(points, out_prefix, args, formats)
        plot_lambda_consistency(points, out_prefix, args, formats)
        selected = plot_gate_heatmap(points, out_prefix, args, formats)

    data_columns = [
        "torf_id", "gorf_id", "gene_id", "gene_name", "overlap_type", "classification",
        "gene_level_orf_class", "gene_level_pure_intorf_eligible",
        "cds_p0", "cds_plus1", "cds_plus2", "ternary_x", "ternary_y", "lambda_hat",
        "lambda_proj", "signed_branch_distance", "q_BH", "core_reads", "n_active_core_codons",
        "active_core_codon_frac", "target_residual_supported_frac", "gate_bh", "gate_direction",
        "gate_geometry", "gate_active_breadth", "gate_residual_breadth",
    ]
    write_tsv(f"{out_prefix}.plot_data.tsv", plot_data_rows(points, args), data_columns)

    waterfall_rows = [
        {"stage": stage, "count": count, "fraction_of_all": f"{count / candidate_total:.10g}"}
        for stage, count in waterfall
    ]
    write_tsv(f"{out_prefix}.gate_waterfall.tsv", waterfall_rows, ["stage", "count", "fraction_of_all"])

    write_tsv(
        f"{out_prefix}.gate_combinations.tsv",
        combinations,
        ["combination", "count", "direction", "geometry", "active_breadth", "residual_breadth"],
    )
    write_tsv(
        f"{out_prefix}.length_stratification.tsv",
        length_rows,
        ["length_bin_aa", "eligible", "bh_significant", "primary_credible",
         "ultra_short_exploratory", "bh_significant_frac", "primary_credible_frac"],
    )

    selected_rows = []
    for point in selected:
        selected_rows.append({
            "gene_name": point.gene_name,
            "gene_id": point.gene_id,
            "gorf_id": point.gorf_id,
            "torf_id": point.torf_id,
            "overlap_type": point.overlap_type,
            "classification": point.classification,
            "q_BH": f"{point.q_bh:.10g}",
            "lambda_hat": f"{point.lambda_hat:.10g}",
            "distance_to_segment": f"{point.distance:.10g}",
            "active_core_codons": point.n_active,
            "active_core_frac": f"{point.active_frac:.10g}",
            "target_residual_supported_frac": f"{point.residual_frac:.10g}",
        })
    write_tsv(
        f"{out_prefix}.gate_margin_candidates.tsv", selected_rows,
        ["gene_name", "gene_id", "gorf_id", "torf_id", "overlap_type", "classification", "q_BH",
         "lambda_hat", "distance_to_segment", "active_core_codons", "active_core_frac",
         "target_residual_supported_frac"],
    )

    manifest = [
        ("program", PROGRAM),
        ("version", VERSION),
        ("input_dm_results", input_path),
        ("candidate_total_raw", stats["candidate_total_raw"]),
        ("candidate_after_gene_context", stats["candidate_after_gene_context"]),
        ("gene_context_excluded_rows", stats["gene_context_excluded_rows"]),
        ("gene_level_pure_intorf_only", int(args.gene_level_pure_intorf_only)),
        ("eligible_geometry_rows", len(points)),
        ("fdr_threshold", args.fdr_threshold),
        ("lambda_min", args.lambda_min),
        ("lambda_abs_diff_max", args.lambda_abs_diff_max),
        ("lambda_rel_diff_max", args.lambda_rel_diff_max),
        ("distance_to_segment_max", args.distance_to_segment_max),
        ("min_active_core_codons", args.min_active_core_codons),
        ("min_active_core_frac", args.min_active_core_frac),
        ("min_target_residual_frac_strictly_greater_than", args.min_target_residual_frac),
        ("coordinate_contract_F1", "CDS_P0=pi_obs_2;CDS_+1=pi_obs_0;CDS_+2=pi_obs_1"),
        ("coordinate_contract_F2", "CDS_P0=pi_obs_1;CDS_+1=pi_obs_2;CDS_+2=pi_obs_0"),
        ("highlight_genes", ",".join(args.highlight_gene)),
        ("include_diagnostics", int(args.include_diagnostics)),
        ("formats", ",".join(formats)),
    ]
    write_tsv(f"{out_prefix}.manifest.tsv", ({"key": k, "value": v} for k, v in manifest), ["key", "value"])

    print(f"Loaded candidates: {stats['candidate_total_raw']:,}")
    if args.gene_level_pure_intorf_only:
        print(f"Excluded by gene-level CDS context: {stats['gene_context_excluded_rows']:,}")
        print(f"Candidates after gene-level context: {stats['candidate_after_gene_context']:,}")
    print(f"Eligible geometry rows: {len(points):,}")
    print(f"Primary credible: {sum(p.classification == 'credible_extra_ORF_like_signal' for p in points):,}")
    print(f"Wrote visualization prefix: {out_prefix}")


if __name__ == "__main__":
    main()
