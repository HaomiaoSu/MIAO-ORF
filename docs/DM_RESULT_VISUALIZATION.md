# DM result visualization

`src/miao_orf/visualize_intorf_dm_results.py` is a read-only reporting module.
It does not refit the DM model, recompute p-values, or modify caller output.

## Coordinate contract

Every observed phase composition is rotated into annotated-CDS coordinates:

- F1: `(CDS P0, CDS +1, CDS +2) = (pi_obs_2, pi_obs_0, pi_obs_1)`
- F2: `(CDS P0, CDS +1, CDS +2) = (pi_obs_1, pi_obs_2, pi_obs_0)`

The annotated-CDS host template is therefore shared by F1 and F2. The two
alternative translation endpoints are cyclic permutations of the same
translation template. Since `pi(lambda)=(1-lambda)H+lambda*T`, the expected
trajectories are straight line segments in the ternary simplex.

## Default end-user outputs

- `model_geometry_ternary`: all eligible candidates in CDS coordinates, the
  F1/F2 H-to-T segments, lambda landmarks, and the distance corridor.
- `model_geometry_ternary_by_frame`: separate F1/F2 views so the smaller F2
  candidate set is not hidden by the larger F1 set.
- `gate_waterfall`: sequential counts from all candidates to primary credible.
- `gate_combinations`: common combinations of post-BH gate passes/failures.
- `breadth_gate_plane`: active-codon breadth versus target-residual breadth.
- `effect_significance`: fitted lambda versus BH significance, split by frame.
- `length_stratification`: eligible, significant, and credible results by
  peptide-length bin.
- `depth_effects`: lambda and BH significance versus core P-site depth.
- `plot_data.tsv`: auditable transformed coordinates and gate flags.
- `manifest.tsv`: input path, thresholds, and coordinate definitions.

These are run-level result summaries intended for software users. They do not
include article-specific highlighted genes or benchmark-only comparisons.

## Optional diagnostics

With `--include-diagnostics`, the module additionally writes:

- `branch_aligned_geometry`: fitted codon-level DM lambda versus signed
  distance from the assigned segment.
- `lambda_gain_vs_drop`: direct view of gain/drop consistency.
- `gate_margin_heatmap`: candidates nearest the BH boundary, plus any explicit
  `--highlight-gene` values.

These diagnostics are intended for model development or detailed review and
are not produced by the default software command.

The main ternary plot can display geometry gates directly. BH significance,
active-codon breadth, and target-residual breadth are not ternary-coordinate
boundaries and are represented by point classes or companion plots.

The ternary point is the pooled observed phase composition, whereas
`lambda_hat` is the maximum-likelihood estimate from all per-codon
Dirichlet-multinomial count vectors. Therefore the orthogonal ternary projection
(`lambda_proj`) and fitted `lambda_hat` are not interchangeable. Mixture ticks
on the ternary segment are reference compositions, not the `lambda_hat > 0.05`
decision boundary.

## ES-Ribo-Rep1 command

```bash
bash examples/commands/visualization/01_visualize_formal_dm_ES-Ribo-Rep1.sh
```
