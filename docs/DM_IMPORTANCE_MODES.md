# DM adaptive-importance modes

`ribo_intorf_dm_caller.py` provides two formal compute modes for resolving the
small-p tail after the ordinary IID bootstrap becomes too coarse:

| Mode | CLI | Calculation | Recommended use |
|---|---|---|---|
| Accurate | `--importance-mode accurate` | Fixed `--importance-reps` importance sample | Final analyses, stable continuous ranking and ROC |
| Fast | `--importance-mode fast` | Independent pilot allocates a separate confirmation sample | Screening and routine runs where final calls matter more than fine p-value ranking |

Both modes use the same biological DM model, host/translation templates, LRT,
importance proposals, diagnostic gates and FDR procedure. They differ only in
Monte Carlo allocation. The fast-mode pilot is not included in the reported p
value; it only determines the size of an independent confirmation sample. This
avoids treating a naively optional-stopped running mean as a fixed-sample p
estimate.

## Recommended commands

Accurate mode is the default once the formal adaptive engine is selected:

```bash
python ribo_intorf_dm_caller.py \
  ... \
  --bootstrap-engine adaptive_importance \
  --importance-mode accurate \
  --importance-reps 9999 \
  --importance-iid-exceedance-threshold 10 \
  --importance-min-tail-ess 30 \
  --importance-max-relative-se 0.25
```

Fast mode:

```bash
python ribo_intorf_dm_caller.py \
  ... \
  --bootstrap-engine adaptive_importance \
  --importance-mode fast \
  --importance-reps 9999 \
  --importance-pilot-reps 999 \
  --importance-confirm-min-reps 999 \
  --importance-allocation-safety-factor 1.25 \
  --importance-iid-exceedance-threshold 10 \
  --importance-min-tail-ess 30 \
  --importance-max-relative-se 0.25
```

For the bundled ES-Ribo command, the same selection is available without
editing the script:

```bash
# Accurate is the default.
bash examples/commands/04_run_intorf_dm_ES-Ribo-Rep1.sh

# Fast mode writes to a separate output directory.
DM_MODE=fast \
  bash examples/commands/04_run_intorf_dm_ES-Ribo-Rep1.sh
```

Selecting the formal `adaptive_importance` engine automatically allows a
diagnostically acceptable importance estimate to contribute to `p_final`.
`--importance-use-for-p-final` remains accepted for old trial commands but is
not required by the formal engine.

## Output and fallback

The result table records:

- `importance_compute_mode`;
- IID and importance p values;
- Monte Carlo SE, relative SE and tail ESS;
- pilot, confirmation and total replicate counts in fast mode;
- diagnostic status and whether importance sampling supplied `p_final`.

If the importance estimate fails the configured ESS, relative-error or weight
diagnostics, the caller falls back to the IID component instead of silently
using an unstable tail estimate.

## ES-Ribo-Rep1 benchmark

The formal comparison used the same 500 positive and 500 negative primary
spike-in targets in original and clean backgrounds. Twenty additional R0.1
pairs were supplementary and excluded from the primary ROC.

| Background | Mode | ROC AUC | Average precision | TP | FP |
|---|---|---:|---:|---:|---:|
| Original | Accurate | 0.965 | 0.966 | 471 | 35 |
| Original | Fast | 0.956 | 0.956 | 471 | 35 |
| Clean | Accurate | 0.922 | 0.930 | 430 | 51 |
| Clean | Fast | 0.902 | 0.907 | 430 | 51 |

The two modes produced identical default calls in this benchmark. Fast mode
used about 79% fewer importance replicates in the separate 200-candidate
allocation validation, but its smaller confirmation samples added continuous
p-value ranking noise, especially in the clean background. Therefore accurate
mode remains the default and should be used for final ranking-sensitive work.

Benchmark outputs are under:

`05_spikein/ES-Ribo-Rep1/source_informed/benchmark_520_dual/fixed_vs_two_stage_importance/`

## Backward compatibility

The legacy engine names `adaptive_importance_trial`,
`adaptive_two_stage_importance_trial` and `iid_plus_importance_trial` remain
accepted. New commands should use `adaptive_importance` with an explicit mode.
