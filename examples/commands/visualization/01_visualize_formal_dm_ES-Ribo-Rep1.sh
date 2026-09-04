#!/usr/bin/env bash
set -euo pipefail

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
RESULTS_ROOT="${MIAO_ORF_RESULTS_ROOT:-/path/to/miao-orf-results}"
INPUT="$RESULTS_ROOT/04_intorf_dm/ES-Ribo-Rep1/accurate/ES-Ribo-Rep1.intorf_dm_results.tsv"
OUT_DIR="$RESULTS_ROOT/05_visualization/ES-Ribo-Rep1/accurate"
OUT_PREFIX="$OUT_DIR/ES-Ribo-Rep1.formal_adaptive_importance"
CONDA=${CONDA_EXE:-conda}

mkdir -p "$OUT_DIR"

"$CONDA" run --no-capture-output -n bio \
  python "$PROJECT/src/miao_orf/visualize_intorf_dm_results.py" \
    --dm-results "$INPUT" \
    --out-prefix "$OUT_PREFIX" \
    --fdr-threshold 0.05 \
    --lambda-min 0.05 \
    --lambda-abs-diff-max 0.10 \
    --lambda-rel-diff-max 0.30 \
    --distance-to-segment-max 0.10 \
    --min-active-core-codons 5 \
    --min-active-core-frac 0.15 \
    --min-target-residual-frac 1/3 \
    --formats png,pdf \
    --dpi 220 \
  2>&1 | tee "$OUT_DIR/run.log"
