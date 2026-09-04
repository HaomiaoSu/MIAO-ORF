#!/usr/bin/env bash
set -euo pipefail

# Formal DM caller. Accurate is the default; set DM_MODE=fast for the validated
# lower-cost Monte Carlo allocation. Both modes use the same model and gates.
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RESULTS_ROOT="${MIAO_ORF_RESULTS_ROOT:-/path/to/miao-orf-results}"
CALLER="$PROJECT/src/miao_orf/ribo_intorf_dm_caller.py"
PSITE="${MIAO_ORF_PSITE_BAM:-$RESULTS_ROOT/02_psite/ES-Ribo-Rep1/ES-Ribo-Rep1.psite.bam}"
TORF="$RESULTS_ROOT/01_orfscan/gencode_v48.torf.tsv"
BACKGROUND="$RESULTS_ROOT/03_metagene_qc/ES-Ribo-Rep1/ES-Ribo-Rep1.dm_background.tsv"
DM_MODE="${DM_MODE:-accurate}"

case "$DM_MODE" in
  accurate|fast) ;;
  *)
    echo "ERROR: DM_MODE must be accurate or fast, got: $DM_MODE" >&2
    exit 1
    ;;
esac

OUT_DIR="$RESULTS_ROOT/04_intorf_dm/ES-Ribo-Rep1/$DM_MODE"
OUT_PREFIX="$OUT_DIR/ES-Ribo-Rep1"
mkdir -p "$OUT_DIR"

for required in "$CALLER" "$PSITE" "$PSITE.bai" "$TORF" "$BACKGROUND"; do
  [[ -s "$required" ]] || {
    echo "ERROR: missing or empty required input: $required" >&2
    exit 1
  }
done

time conda run --no-capture-output -n bio python "$CALLER" \
  --psite-bam "$PSITE" \
  --torf "$TORF" \
  --dm-background "$BACKGROUND" \
  --out-prefix "$OUT_PREFIX" \
  --pi-method codon_equal \
  --min-a0 1.0 \
  --min-intorf-aa 6 \
  --primary-min-intorf-aa 10 \
  --min-core-codons 5 \
  --min-active-core-codons 3 \
  --min-credible-active-core-codons 5 \
  --min-credible-active-core-frac 0.15 \
  --min-credible-target-residual-frac 1/3 \
  --min-core-reads 15 \
  --exclude-start-codons 1 \
  --exclude-stop-codons 1 \
  --bootstrap-gate-p 0.20 \
  --bootstrap-reps 999 \
  --bootstrap-engine adaptive_importance \
  --importance-mode "$DM_MODE" \
  --importance-reps 9999 \
  --importance-iid-exceedance-threshold 10 \
  --importance-etas 0,0.05,0.1,0.2,0.4,0.7,1.0 \
  --importance-min-tail-ess 30 \
  --importance-max-relative-se 0.25 \
  --disable-block-bootstrap \
  --seed 20260821 \
  --fdr-threshold 0.05 \
  --lambda-min 0.05 \
  --template-separation-min 0.05 \
  --lambda-grid-size 101 \
  --lambda-abs-diff-max 0.10 \
  --lambda-rel-diff-max 0.30 \
  --lambda-rel-eps 0.01 \
  --distance-to-segment-max 0.10 \
  --review-gate-window 0.05 \
  --candidate-dedup gorf \
  --require-gorf-id \
  --member-preview-n 20 \
  --workers 8 \
  --mp-chunksize 1 \
  2>&1 | tee "$OUT_DIR/run.log"

echo "DONE"
echo "  mode: $DM_MODE"
echo "  result: $OUT_PREFIX.intorf_dm_results.tsv"
