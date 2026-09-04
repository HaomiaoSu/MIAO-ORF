#!/usr/bin/env bash
set -euo pipefail
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RESULTS_ROOT="${MIAO_ORF_RESULTS_ROOT:-/path/to/miao-orf-results}"

# Step 03: run pooled-length metagene QC and estimate the DM background.
# --rl is intentionally omitted so all retained lengths (27-30 nt) are pooled.
mkdir -p "$RESULTS_ROOT/03_metagene_qc/ES-Ribo-Rep1"

python3 "$PROJECT/src/miao_orf/ribo_metagene_qc.py" \
  --psite-bam "$RESULTS_ROOT/02_psite/ES-Ribo-Rep1/ES-Ribo-Rep1.psite.bam" \
  --torf "$RESULTS_ROOT/01_orfscan/gencode_v48.torf.tsv" \
  --out-prefix "$RESULTS_ROOT/03_metagene_qc/ES-Ribo-Rep1/ES-Ribo-Rep1" \
  --pi-method codon_equal \
  --min-a0 1.0 \
  --workers 8 \
  2>&1 | tee "$RESULTS_ROOT/03_metagene_qc/ES-Ribo-Rep1/run.log"
