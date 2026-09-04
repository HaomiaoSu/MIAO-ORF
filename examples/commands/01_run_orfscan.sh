#!/usr/bin/env bash
set -euo pipefail
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA_ROOT="${MIAO_ORF_DATA_ROOT:-/path/to/miao-orf-data}"
RESULTS_ROOT="${MIAO_ORF_RESULTS_ROOT:-/path/to/miao-orf-results}"

# Step 01: scan GENCODE v48 transcripts and build tORF/gORF outputs.
mkdir -p "$RESULTS_ROOT/01_orfscan/per_chrom"

python3 "$PROJECT/src/miao_orf/orf_scan_transcriptome.py" \
  --gtf "$DATA_ROOT/database/gencode.v48.annotation.gtf" \
  --fa "$DATA_ROOT/database/pri_hg38.fa" \
  --out-prefix "$RESULTS_ROOT/01_orfscan/gencode_v48" \
  --primary-only \
  --by-chrom \
  --perchrom-outdir "$RESULTS_ROOT/01_orfscan/per_chrom" \
  --min-aa 6 \
  --workers 8 \
  --mp-chunksize 50
