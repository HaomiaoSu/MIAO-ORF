#!/usr/bin/env bash
set -euo pipefail
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA_ROOT="${MIAO_ORF_DATA_ROOT:-/path/to/miao-orf-data}"
RESULTS_ROOT="${MIAO_ORF_RESULTS_ROOT:-/path/to/miao-orf-results}"

# Step 02: convert the mapped Ribo-seq BAM into a 1-nt P-site BAM.
mkdir -p "$RESULTS_ROOT/02_psite/ES-Ribo-Rep1"

python3 "$PROJECT/src/miao_orf/psite-caller.py" \
  --bam "$DATA_ROOT/Ribo/ES-Ribo-Rep1.mapped.bam" \
  --offsets "$DATA_ROOT/Ribo/ES-Ribo-Rep1.mapped.para.py" \
  --ribotish-quality "$DATA_ROOT/Ribo/ES-Ribo-Rep1.mapped_qual.txt" \
  --out-prefix "$RESULTS_ROOT/02_psite/ES-Ribo-Rep1/ES-Ribo-Rep1" \
  --min-mapq 20 \
  --require-unique \
  --workers 8 \
  --bgzf-threads 2 \
  --merge \
  --no-bed \
  --no-bedgraph
