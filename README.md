# MIAO

**MIAO** (**M**ixture-aware **I**nference of **A**lternative-frame **O**RFs) detects translated alternative-frame internal open reading frames (intORFs) from ribosome-profiling data. Its Dirichlet–multinomial model tests whether codon-level phase counts are better explained by host-CDS translation alone or by a mixture of host and alternative-frame translation, and estimates the intORF-like mixture fraction, λ.

The project uses distinct names in different contexts:

- method and manuscript: **MIAO**;
- GitHub repository and release: **MIAO-ORF**;
- Python distribution and command: `miao-orf`;
- Python import package: `miao_orf`.

Version **1.0.0** is the first public release. Its formal statistical settings are recorded in [`config/releases/benchmark_520_dual_v1.0.0.json`](config/releases/benchmark_520_dual_v1.0.0.json).

## What MIAO provides

- transcriptome-wide candidate ORF scanning;
- one-nucleotide P-site BAM generation;
- codon-level metagene QC and sample-specific phase calibration;
- host-only versus host-plus-intORF Dirichlet–multinomial inference;
- fitted λ, false-discovery-rate control and evidence gates;
- model-allocated intORF abundance in pFPKM;
- post-processing for same-gene annotated-CDS N-terminal reuse;
- replicate-aware summaries and publication-oriented visualizations;
- an optional codon-level P-site exporter for selected intORFs.

MIAO preserves complete raw DM and abundance tables. Final count tables and figures use the gene-level pure-intORF eligibility annotation without changing the transcript-level ORF type or the original DM statistics.

## Requirements

- Linux or WSL2;
- Python 3.10–3.12;
- `samtools` 1.18 or later;
- an indexed genome FASTA, matching GTF annotation and coordinate-sorted/indexed Ribo-seq BAM.

## Installation

Create the reproducible Conda environment:

```bash
git clone https://github.com/HaomiaoSu/MIAO-ORF.git
cd MIAO-ORF
conda env create -f environment.yml
conda activate miao-orf-1.0.0
miao-orf --version
```

Alternatively, install into an existing environment that already supplies `samtools>=1.18`:

```bash
python -m pip install .
miao-orf --version
```

## P-site offsets

MIAO accepts either explicit length-specific P-site offsets or offsets estimated by Ribo-TISH 0.2.8.

For explicit offsets:

```bash
miao-orf \
  --sample SAMPLE \
  --out-root results \
  --gtf annotation.gtf \
  --fa genome.fa \
  --bam SAMPLE.bam \
  --length-offsets 28:12,29:12,30:12
```

For Ribo-TISH-assisted selection, first run `ribotish quality` and provide its `mapped.para.py` and `mapped_qual.txt` outputs:

```bash
miao-orf \
  --sample SAMPLE \
  --out-root results \
  --gtf annotation.gtf \
  --fa genome.fa \
  --bam SAMPLE.bam \
  --offsets SAMPLE_mapped.para.py \
  --ribotish-quality SAMPLE_mapped_qual.txt
```

MIAO does not re-estimate offsets. It imports Ribo-TISH's principal offsets, re-indexes the exact 5′-matched phase counts from Ribo-TISH's fixed 12-nt QC anchor to each length-specific offset, retains lengths with corrected F0 greater than two thirds and, when needed, selects the contiguous passing length block with the greatest matched-read support. The Ribo-TISH `m0` mismatch group is excluded. Every run writes an auditable offset-selection table before placing one P-site per retained footprint.

The standalone audit command is:

```bash
miao-orf-offsets --help
```

## Pipeline stages

The unified command can run the full workflow or any contiguous subset:

1. `orfscan` — scan the annotated transcriptome and generate reference candidate tables;
2. `psite` — construct a one-nucleotide P-site BAM;
3. `qc` — generate metagene and codon-level QC and estimate the sample background;
4. `dm` — perform host-only versus mixture inference;
5. `abundance` — calculate model-allocated intORF pFPKM;
6. `context` — annotate same-gene CDS N-terminal reuse;
7. `codon` — optionally export selected-intORF codon counts;
8. `visualize` — generate final sample-level figures and summaries.

Use `--from-stage`, `--to-stage` or `--only-stage` to control execution. Completed stages have content-validated completion records and can be resumed safely. Run `miao-orf --help` for all options.

Replicate batches analyze each BAM independently, reuse the reference scan and then integrate evidence by `gorf_id + overlap_type`. Reads are not pooled and replicate P values are not combined. See [`docs/MANUAL.md`](docs/MANUAL.md) for the complete interface and [`docs/MIAO_ORF_PIPELINE.md`](docs/MIAO_ORF_PIPELINE.md) for a concise workflow reference.

## Statistical interpretation

The estimated λ is the fitted intORF-like fraction of analyzed overlap-core P-sites, not a hard assignment of individual reads. MIAO calls are therefore evidence for an alternative-frame translated component within a host CDS. They do not by themselves establish protein stability, biological function or the exact molecular mechanism of initiation.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Please cite both the software release and the associated MIAO publication when available.

## License

MIAO-ORF is released under the MIT License. See [`LICENSE`](LICENSE).

## Author

Haomiao Su — Molecular Biology Research Center, School of Life Sciences, Hunan Province Key Laboratory of Basic and Applied Hematology, Central South University, Hunan, China. Contact: `suhaomiao@csu.edu.cn`.
