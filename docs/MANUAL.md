# MIAO User Manual

## 1. Purpose

MIAO identifies candidate internal open reading frames that overlap an
annotated coding sequence in an alternative frame (`intORF_altframe`). It uses
Ribo-seq P-site data and a Dirichlet-multinomial (DM) mixture model to test
whether an overlap region contains an additional ORF-like translation
component beyond the annotated CDS signal.

MIAO abbreviates **Mixture-aware Inference of Alternative-frame ORFs**.

The default end-user workflow has seven stages:

```text
Transcriptome ORF scan
        |
        v
Mapped Ribo-seq BAM -> 1-nt P-site BAM
        |
        v
Metagene QC and DM background estimation
        |
        v
intORF DM caller
        |
        v
Model-allocated intORF abundance (pFPKM)
        |
        v
Gene-level annotated-CDS N-terminal context
        |
        v
Result visualization
```

The repository-level entry point is `miao_orf.py`. It orchestrates the seven
default programs, validates inputs and outputs, supports resuming, and records
the exact commands. An eighth `codon` stage is an optional internal ORF P-site
exporter and is disabled by default. The runner does not duplicate scientific
logic.

This manual covers the formal end-user workflow. Benchmark and manuscript
reproducibility materials are distributed separately from the user software.

## 2. Current scope

The current release is designed for:

- a GENCODE/Ensembl-style GTF;
- a coordinate-sorted and indexed genome-aligned Ribo-seq BAM;
- either explicit read-length/P-site-offset pairs or matching Ribo-TISH
  `mapped_qual.txt` and `mapped.para.py` files;
- analysis focused on `intORF_altframe` candidates;
- formal candidate deduplication by `gorf_id`;
- exclusion of annotated-stop-confounded candidates by default;
- analyses without TI-seq;
- a default P-site subset selected by corrected frame-0 proportion > 2/3,
  followed by the read-richest contiguous passing-length block, with MAPQ >= 20,
  uniquely mapped, 5'-matched reads.

MIAO does not currently infer P-site offsets. An internal inference mode
may be added later, but the formal workflow currently requires a Ribo-TISH
quality/offset pair.

## 3. Runtime environment

### 3.1 Recommended platform

Linux or WSL2 is recommended. P-site BAM merging, sorting, and indexing require
`samtools`.

### 3.2 Dependencies

- Python >= 3.10
- samtools
- pyfaidx
- pysam
- numpy
- scipy
- matplotlib

Example conda environment:

```bash
conda create -n miao-orf python=3.12
conda activate miao-orf
conda install -c conda-forge -c bioconda \
  numpy scipy matplotlib pysam pyfaidx samtools
```

If the existing project environment is available:

```bash
conda activate bio
```

Verify the environment:

```bash
python3 --version
samtools --version
python3 -c "import pyfaidx, pysam, numpy, scipy, matplotlib; print('dependencies OK')"
```

## 4. Input requirements

### 4.1 GTF annotation

The annotation should:

- use GENCODE/Ensembl-style attributes;
- provide `gene_id` and `transcript_id` for exon and CDS records;
- use the same genome build and contig names as the FASTA and BAM;
- preferably be a complete GENCODE annotation GTF.

Example:

```text
/path/to/miao-orf-data/database/gencode.v48.annotation.gtf
```

### 4.2 Genome FASTA

The FASTA must match the GTF genome build and have a `.fai` index:

```bash
samtools faidx /path/to/genome.fa
```

Example:

```text
/path/to/miao-orf-data/database/pri_hg38.fa
/path/to/miao-orf-data/database/pri_hg38.fa.fai
```

### 4.3 Mapped Ribo-seq BAM

The BAM must be:

- coordinate-sorted;
- indexed;
- aligned to the same genome as the GTF and FASTA;
- annotated with MD tags on aligned reads.

Check and index the BAM:

```bash
samtools quickcheck -v sample.mapped.bam
samtools index sample.mapped.bam
```

If MD tags are missing, use the matching genome FASTA:

```bash
samtools calmd -b sample.mapped.bam genome.fa > sample.with_md.bam
samtools index sample.with_md.bam
```

### 4.4 P-site length and offset inputs

Choose exactly one of two formal modes.

**Traditional explicit mode** supplies the accepted read lengths and offsets
directly and does not require Ribo-TISH files:

```bash
--length-offsets 28:12,29:12,30:12
```

**Automatic Ribo-TISH mode** supplies the matching `mapped_qual.txt` and
`mapped.para.py` generated for the same BAM. The P-site caller reads top-level
integer entries from `para.py`:

```python
offdict = {
    27: 11,
    28: 12,
    29: 12,
    30: 13,
    "m0": {...},
}
```

Only the main top-level offsets are used. In `mapped_qual.txt`, Ribo-TISH writes
five dictionaries for reads whose aligned 5' nucleotide matches the reference,
then five dictionaries for its `m0` (5'-mismatch) group. MIAO safely
parses the first group, obtains frame counts from its fourth dictionary, and
does not use the nested `m0` offsets or the second quality group.

Ribo-TISH computes these QC frame counts using a fixed 12-nt anchor. Before
screening, MIAO rotates them to the frame implied by each
length-specific main offset. The default first retains lengths whose corrected
P0 / (P0 + P1 + P2) is strictly greater than 2/3, then keeps the contiguous
passing-length block with the largest total number of matched reads. This avoids
admitting isolated periodic lengths when a much better-supported neighboring
length block is present.

The companion `mapped_qual.pdf` is useful for human visual QC only. It is not
an MIAO input and is never parsed: plotted proportions are rounded (for
example, 27 nt may display as 0.60 even when its underlying proportion is
0.596203), while selection must use the exact counts in `mapped_qual.txt`.

### 4.5 Existing intermediate files

When starting after an upstream stage, provide the corresponding file:

- `--torf`: an existing `*.torf.tsv`;
- `--psite-bam`: an existing indexed 1-nt P-site BAM;
- `--dm-background`: a validated `*.dm_background.tsv`;
- `--dm-results`: an existing `*.intorf_dm_results.tsv`.

## 5. Quick start

Enter the repository and activate the environment:

```bash
cd /path/to/MIAO
conda activate bio
```

`--out-root` is the shared analysis root, not a sample directory. Do not append
the value of `--sample`. Reference-level ORFscan outputs are written once under
`01_orfscan/`, while sample-dependent outputs acquire their own sample
subdirectory in later stages.

First print and validate the planned commands:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --gtf /path/to/miao-orf-data/database/gencode.v48.annotation.gtf \
  --fa /path/to/miao-orf-data/database/pri_hg38.fa \
  --bam /path/to/miao-orf-data/Ribo/ES-Ribo-Rep1.mapped.bam \
  --offsets /path/to/miao-orf-data/Ribo/ES-Ribo-Rep1.mapped.para.py \
  --ribotish-quality /path/to/miao-orf-data/Ribo/ES-Ribo-Rep1.mapped_qual.txt \
  --workers 8 \
  --dry-run
```

After checking the paths and parameters, remove `--dry-run` to execute the
workflow:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --gtf /path/to/miao-orf-data/database/gencode.v48.annotation.gtf \
  --fa /path/to/miao-orf-data/database/pri_hg38.fa \
  --bam /path/to/miao-orf-data/Ribo/ES-Ribo-Rep1.mapped.bam \
  --offsets /path/to/miao-orf-data/Ribo/ES-Ribo-Rep1.mapped.para.py \
  --ribotish-quality /path/to/miao-orf-data/Ribo/ES-Ribo-Rep1.mapped_qual.txt \
  --workers 8
```

To produce both formal DM inference modes while reusing the same ORF scan,
P-site BAM, and QC background, add:

```bash
--dm-mode both
```

When shell activation is unavailable (for example, from a non-interactive
launcher), invoke the complete environment explicitly so that both Python
dependencies and `samtools` are on `PATH`:

```bash
conda run --no-capture-output -n bio python miao_orf.py ...
```

### 5.1 Replicate batch mode

Use `--batch` to process biological replicates independently while sharing the
reference-level ORFscan output. The TSV must contain these six tab-separated
columns:

| replicate_id | bam | input_mode | length_offsets | ribotish_para | ribotish_quality |
|---|---|---|---|---|---|
| Rep1 | data/Rep1.bam | explicit | 28:12,29:12,30:12 |  |  |
| Rep2 | data/Rep2.bam | ribotish |  | qc/Rep2.mapped.para.py | qc/Rep2.mapped_qual.txt |

Paths that are not absolute are resolved relative to the batch TSV. Each row
must choose `explicit` or `ribotish`. Do not populate fields belonging to the
other mode.

Run and first validate the complete batch with:

```bash
python3 miao_orf.py \
  --batch replicates.tsv \
  --out-root /path/to/miao-orf-results/pipeline \
  --gtf /path/to/annotation.gtf \
  --fa /path/to/genome.fa \
  --workers 8 \
  --dry-run
```

Remove `--dry-run` to execute. Replicates run sequentially so one sample cannot
oversubscribe the machine alongside another. ORFscan is generated once and is
then reused; BAMs and read counts are never pooled. A failed replicate is
recorded without preventing later rows from running, and rerunning the same
command skips outputs already complete.

Completed batches write timestamped files under `OUT_ROOT/batch_runs/`:

- `*.batch_summary.tsv`: replicate/mode status, selected offsets, P-site
  alignment count, QC frame-0 proportion and A0, candidate count, primary
  credible-call count, output paths and errors;
- `*.batch_run.json`: batch provenance and the same per-replicate records.

When the selected range includes `dm`, the runner then calls the independent
`src/miao_orf/integrate_replicates.py` module once per DM mode. Outputs are
written under `OUT_ROOT/06_replicate_integration/MODE/`:

- `*.replicate_long.tsv`: one observed candidate/replicate row retaining the
  original classification, q values, lambda and core-read evidence;
- `*.consensus.tsv`: one row per union candidate plus per-replicate evidence
  columns and transparent support counts;
- `*.summary.tsv` and `*.manifest.json`: policy, replicate availability,
  category counts, inputs and outputs.

After successful integration, the runner automatically calls the independent
`src/miao_orf/visualize_replicate_integration.py` module. It writes PNG and
PDF figures beside the integration tables (controlled by the existing
`--formats` and `--dpi` options):

- `*.call_counts.*`: observed, BH-significant and primary-credible candidate
  counts for every replicate;
- `*.support_combinations.*`: UpSet-style primary-call support combinations;
- `*.lambda_concordance.*` and `*.lambda_concordance_all_eligible.*`:
  pairwise replicate lambda estimates arranged as dynamic lower-triangle
  matrices, with per-replicate lambda distributions on the diagonal. The
  former retains the primary-union/both-versus-one-primary distinction; the
  latter shows every candidate observed with finite lambda in both samples.
  The layout follows the actual replicate count rather than assuming four
  samples;
- `*.lambda_heatmap.*` and `*.lambda_matrix.tsv`: absolute `lambda_hat` values
  (fixed 0--1 scale) for candidates primary credible in any replicate, with
  missing observations shown separately;
- `*.primary_reproducibility.*`: candidate-by-replicate primary-call heatmap,
  with unavailable replicates distinguished from negative observations;
- when abundance tables are available, `*.pfpkm_correlation.*`,
  `*.pfpkm_correlation_matrix.*`, and `*.pfpkm_correlations.tsv`: pairwise
  Pearson/Spearman concordance on `log2(1 + intORF_pFPKM)` with the scatter
  panels arranged as a lower-triangle matrix, per-replicate distributions on
  the diagonal, and the shared candidate count retained for every comparison;
- when abundance tables are available, `*.pfpkm_heatmap.*` and
  `*.pfpkm_matrix.tsv`: absolute abundance for the primary-union candidate set;
  missing/non-quantifiable values remain missing rather than being replaced by
  zero. No within-candidate row-z-score panel is produced because replicate
  batches are interpreted as within-group measurements;
- `*.plot_summary.tsv` and `*.plot_manifest.json`: plotted counts, parameters,
  source tables and figure paths.

The figures are descriptive QC. They do not pool reads, combine p values,
change calls or replace the evidence-preserving TSV outputs.

Integration uses `gorf_id + overlap_type`, matching formal DM deduplication. It
does not pool reads, refit the model or combine p values. By default a candidate
is `reproducible_primary_credible` only when it is a primary credible call in at
least two available replicate results and in at least half of those results.
Configure these descriptive consensus thresholds with:

```bash
--integration-min-replicates 2 --integration-min-fraction 0.5
```

If only one replicate result is available, a primary call is explicitly marked
`single_replicate_primary_credible`, never reproducible. Failed or unavailable
replicates are retained in the manifest and are not silently treated as
negative observations.

Integration is intentionally strict about provenance. Every available result
must have its matching `*.run_manifest.json`. The integration gate requires the
same software and schema versions, statistical engine, benchmark certification,
caller source hash, semantic DM parameters, tORF input identities and output
schema. It also verifies that each result still matches the size and SHA-256
recorded by its own run manifest. This prevents different references or caller
contracts from being summarized as biological replicate disagreement.

`--batch` is mutually exclusive with `--sample` and with replicate-specific
single-sample arguments. Mixed explicit/Ribo-TISH rows are accepted but clearly
labelled; using one input mode consistently within a batch remains preferable.
The current batch contract requires the selected stage range to include
`psite`, because the TSV describes replicate-level BAM and offset inputs.

## 6. Starting at any stage

Valid stage names are:

```text
orfscan, psite, qc, dm, abundance, context, visualize
```

`codon` is an optional internal ORF P-site export stage. It is skipped by
default and runs only when explicitly selected or enabled with
`--export-orf-psites`.

Use a contiguous stage range:

```bash
--from-stage STAGE --to-stage STAGE
```

Run exactly one stage:

```bash
--only-stage STAGE
```

### 6.1 Inputs required for each starting point

| Starting stage | Required external inputs |
|---|---|
| `orfscan` | `--gtf`, `--fa`; if P-site calling is included, also `--bam` plus either `--length-offsets` or the Ribo-TISH pair |
| `psite` | `--bam` plus either `--length-offsets` or `--offsets` and `--ribotish-quality`; also an existing/default `--torf` if QC or DM is included |
| `qc` | `--torf`, `--psite-bam` |
| `dm` | `--torf`, `--psite-bam`, `--dm-background` |
| `abundance` | `--dm-results`, `--psite-bam` |
| `context` | `--dm-results`, `--torf` |
| `codon` | `--gene-context-results`, `--torf`, `--psite-bam`, `--dm-run-manifest` |
| `visualize` | `--gene-context-results` |

Start at QC and continue through visualization:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --from-stage qc \
  --to-stage visualize \
  --torf /path/to/miao-orf-results/01_orfscan/gencode_v48.torf.tsv \
  --psite-bam /path/to/miao-orf-results/02_psite/ES-Ribo-Rep1/ES-Ribo-Rep1.psite.bam
```

Run only the formal DM caller:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --only-stage dm \
  --torf /path/to/miao-orf-results/01_orfscan/gencode_v48.torf.tsv \
  --psite-bam /path/to/miao-orf-results/02_psite/ES-Ribo-Rep1/ES-Ribo-Rep1.psite.bam \
  --dm-background /path/to/miao-orf-results/03_metagene_qc/ES-Ribo-Rep1/ES-Ribo-Rep1.dm_background.tsv \
  --dm-mode accurate
```

Quantify abundance only from an existing DM table and its sample-matched
1-nt P-site BAM:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --only-stage abundance \
  --psite-bam /path/to/sample.psite.bam \
  --dm-results /path/to/sample.intorf_dm_results.tsv
```

Regenerate only the user-facing plots:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --only-stage visualize \
  --gene-context-results /path/to/sample.gene_cds_context.tsv
```

## 7. Resume behavior and provenance

### 7.1 Resume by default

A stage is skipped only when its primary outputs pass content validation and
its atomic `<first-output>.complete.json` record still matches the exact command,
stage-program hash, required-input identities and output identities. Re-running
the same command therefore continues at the first missing, malformed or stale
stage. Non-empty legacy outputs without a completion record are reported as
`STALE`, not silently trusted.

Inspect expected outputs:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --from-stage dm --to-stage visualize \
  --status
```

### 7.2 Force a rerun

```bash
--force
```

`--force` invalidates the old completion record and reruns all selected stages.
The top-level runner does not recursively delete old directories, but individual
stage programs may replace outputs with the same prefix. Confirm the output
prefix before using this option.

### 7.3 Run records

Each executed stage writes a `run.log`. The top-level runner also writes:

```text
OUT_ROOT/pipeline_runs/SAMPLE/TIMESTAMP.pipeline_run.json
```

The JSON record contains:

- start and finish times;
- sample and reference labels;
- `accurate` or `fast` DM mode;
- the exact argument vector for every stage;
- a shell-readable copy of every command;
- expected outputs;
- stage completion-record path and validation state;
- completed, skipped, or failed status.

## 8. Default output layout

```text
OUT_ROOT/
|-- 01_orfscan/REFERENCE.*
|-- 02_psite/SAMPLE/SAMPLE.psite.bam
|-- 03_metagene_qc/SAMPLE/SAMPLE.dm_background.tsv
|-- 03_metagene_qc/SAMPLE/SAMPLE.codon_frame_ternary.tsv
|-- 03_metagene_qc/SAMPLE/SAMPLE.codon_frame_ternary.{png,pdf}
|-- 04_intorf_dm/SAMPLE/accurate/SAMPLE.intorf_dm_results.tsv
|-- 04_intorf_dm/SAMPLE/accurate/SAMPLE.intorf_abundance.tsv
|-- 04_intorf_dm/SAMPLE/accurate/SAMPLE.gene_cds_context.tsv
|-- 04_intorf_dm/SAMPLE/accurate/SAMPLE.orf_psite_codons.tsv.gz       # optional
|-- 04_intorf_dm/SAMPLE/accurate/SAMPLE.orf_psite_ternary.tsv         # optional
|-- 04_intorf_dm/SAMPLE/accurate/SAMPLE.orf_psite_summary.tsv         # optional
|-- 05_visualization/SAMPLE/accurate/SAMPLE.*
`-- pipeline_runs/SAMPLE/TIMESTAMP.pipeline_run.json
```

For multiple samples, reuse the same `OUT_ROOT`. A completed reference scan is
then reused automatically, while each new `--sample` receives separate P-site,
QC, DM, abundance, context, visualization, and run-record outputs. Optional ORF
P-site tables appear only when explicitly enabled.

With `--dm-mode fast`, the DM and visualization directories use `fast`
instead of `accurate`. Every stage output prefix can also be overridden with a
corresponding `--*-out-prefix` option.

## 9. Stage 1: transcriptome ORF scan

Program:

```text
src/miao_orf/orf_scan_transcriptome.py
```

Top-level defaults:

- minimum ORF length: 6 aa;
- start codon: ATG only;
- primary assembly contigs only;
- chromosome-partitioned execution;
- worker count inherited from `--workers`.

Enable near-cognate CTG/TTG/GTG starts with:

```bash
--near-cognate
```

Primary outputs:

| File | Description |
|---|---|
| `*.torf.tsv` | Transcript-level ORFs, annotations, and gORF membership |
| `*.torf.faa` | tORF peptide FASTA |
| `*.gorf.tsv` | Genomic ORF summary |
| `*.gorf_members.tsv` | tORF-to-gORF membership |
| `*.gorf.faa` | gORF peptide FASTA |
| `*.gorf_validation.tsv` | gORF validation summary |

The DM caller focuses on `intORF_altframe` rows in `*.torf.tsv`. Altframe
intORFs that share an annotated stop codon elsewhere in the same gene are
marked as stop-confounded. They remain auditable in scan output but are
excluded from formal DM analysis by default.

## 10. Stage 2: P-site BAM

Program:

```text
src/miao_orf/psite-caller.py
```

Defaults:

- accept exactly one formal offset-input mode: traditional explicit
  `--length-offsets`, or automatic Ribo-TISH `--offsets` plus
  `--ribotish-quality`;
- safely parse the sample-matched Ribo-TISH `mapped.para.py` and
  `mapped_qual.txt` without executing either file when automatic mode is used;
- use the main offset group and discard reads with a mismatch at the aligned
  5' nucleotide (the Ribo-TISH `m0` group is not used);
- rotate Ribo-TISH QC frame counts from its fixed 12-nt anchor to each
  length-specific offset, then retain lengths with corrected frame-0
  proportion strictly greater than `--min-frame0-proportion 0.6666666667`;
- among passing lengths, retain the contiguous block with the largest total
  matched-read support (`--length-selection-policy dominant_contiguous`);
- MAPQ >= 20;
- NH=1 when the tag is present;
- reads with an unmatched or unusable 5' end are discarded;
- output is a coordinate-sorted and indexed 1-nt BAM;
- BED and bedGraph output are disabled by the top-level runner to save space.

Primary outputs:

```text
SAMPLE.psite.bam
SAMPLE.psite.bam.bai
SAMPLE.offset_selection.tsv
```

`SAMPLE.offset_selection.tsv` records the selected input mode. In automatic
mode it reports every observed length, its main offset, raw and corrected
P0/P1/P2 counts, corrected P0 proportion, length-selection policy, selection
status and rejection reason. In traditional mode it records every explicit
pair. Use `--keep-lengths 28,29,30` only as an additional automatic-mode
whitelist. The older `--offsets para.py --keep-lengths ...` command remains
supported solely to reproduce frozen analyses.

The offset decision can be inspected without reading the BAM:

```bash
python3 src/miao_orf/ribotish_offsets.py \
  --offsets SAMPLE.mapped.para.py \
  --ribotish-quality SAMPLE.mapped_qual.txt \
  --out SAMPLE.offset_selection.tsv
```

The DM caller consumes this P-site BAM and does not infer a new P-site from the
original full-length read.

Enable optional text tracks with:

```bash
--emit-bed
--emit-bedgraph
```

## 11. Stage 3: metagene QC and DM background

Program:

```text
src/miao_orf/ribo_metagene_qc.py
```

This stage:

- calculates MANE CDS start/stop metagene profiles;
- reports P0/P1/P2 distributions;
- writes a codon-level phase-composition QC on a fixed 5% ternary grid, with
  bubble area equal to the percentage of A0-input codons at each grid point;
- generates a scaled CDS profile;
- estimates the translated CDS template `pi_trans`;
- estimates the Dirichlet-multinomial concentration `A0`;
- validates `A0 >= --min-a0`.

Current defaults:

```text
pi_method = codon_equal
min_A0 = 1.0
```

Primary outputs:

| File | Purpose |
|---|---|
| `*.metagene.pdf` | Start, stop, and frame QC plots |
| `*.codon_frame_ternary.png`, `*.codon_frame_ternary.pdf` | Normalized codon-level P0/P1/P2 composition QC |
| `*.codon_frame_ternary.tsv` | Plot-ready 5% grid counts and linear sample percentages |
| `*.qc_summary.tsv` | Run-level QC summary |
| `*.template.tsv` | Template estimation details |
| `*.dm_background.tsv` | Validated background for the DM caller |
| `*.cds_metrics.tsv` | CDS-level QC measurements |

A very small A0 indicates excessive phase-pattern dispersion, insufficient
coverage, or a sample that is poorly described by the current model. Do not
bypass the A0 check only to force a DM result.

The ternary QC complements the conventional metagene profile: metagene plots
show the average periodic signal, whereas the ternary plot preserves
codon-level phase heterogeneity. Percentages are encoded once, by bubble area
on a linear scale; the figure intentionally does not use a logarithmic color
scale. Its header reports the exact codon-equal mean P0/P1/P2 composition; these
means are calculated before 5% grid projection and are also written to
`qc_summary.tsv`. A white star with a dark outline marks that same mean directly
inside the triangle. The fixed grid keeps samples directly comparable without
adding read-depth-stratified panels to the default report.

## 12. Stage 4: intORF DM caller

Program:

```text
src/miao_orf/ribo_intorf_dm_caller.py
```

### 12.1 Model

For each active intORF codon, the caller records P-site counts in three phases:

```text
H0: lambda = 0
H1: lambda > 0

pi(lambda) = (1 - lambda) * pi_host + lambda * pi_trans
counts_j ~ DM(n_j, A0 * pi(lambda))
```

Here:

- `pi_host` is the annotated-CDS template in intORF-relative coordinates;
- `pi_trans` is the translated-ORF template from QC;
- `lambda_hat` is the fitted strength of the additional intORF-like component;
- `A0` is taken from the validated QC background.

`lambda_hat` is a model-based fraction of the overlap-region P-site component.
It is not an absolute protein abundance or translation-efficiency estimate.
The caller also reports `lambda_profile_ci95_low` and
`lambda_profile_ci95_high`, a likelihood-ratio profile interval bounded to
`[0, 1]`. This is a frequentist uncertainty interval, not a posterior credible
interval.

For candidates with a fitted model, the caller also reports fractional
model-expected component counts in the analyzed core:

```text
model_expected_intorf_core_reads = core_reads * lambda_hat
model_expected_host_cds_core_reads = core_reads * (1 - lambda_hat)
```

The corresponding `phase0`, `phase1`, and `phase2` columns partition each
component total using `pi_trans_*` or `pi_host_*`, respectively. Phase indices
are in intORF-relative coordinates. These values are fitted expectations, not
integer read-level assignments: an individual read in an overlapping interval
cannot generally be attributed uniquely to the host CDS or intORF. Rows without
an eligible fitted model leave these fields empty.

### 12.2 Formal inference modes

The default is:

```bash
--dm-mode accurate
```

`accurate` uses the formal adaptive-importance configuration and is intended
for final results.

```bash
--dm-mode fast
```

`fast` uses the validated faster execution policy and is appropriate for
screening, parameter checks, and large preliminary runs. Both modes use the
same biological model and final gates; they differ in extreme-tail p-value
computation cost.

```bash
--dm-mode both
```

`both` runs the shared upstream stages once, then writes independent
`accurate` and `fast` DM result directories. It is intended for validation or
side-by-side mode comparison; `accurate` remains the normal single-mode
default.

### 12.3 Minimum evidence for inference

Default eligibility requirements are:

- at least 5 core codons;
- at least 3 active core codons;
- at least 15 core P-site reads.

Candidates below these requirements are labelled `insufficient_data`. This is
an indeterminate low-information state, not proof of absent translation.

### 12.4 BH correction and final gates

The DM `p_final` is the only p-value used for BH correction. The formal caller
does not combine DM p-values with uniformity, stop, or noise p-values.

After BH correction, a primary credible call must pass:

1. `q_BH < 0.05`;
2. `lambda_hat > 0.05` and the expected effect direction;
3. gain/drop and H-to-T mixture-geometry consistency;
4. at least 5 active core codons;
5. active-core fraction >= 0.15;
6. target-residual-supported fraction strictly greater than 1/3.

These gates do not replace or recompute the DM p-value. They classify whether
a statistically significant pattern is broad and mechanistically compatible
with the expected CDS/intORF mixture.

### 12.5 Main classifications

| `classification` | Interpretation |
|---|---|
| `credible_extra_ORF_like_signal` | BH and all primary gates pass |
| `host_only_supported` | DM does not reach BH significance or does not enter formal tail calculation |
| `atypical_pattern` | DM-significant, but direction or mixture geometry is inconsistent |
| `localized_core_signal` | Signal is supported by too few codons |
| `low_target_frame_residual_breadth` | Coverage is present, but too few active codons specifically support the target frame |
| `insufficient_data` | Minimum evidence for inference is not available |
| `ultra_short_exploratory` | Below the primary length threshold but credible in the separate exploratory family |

### 12.6 Important result columns

| Column | Meaning |
|---|---|
| `torf_id`, `gorf_id` | Transcript-level and genomic ORF identifiers |
| `classification` | Final result class |
| `primary_credible_call` | Formal primary credible-call indicator |
| `p_final`, `q_BH`, `q_BY` | DM p-value and multiple-testing corrections |
| `lambda_hat` | Codon-level DM estimate of the additional component |
| `lambda_profile_ci95_low`, `lambda_profile_ci95_high` | Bounded 95% likelihood-ratio profile interval for `lambda_hat` |
| `orf_region_reads` | Observed P-sites across the complete candidate ORF interval |
| `core_reads`, `analyzed_core_reads` | P-sites in the analyzed core; the second name is the explicit public-schema alias |
| `observed_phase{0,1,2}_core_reads`, `observed_phase{0,1,2}_core_fraction` | Observed core composition in intORF-relative coordinates |
| `model_expected_host_cds_fraction`, `model_expected_intorf_fraction` | Fitted component fractions, equal to `1-lambda_hat` and `lambda_hat` |
| `model_expected_host_cds_core_reads`, `model_expected_intorf_core_reads` | Fractional model-expected host-CDS and intORF-like P-site counts in the analyzed core |
| `model_expected_host_cds_phase{0,1,2}_core_reads`, `model_expected_intorf_phase{0,1,2}_core_reads` | Model-expected component counts partitioned by intORF-relative phase |
| `dm_evidence_score`, `dm_fdr_score` | Finite `-log10(p_final)` and `-log10(q_BH)` scores |
| `qc_status`, `qc_flags`, `filter_reason` | Compact fit/QC outcome and explicit exclusion reason |
| `background_template_sensitivity_status` | Whether the fitted point is geometrically consistent with the host-to-translated template segment |
| `n_active_core_codons` | Core codons with at least one P-site |
| `active_core_codon_frac` | Fraction of core codons that are active |
| `target_vs_unused_z05_frac_active` | Fraction of active codons supporting the target frame over the unused frame |
| `distance_to_mixture_segment` | Distance from observed phase composition to the H-to-T segment |
| `mixture_geometry_consistent` | Gain/drop consistency indicator |

Primary output files:

```text
*.intorf_dm_results.tsv
*.summary.tsv
*.template.tsv
*.prefilter_summary.tsv
*.run_manifest.json
```

The JSON run manifest records the software/schema/engine versions, all actual
parameters, runtime, input file identities, caller source hash, output schema
and read-count semantics. The certified benchmark defaults and artifact hashes
are frozen in `config/releases/benchmark_520_dual_v1.0.0.json`.

## 13. Stage 5: model-allocated intORF abundance

Program:

```text
src/miao_orf/quantify_intorf_abundance.py
```

This post-processing stage reads the DM result table and the same sample 1-nt
P-site BAM used by the caller. It does not refit the DM model and does not
change calls, p-values, FDR results, gates, or classifications.

For each row with a fitted intORF component it calculates:

```text
effective_core_nt = 3 * n_core_codons
usable_library_psites = mapped alignments in the 1-nt P-site BAM
intorf_pFPKM = 1e9 * model_expected_intorf_core_reads
                / (effective_core_nt * usable_library_psites)
```

The numerator and effective length therefore describe the same DM-analyzed
core. The library denominator automatically inherits the read-length subset,
P-site offsets, MAPQ, uniqueness, and 5'-match policy already applied when the
1-nt P-site BAM was built. This permits sample-to-sample comparison under the
same interpretation while keeping each sample's own usable library size.
Samples should still use compatible library preparation, alignment, and
P-site-selection policies. Library-size normalization does not erase a
systematic difference caused by selecting biologically different footprint
length subsets; the per-sample manifest therefore states that the subset is
inherited from its supplied P-site BAM.

The likelihood-profile limits for `lambda_hat` are propagated to
`intorf_pFPKM_ci95_low` and `intorf_pFPKM_ci95_high`. The table also includes
observed-core pFPKM, host-component pFPKM, intORF P-site RPM, and the fitted
intORF/host ratio. These are model-allocated translation-signal abundance
measures; they are not protein concentration and are not translation
efficiency because no RNA abundance denominator is used.

Primary outputs:

```text
*.intorf_abundance.tsv
*.intorf_abundance_summary.tsv
*.intorf_abundance_manifest.json
```

Rows without an eligible fitted DM model are retained with
`abundance_status=not_model_quantifiable` and blank allocated abundance rather
than being silently discarded. The manifest records the exact formula, input
identities, mapped P-site denominator, inherited length-subset semantics, and
output schema.

The script can also be run independently:

```bash
python3 src/miao_orf/quantify_intorf_abundance.py \
  --sample ES-Ribo-Rep1 \
  --dm-results /path/to/sample.intorf_dm_results.tsv \
  --psite-bam /path/to/sample.psite.bam \
  --out-prefix /path/to/sample
```

## 14. Stage 6: gene-level annotated-CDS N-terminal context

Program:

```text
src/miao_orf/annotate_gene_cds_context.py
```

This postprocessing stage addresses cases in which an ORF is correctly typed
as `intORF_altframe` on its own transcript but reuses the N terminus of an
annotated CDS from another transcript of the same gene. It does not change the
transcript-level ORF type, DM fit, p/q values, lambda, gates, or primary call.

A match requires all of the following: exact versioned `gene_id`, same strand,
the candidate start at an annotated-CDS codon boundary, and an identical
translation-oriented genomic coordinate path from the candidate start. A
simple genomic interval or BED overlap is not sufficient. The default reported
match threshold is five complete coordinate-identical N-terminal codons.

Primary classes are:

- `pure_intorf_no_annotated_cds_nterm_reuse`;
- `annotated_cds_nterm_reuse_with_splice_derived_cterm`;
- `annotated_cds_full_peptide_alternative_coordinate_path`;
- `annotated_cds_derived_full_coordinate_path`.

`gene_level_pure_intorf_eligible=1` provides the unified gene-level subset for
downstream novelty summaries. The original transcript-level classification is
retained in the same row for auditability. Coordinate-prefix nt/codon counts,
peptide-prefix length, best and tied annotated transcripts, the break reason,
summary counts, and a SHA-256 input manifest are written to:

```text
*.gene_cds_context.tsv
*.gene_cds_context_only.tsv
*.gene_cds_context_summary.tsv
*.gene_cds_context_manifest.json
```

Existing HPC results can be annotated without rerunning ORFscan, P-site
calling, QC, or DM:

```bash
python3 src/miao_orf/annotate_gene_cds_context.py \
  --torf /path/to/reference.torf.tsv \
  --input-tsv /path/to/sample.intorf_dm_results.tsv \
  --out-prefix /path/to/sample
```

The replicate batch workflow applies the same annotator to each sample through
the normal stage plan and to the final consensus table after integration. Raw
DM and integration tables remain complete. Final single-sample and replicate
count tables, lambda/pFPKM matrices, and figures automatically retain only
`gene_level_pure_intorf_eligible=1`.

## 15. Optional internal tool: selected-ORF codon-level P-site export

Program:

```text
src/miao_orf/export_orf_psites.py
```

The DM result table contains candidate-level totals such as
`observed_phase0_core_reads`, but it intentionally does not retain every
per-codon phase vector. Therefore those original result columns cannot be used
to reconstruct a ternary distribution for an individual ORF.

This postprocessing stage reconstructs those data from the matching 1-nt
P-site BAM and tORF genomic blocks. It reads the exact start/stop codon
exclusions and optional read-length contract from the matching DM run manifest.
Before writing outputs, every candidate is required to reproduce the saved core
read total, active/core codon counts, and P0/P1/P2 totals exactly. A mismatch
stops the stage, which protects against accidental use of another sample's BAM,
tORF table, or DM contract.

The normal pipeline does not run this tool. `--export-orf-psites` enables it in
a larger stage range, while `--only-stage codon` runs it alone. Without IDs the
top-level runner selects rows satisfying both `primary_credible_call=1` and
`gene_level_pure_intorf_eligible=1`. The standalone module defaults to primary
credible calls. Explicit IDs override credibility and gene-level filters so a
specific failed or diagnostic candidate can be inspected.

Primary outputs are:

| File | Purpose |
|---|---|
| `*.orf_psite_codons.tsv.gz` | One row per plotted active core codon, with genomic positions, raw P0/P1/P2 counts, fractions, and percentages |
| `*.orf_psite_ternary.tsv` | One row per candidate/exact reduced P0:P1:P2 ratio, with its unrounded coordinate, codon count, and percentage within that candidate |
| `*.orf_psite_summary.tsv` | One row per selected ORF, including exact codon-equal mean P0/P1/P2 percentages and vertex percentages |
| `*.orf_psite_manifest.json` | Input identities, inherited DM parameters, selection, semantics, and row counts |

The ternary percentages are codon-equal: each plotted codon contributes one
phase-composition point regardless of read depth. Coordinates are exact: codons
are grouped only when their reduced integer P0:P1:P2 ratios match, with no
binning or coordinate rounding. The raw per-codon read counts
remain available in the compressed long table. Because the intORF overlaps a
host CDS, these observed counts describe the overlap mixture and do not assign
individual reads uniquely to the intORF or host-CDS component.

Run the standalone module on completed results without rerunning ORFscan, P-site
calling, QC, or DM:

```bash
python3 src/miao_orf/export_orf_psites.py \
  --input-results /path/to/sample.gene_cds_context.tsv \
  --torf /path/to/reference.torf.tsv \
  --psite-bam /path/to/sample.psite.bam \
  --dm-run-manifest /path/to/sample.run_manifest.json \
  --out-prefix /path/to/sample \
  --sample SAMPLE \
  --gene-level-pure-intorf-only
```

Inspect one ORF (the identifier can instead be a `candidate_key` or `torf_id`):

```bash
python3 src/miao_orf/export_orf_psites.py \
  --input-results /path/to/sample.gene_cds_context.tsv \
  --torf /path/to/reference.torf.tsv \
  --psite-bam /path/to/sample.psite.bam \
  --dm-run-manifest /path/to/sample.run_manifest.json \
  --out-prefix /path/to/inspection \
  --orf-id GORF_ID
```

For batch inspection, repeat `--orf-id`, pass comma-separated IDs, or use
`--orf-list ids.txt`. The list may contain one ID per line or TSV columns named
`orf_id`, `candidate_key`, `gorf_id`, and/or `torf_id`. Use `--selection all`
only when an unfiltered bulk export is intentionally required.

## 16. Stage 7: result visualization

Program:

```text
src/miao_orf/visualize_intorf_dm_results.py
```

The top-level runner supplies the context-annotated caller output and enables
gene-level pure-intORF filtering. The visualization module does not refit the
model, recalculate p-values, or alter classifications. Its manifest records the
raw row count, the context-excluded count, and the final plotted row count.

Default end-user plots:

- `model_geometry_ternary`: all eligible candidates in CDS coordinates;
- `model_geometry_ternary_by_frame`: separate F1 and F2 ternary views;
- `gate_waterfall`: sequential candidate counts after each gate;
- `gate_combinations`: common post-BH gate combinations;
- `breadth_gate_plane`: active breadth versus target-residual breadth;
- `effect_significance`: `lambda_hat` versus BH significance;
- `length_stratification`: results by peptide-length bin;
- `depth_effects`: relationships between read depth and model results.

Default colors:

- blue: credible;
- orange: atypical geometry;
- pink: localized core signal;
- yellow: low residual breadth;
- gray: all eligible background candidates.

Developer diagnostics are disabled by default. Enable them explicitly with:

```bash
--include-diagnostics
```

Default software output does not include highlighted DEDD2 panels, PRICE
comparisons, ROC/PR curves, or article-specific benchmark figures.

## 17. Frequently used options

### 17.1 Compute resources

```bash
--workers 8
```

ORF scan, P-site calling, QC, and DM use this default worker count. The DM
multiprocessing chunksize can be set separately:

```bash
--dm-mp-chunksize 1
```

### 17.2 ORF length thresholds

```bash
--min-aa 6
--min-intorf-aa 6
--primary-min-intorf-aa 10
```

- `min-aa`: minimum ORF length written by transcriptome scanning;
- `min-intorf-aa`: minimum length eligible for DM analysis;
- `primary-min-intorf-aa`: minimum length in the primary FDR family.

Candidates below the primary threshold should not be combined with primary
credible calls in the main result count.

### 17.3 Read lengths

For traditional input, specify the complete mapping directly:

```bash
--length-offsets 28:12,29:12,30:12
```

Do not combine this option with Ribo-TISH inputs. For automatic input, read
lengths are selected from the Ribo-TISH quality/offset pair. The threshold can
be changed with:

```bash
--min-frame0-proportion 0.6666666667
```

The comparison is strict (`P0 proportion > threshold`). The default
`--length-selection-policy dominant_contiguous` then retains the contiguous
passing-length block with the greatest total matched-read count. To retain every
passing length, including isolated lengths, use:

```bash
--length-selection-policy all_passing
```

To impose an additional comma-separated whitelist:

```bash
--keep-lengths 28,29,30
```

The older `--offsets para.py --keep-lengths ...` form bypasses periodicity
screening and exists only to reproduce frozen legacy analyses.

### 17.4 Advanced pass-through arguments

Options not directly exposed by the top-level runner can be appended to a
stage:

```bash
--dm-extra='--importance-max-relative-se 0.20'
--visualize-extra='--max-background-points 50000'
```

Use the `--option='...'` form when the value begins with `--`, so it is not
parsed as another top-level option.

## 18. Multiple samples

- Run reference-level ORF scanning once per GTF/FASTA combination.
- Generate a separate P-site BAM for each Ribo-seq sample.
- Perform metagene QC and background estimation separately for each sample.
- Run the DM caller separately for each sample.
- Do not reuse one sample's `dm_background.tsv` for another sample unless that
  choice has been explicitly validated.
- Perform multi-sample integration downstream of the single-sample caller.

## 19. Troubleshooting

### 19.1 Missing MD tag

The mapped BAM lacks MD tags. Run `samtools calmd` with the matching reference
FASTA, then sort and re-index the resulting BAM.

### 19.2 BAM index is older than the BAM

The BAM modification time is newer than its index. Rebuild it:

```bash
samtools index -f sample.bam
```

If the warning remains, check WSL/Windows filesystem times and whether the
program is reading `sample.bam.bai` or `sample.bai`.

### 19.3 No allowed read lengths

No main-offset length has corrected frame-0 proportion strictly above the
configured threshold, or an optional `--keep-lengths` whitelist removed all
passing lengths. Inspect `*.offset_selection.tsv`; also confirm that
`mapped_qual.txt` and `mapped.para.py` came from the same BAM.

In traditional mode, confirm that every item uses `LENGTH:OFFSET`, that lengths
are unique, and that each offset satisfies `0 <= offset < length`.

### 19.4 Per-contig P-site parts remain

The formal top-level command enables merging. Interrupted runs may leave the
parts directory for diagnosis. Treat the final sorted and indexed
`*.psite.bam` as the primary output and remove parts only after confirming that
the final BAM is complete.

### 19.5 QC does not produce a usable DM background

Inspect `qc_summary.tsv`, the metagene PDF, frame composition, and A0. Do not
lower `min-a0` only to force the workflow to continue; a very low A0 indicates
that the current DM model may not be reliable for that sample.

### 19.6 Accurate DM mode is slow

`accurate` mode performs adaptive importance sampling for extreme tails. Its
runtime depends on the number of eligible candidates, candidates near the BH
boundary, the requested importance samples, and worker count. Use
`--dm-mode fast` for preliminary inspection, then retain `accurate` and its
run record for final results.

### 19.7 A candidate is absent from the result

Check, in order:

1. whether the ORF was written to `torf.tsv`;
2. whether it is classified as `intORF_altframe`;
3. whether it is annotated-stop-confounded and therefore excluded by default;
4. whether it passes the DM length threshold;
5. whether it enters the coverage-eligible set;
6. whether gORF deduplication selected another member as its representative.

### 19.8 Does `insufficient_data` mean not translated?

No. It means the current Ribo-seq sample does not provide enough evidence for
a reliable decision. Low coverage may reflect low expression, condition-specific
inactivity, or insufficient sequencing depth.

### 19.9 Why are there fewer credible calls than BH-significant candidates?

BH tests DM significance only. A credible call must additionally pass effect
direction, mixture geometry, active-codon breadth, and target-residual breadth
gates. Significant candidates that fail a gate remain in the output with an
explicit classification.

## 20. Recommended reporting

At minimum, retain and report:

- GENCODE and genome-build versions;
- the matching Ribo-TISH quality and offset files;
- `offset_selection.tsv`, retained read lengths, corrected frame-0 threshold,
  MAPQ threshold, and unique-mapping policy;
- metagene QC, `pi_trans`, and A0;
- `accurate` or `fast` mode;
- ORF length thresholds;
- BH and gate thresholds;
- the number of primary credible calls;
- counts of `insufficient_data`, `atypical_pattern`, `localized_core_signal`,
  and `low_target_frame_residual_breadth`;
- the corresponding `pipeline_run.json`.

Do not retain only a filtered credible-call table. The full
`intorf_dm_results.tsv` is necessary for reviewing low-coverage candidates,
borderline results, and gate behavior.

## 21. Help and issue reports

Top-level help:

```bash
python3 miao_orf.py --help
```

Each stage remains independently runnable:

```bash
python3 src/miao_orf/orf_scan_transcriptome.py --help
python3 src/miao_orf/psite-caller.py --help
python3 src/miao_orf/ribo_metagene_qc.py --help
python3 src/miao_orf/ribo_intorf_dm_caller.py --help
python3 src/miao_orf/quantify_intorf_abundance.py --help
python3 src/miao_orf/annotate_gene_cds_context.py --help
python3 src/miao_orf/visualize_intorf_dm_results.py --help
```

When reporting a problem, include:

- the exact command or `pipeline_run.json`;
- the relevant stage `run.log`;
- Python, dependency, and samtools versions;
- at least 20 lines around the error;
- confirmation that required inputs and indexes exist and are non-empty.
