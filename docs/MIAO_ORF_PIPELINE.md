# MIAO top-level pipeline

For the complete English user manual, see [`MANUAL.md`](MANUAL.md).

`miao_orf.py` is a thin runner for seven default stages plus one optional
internal export stage. It does not copy or reimplement their scientific logic.

1. `orfscan`: `orf_scan_transcriptome.py`
2. `psite`: `psite-caller.py`
3. `qc`: `ribo_metagene_qc.py`
4. `dm`: `ribo_intorf_dm_caller.py`
5. `abundance`: `quantify_intorf_abundance.py`
6. `context`: `annotate_gene_cds_context.py`
7. `visualize`: `visualize_intorf_dm_results.py`

Optional internal stage: `codon` (`export_orf_psites.py`). It is not
part of a default run and must be enabled explicitly.

Run the top-level program in the same Python/conda environment that contains
the dependencies of the stage programs. The current project uses the `bio`
environment.

```bash
cd /path/to/MIAO
conda activate bio
python3 miao_orf.py --help
```

## Full run

`--out-root` must name the shared analysis root. Do not add the sample name to
it; `--sample` creates the sample-level paths. This keeps `01_orfscan` reusable
across all samples.

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

By default, P-site calling uses only the Ribo-TISH 5'-matched group and keeps
read lengths whose length-specific offset produces corrected frame-0
proportion strictly greater than 2/3. It then retains the contiguous passing-
length block with the largest total matched-read support. The selected lengths,
offsets, policy and rejection reasons are recorded in `*.offset_selection.tsv`.
The Ribo-TISH `mapped_qual.pdf` may be reviewed manually but is never used as a
pipeline input; all decisions use exact counts from `mapped_qual.txt`.

Alternatively, bypass Ribo-TISH parsing and provide a traditional explicit
mapping with `--length-offsets 28:12,29:12,30:12`. The explicit and automatic
input modes are mutually exclusive and both are recorded in the offset audit.

The default formal DM mode is `accurate`. Use `--dm-mode fast` for the faster
validated alternative, or `--dm-mode both` to run the shared upstream stages
once and then produce independent `accurate` and `fast` DM outputs.

## Replicate batch

Prepare a tab-separated file with columns `replicate_id`, `bam`, `input_mode`,
`length_offsets`, `ribotish_para`, and `ribotish_quality`, then run:

```bash
python3 miao_orf.py \
  --batch replicates.tsv \
  --out-root /path/to/miao-orf-results/pipeline \
  --gtf /path/to/annotation.gtf \
  --fa /path/to/genome.fa \
  --workers 8
```

The runner executes replicates sequentially, shares ORFscan, never pools BAMs,
continues after a replicate failure, and writes `batch_runs/*.batch_summary.tsv`
plus a JSON batch record. Relative row paths are resolved from the TSV folder.
If DM is selected, it then invokes `integrate_replicates.py` and writes long,
consensus, summary and manifest outputs under `06_replicate_integration/MODE/`.
Consensus is based on `gorf_id + overlap_type` replicate support; reads and
p-values are never pooled or combined. Successful integration automatically
applies the same gene-level annotated-CDS N-terminal context annotator to the
consensus table. The annotator requires exact versioned gene, strand, genomic
translation path and codon phase; interval overlap alone is not accepted. It
preserves transcript-level ORF types and all DM statistics. Integration then
invokes `visualize_replicate_integration.py` on the context-annotated consensus.
All final replicate call counts, support summaries, lambda/pFPKM matrices and
figures retain only `gene_level_pure_intorf_eligible=1`; the unfiltered
integration tables remain unchanged for audit. The visualizer writes PNG/PDF call-count,
support-combination, pairwise-lambda-concordance, absolute-lambda heatmap and
primary-reproducibility figures plus a plot summary and manifest in the same
directory. Both lambda concordance variants use dynamic lower-triangle scatter
matrices with per-replicate distributions on the diagonal and therefore adapt
to the actual number of input samples. When per-replicate abundance tables are
available, it also writes a matching lower-triangle pFPKM scatter matrix, a
Spearman correlation matrix, an absolute-abundance pFPKM heatmap, and auditable
lambda/pFPKM matrices and correlation tables. The heatmap intentionally omits
within-candidate row-z scores because replicate batches are treated as
within-group measurements. pFPKM concordance uses
`log2(1 + pFPKM)` and candidates primary credible in any replicate; missing
values are not converted to zero. Unavailable replicates are shown separately
from observed negative calls. Use `--formats` and `--dpi` to control figure
output.

Before integration, every result is checked against its DM `run_manifest.json`.
Software/schema/engine versions, caller source, semantic parameters, tORF input
identities and output schema must agree, and each result must still match the
size/hash recorded by its own manifest.

For non-interactive launchers, use
`conda run --no-capture-output -n bio python miao_orf.py ...`; invoking only
the environment's Python executable is insufficient if `samtools` is not also
on `PATH`.

## Start at any stage

The selected stages must be contiguous. `--only-stage` is shorthand for using
the same value for `--from-stage` and `--to-stage`.

Start at metagene QC and continue through visualization:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --from-stage qc \
  --torf /path/to/miao-orf-results/01_orfscan/gencode_v48.torf.tsv \
  --psite-bam /path/to/miao-orf-results/02_psite/ES-Ribo-Rep1/ES-Ribo-Rep1.psite.bam
```

Run only the formal DM caller from existing inputs:

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

Run only visualization from an existing context-annotated result table:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --only-stage visualize \
  --gene-context-results /path/to/sample.gene_cds_context.tsv
```

Run only model-allocated abundance quantification from existing outputs:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --only-stage abundance \
  --psite-bam /path/to/sample.psite.bam \
  --dm-results /path/to/sample.intorf_dm_results.tsv
```

The same component is independently runnable as
`src/miao_orf/quantify_intorf_abundance.py`. It reports
`intorf_pFPKM = 1e9 * model_expected_intorf_core_reads /
(3 * n_core_codons * mapped_1nt_psite_alignments)`. The length subset and all
mapping filters are inherited from the supplied 1-nt P-site BAM. This is a
read-only normalization of fitted DM component counts; it does not refit the
model or change calls, p-values, FDR, thresholds, or classifications.

Run only the gene-level annotated-CDS N-terminal context stage on an existing
DM result and the matching full tORF table:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --only-stage context \
  --torf /path/to/reference.torf.tsv \
  --dm-results /path/to/sample.intorf_dm_results.tsv
```

This creates a full annotated result table, compact context table, summary and
manifest beside the DM output. `gene_level_pure_intorf_eligible` is the unified
gene-level novelty flag; the original `intORF_altframe` and DM call columns are
not overwritten.

Run the optional ORF P-site exporter from existing final calls. Coordinates are
exact normalized P0/P1/P2 compositions; this tool does not bin finite codons
onto the sample-level 5% QC grid:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --only-stage codon \
  --torf /path/to/reference.torf.tsv \
  --psite-bam /path/to/sample.psite.bam \
  --gene-context-results /path/to/sample.gene_cds_context.tsv \
  --dm-run-manifest /path/to/sample.run_manifest.json
```

`--only-stage codon` is an explicit opt-in. In a larger stage range, use
`--export-orf-psites`. With no IDs the runner exports final primary credible,
gene-level pure intORFs. Add `--orf-psite-id ID` once or repeatedly (comma-
separated values are accepted), or `--orf-psite-list FILE`, to inspect selected
ORFs regardless of their credibility flag. IDs may be `candidate_key`,
`gorf_id`, or `torf_id`. The exporter verifies recomputed core and P0/P1/P2
totals against the DM table before writing any result.

## Resume, inspect and override

Validated outputs are skipped only when their atomic `*.complete.json` record
still matches the stage command, source program, required inputs and outputs.
Rerunning the same command resumes at the first stale or incomplete stage.
`--force` invalidates completion before explicitly rerunning selected stages.
The runner does not recursively delete previous outputs itself.

Inspect output status:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --from-stage dm --to-stage visualize \
  --status
```

Validate external inputs and print exact commands without running them:

```bash
python3 miao_orf.py \
  --sample ES-Ribo-Rep1 \
  --out-root /path/to/miao-orf-results/pipeline \
  --only-stage visualize \
  --dm-results /path/to/sample.intorf_dm_results.tsv \
  --dry-run
```

Advanced stage options can be appended without changing the runner. Use the
equals form when the value begins with `--`:

```bash
--dm-extra='--importance-max-relative-se 0.20'
```

## Output layout

Unless an output-prefix override is supplied, the runner uses:

```text
OUT_ROOT/
├── 01_orfscan/REFERENCE.*
├── 02_psite/SAMPLE/SAMPLE.psite.bam
├── 03_metagene_qc/SAMPLE/SAMPLE.dm_background.tsv
├── 03_metagene_qc/SAMPLE/SAMPLE.codon_frame_ternary.tsv
├── 03_metagene_qc/SAMPLE/SAMPLE.codon_frame_ternary.{png,pdf}
├── 04_intorf_dm/SAMPLE/{accurate|fast}/SAMPLE.intorf_dm_results.tsv
├── 04_intorf_dm/SAMPLE/{accurate|fast}/SAMPLE.intorf_abundance.tsv
├── 04_intorf_dm/SAMPLE/{accurate|fast}/SAMPLE.gene_cds_context.tsv
├── 04_intorf_dm/SAMPLE/{accurate|fast}/SAMPLE.orf_psite_codons.tsv.gz       # optional
├── 04_intorf_dm/SAMPLE/{accurate|fast}/SAMPLE.orf_psite_ternary.tsv         # optional
├── 04_intorf_dm/SAMPLE/{accurate|fast}/SAMPLE.orf_psite_summary.tsv         # optional
├── 05_visualization/SAMPLE/{accurate|fast}/SAMPLE.*
└── pipeline_runs/SAMPLE/TIMESTAMP.pipeline_run.json
```

Each executed stage has a `run.log`. The JSON run record contains the exact
argument vector, a shell-readable command, expected outputs and completion
status for every selected stage, including skipped stages.
