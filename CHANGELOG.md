# Changelog

All notable public changes to MIAO are documented here. The project
uses semantic versioning for public releases; earlier 1.2/1.3 labels were
unpublished development snapshots.

## [Unreleased]

No unreleased changes are currently recorded.

## [1.0.0] - 2026-09-04

### Added

- unified eight-stage `miao-orf` pipeline, including one optional codon-export stage;
- explicit and automatic Ribo-TISH length/offset input modes;
- accurate and fast intORF-DM inference modes;
- observed and fitted component read-count reporting;
- replicate batch execution, evidence-preserving consensus and replicate plots;
- atomic run manifests and versioned public output schema;
- installable package metadata and reproducible Conda environment;
- content-validated stage completion records and safe resume semantics;
- strict DM replicate compatibility validation;
- standalone model-allocated intORF pFPKM quantification from DM results and
  the sample 1-nt P-site BAM;
- abundance provenance, uncertainty propagation and replicate-level pFPKM
  correlation matrices and absolute-abundance heatmaps;
- dynamic lower-triangle pFPKM and lambda concordance layouts that adapt to the
  actual replicate count;
- an auditable postprocessing stage for same-gene annotated-CDS N-terminal
  reuse, using translation-oriented genomic coordinate paths and codon phase;
- context-aware final plots and count tables while retaining complete raw DM,
  abundance and replicate-integration tables for audit.

### Fixed

- accurate fixed-size importance sampling now records its actual Monte Carlo
  draw count in per-candidate output and aggregate provenance summaries.

### Release identity

- the method is named **MIAO**, the repository is **MIAO-ORF**, the command and
  distribution are `miao-orf`, and the Python package is `miao_orf`.
