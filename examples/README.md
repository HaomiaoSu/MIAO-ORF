# Examples

This directory contains sample-specific command and configuration examples.
They are templates rather than required MIAO source files.

- `commands/` contains the earlier stage-by-stage ES-Ribo-Rep1 commands.
- `configs/ES-Ribo-Rep1.mapped.para.py` illustrates the Ribo-TISH parameter
  format. Formal P-site calling should use the sample's matching
  `mapped.para.py` and `mapped_qual.txt` files together.
- If lengths and offsets have already been determined, traditional mode accepts
  them directly as `--length-offsets 28:12,29:12,30:12` and requires neither
  Ribo-TISH file.
- `configs/replicates.tsv` is the six-column template for sequential,
  failure-isolated replicate batch processing.

For new analyses, prefer the repository-level `miao_orf.py` unified entry
point and replace all sample-specific paths. The shell examples use
`MIAO_ORF_DATA_ROOT` and `MIAO_ORF_RESULTS_ROOT`; neither has a machine-specific
default.
