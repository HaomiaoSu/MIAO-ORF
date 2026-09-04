# Contributing to MIAO-ORF

Bug reports and focused improvements are welcome through GitHub Issues and pull requests.

Before opening a pull request:

1. create a branch from `main`;
2. keep scientific-behavior changes separate from documentation or formatting changes;
3. add or update tests for changed behavior;
4. run `python -m unittest discover -s tests -v`;
5. describe any effect on output schemas, statistical thresholds or reproducibility.

Changes to the DM model, evidence gates, default thresholds or release-lock parameters require explicit scientific justification and validation. Do not commit patient-identifying data, controlled-access data, large sequencing files or generated analysis outputs.
