# MIAO 1.0.0 release contract

Version 1.0.0 is the first public MIAO release. Internal labels such as
1.2 and 1.3 were unpublished development snapshots; they are retained only as
development history and do not define the public semantic-version sequence.

The 1.0.0 public source contract was finalized on 2026-09-04.

The release freezes the statistical engine and formal thresholds certified by
the 520-pair, dual-background spike-in benchmark. The machine-readable lock is
`config/releases/benchmark_520_dual_v1.0.0.json`.

The public 1.0 output/provenance contract includes:

- observed P-sites across the full ORF region and analyzed core;
- observed phase 0/1/2 counts and fractions in intORF coordinates;
- fitted host-CDS and intORF-like component fractions and expected read counts;
- a bounded 95% likelihood-ratio profile interval for `lambda_hat`;
- finite p-value and FDR evidence scores;
- compact QC status, flags, filter reason and background-template sensitivity;
- an atomic `<out>.run_manifest.json` containing version, parameters, runtime,
  input identities, source hash, output schema and field semantics.
- accurate fixed-size importance sampling reports the actual simulated draw
  count in `importance_total_reps_used`; exact zero-LRT shortcuts report zero
  simulated draws.

The component read counts are fractional fitted expectations. They are not
hard read-level assignments because reads in the overlapping interval usually
cannot be uniquely attributed to the host CDS or intORF.

The benchmark artifacts remain identified by their truth and result SHA-256
hashes in the release lock. Future changes to the statistical engine or call
gates require a new engine ID and a new benchmark certification.
