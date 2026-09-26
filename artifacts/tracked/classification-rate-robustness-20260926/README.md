# Classification call-rate robustness, 2026-09-26

Compact outputs from recorded runs:

- D6 decision audit: `classification/20260926T144420289486Z-14270294`.
- A4 static-feature proxy fit: `classification/20260926T144248722886Z-fe484e91`.
- A4 proxy decision audit: `classification/20260926T144402127052Z-20212326`.

The `*-rate-curve.csv` files contain complete-pool call-rate sweeps and paired
group-bootstrap uncertainty relative to the submitted count. The
`*-outer-fold-choices.csv` and `*-nested-policy.csv` files show rate selection
outside each evaluated fold. `a4-proxy-all-pool-oof.parquet` has one fold-excluded
score per compound, per UMAP/scaffold split, with null labels retained. The A4
proxy is a new LightGBM readout, not a reconstruction of the submitted equal-seven
ensemble.

See [the analysis](../../../docs/analysis/2026-09-26-classification-rate-robustness.md).
