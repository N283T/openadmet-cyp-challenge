# CheMeleon depth probe

Extraction run: `shared/20260925T065503727832Z-1c269e23`.
Readout run: `shared/20260925T065626160287Z-8017d516`.

- `extraction-summary.json`: depths, dimensions and final-path equivalence.
- `metrics.csv`: pooled per-endpoint metrics.
- `oof.parquet`: keyed held-fold predictions with track, endpoint, depth,
  compound ID, split ID and fold.
- `blend-metrics.csv` and `oof-blends.parquet`: fixed equal-mean OOF comparison
  from `shared/20260925T070604392889Z-6f643137`.

See [analysis](../../../docs/analysis/2026-09-25-chemeleon-intermediate-depth.md).
