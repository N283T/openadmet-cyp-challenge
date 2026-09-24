# Complete-pool D6 transfer artifacts

Project-authored compact outputs from recorded classification runs:

- `c8368f50`: LightGBM complete-pool OOF.
- `2542e03a`: ExtraTrees complete-pool OOF.
- `d67281d0`: LightGBM pair decision evaluation.
- `cc83a23c`: LightGBM versus ExtraTrees decision evaluation.

Each `all-pool-oof.parquet` uses `compound_id`, endpoint, split ID, fold, model,
score, and nullable `y_true`. Each model and split has 6,146 unique molecules.
The CSV files are the final fixed-rate and selected-rate summaries. All sources,
input checksums, configurations, and full outputs remain in the recorded runs.
No hidden labels or third-party model payloads are included.

See [the analysis](../../../docs/analysis/2026-09-24-d6-full-pool-transfer.md).
