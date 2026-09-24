# Historical classification leaderboard benchmark

Finalized outputs from recorded run
`classification/20260924T031257156540Z-4febdf99` using
`configs/analysis/classification-historical-lb-inventory-v1.json`.
The source inputs are seven project-submitted boolean prediction packages and
14 official post-scoring endpoint leaderboard CSV snapshots, all preserved in
the run's input manifest. No hidden labels are included.

- `states.csv`: exact endpoint MCC, timestamp, and positive-call count per state.
- `transitions.csv`: six chronological changes per endpoint.
- `changed-calls.parquet`: molecule IDs of added and removed positive calls.
- `summary.json`: counts and directional summary.

SHA-256:

```text
86bee25264a438ac8c1c0e56cdd68893eca25a928943c632583c5adde209ef02  changed-calls.parquet
691c9e8ce13a6a697ed9bbd406edcfc1401c889e72933557d2179edf3fde0512  states.csv
5bb49001115ee42ec339d79dc064de9da4993e70d02c186421f9acf3f9c09ee6  summary.json
32d378f2419f7b5a06da14fe1a86302535ec5d3cd318bbf6c9d9957199b880f0  transitions.csv
```
