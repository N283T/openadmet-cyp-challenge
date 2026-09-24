# D6 AP-mass90 complete-pool replay

The parquet files contain all 6,146 held-fold molecules on UMAP and scaffold
for the seven-member AP compositions and two September 15 control analogs.
`y_true` is null for unscored molecules. The CSVs compare AP with both its
leading ExtraTrees atom and the September 15 raw-score analog at the historical
49/750 rate and with rate selection on other outer folds. The density table
tests a fixed cross-fitted train-to-blind weighting; it is a negative result.
`extra-vs-sep15-fixed-rate.csv` and `extra-vs-sep15-nested.csv` are the direct
leading-atom check from run `25c958cf`; `label-route-summary.csv` records the
fixed-rate AP/control swap counts by training assay route from run `c55d672a`.
`scoring-mask-sensitivity.csv` checks a hypothetical shift-interval mask from
run `622cb62b`; it is not the known official scoring mask.
`control-ap25-fixed-rate.csv` and `control-ap25-nested.csv` hold the conservative
blend decision from run `2fe73ea9`.
`control-ap25-blind-swaps.csv` lists the eight changed calls in unsubmitted
candidate `20260924T043703823086Z-51106137`.
See
[the analysis](../../../docs/analysis/2026-09-24-d6-ap-full-pool-replay.md)
for interpretation and immutable run IDs.
