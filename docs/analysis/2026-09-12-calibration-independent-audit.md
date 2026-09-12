# Independent calibration audit and conditional uncertainty

Recorded run: `regression/20260912T113437714858Z-d1889f83`.

The native `community_direct_mean` baseline contains exactly6,384 canonical non-anchor endpoint observations, including1,457 CYP2D6 observations. Compound IDs, endpoint names, parent folds, labels and both interval bounds match the canonical measurements and parent assignments. All17 affine variants for both input baselines were reconstructed from saved per-fold slopes/intercepts; predictions agree within1e-12. Their816 endpoint/fold/pooled ST-RAE entries also agree within1e-12 with the exact declared official scorer snapshot used by the nonlinear run.

This verifies affine run `4a9dd6ab` against its inputs and coefficients. It does not treat the community article's alternative metric implementation as interchangeable with that declared scorer. The nonlinear predictions are from `d94acc48`; their2D6 identities, labels, intervals and folds match the same baseline. Independent nonlinear map reconstruction is recorded separately by the trainer owner.

## Paired CYP2D6 group bootstrap

We resampled1,185 parent-assignment groups with replacement,500 times, using seed20260912. Each replicate uses identical sampled rows for all three methods, keeps all rows belonging to each selected group, and recomputes the official scorer including its denominator. The saved predictions and calibration maps remain fixed. Negative score differences favor calibration.

| Fixed calibration | ST-RAE | Difference from identity | Conditional95% percentile interval | Fraction of replicates improved |
|---|---:|---:|---:|---:|
| Identity | 0.823983 | — | — | — |
| Offset, interval loss | 0.822792 | -0.001191 | [-0.004195, +0.001868] | 79.8% |
| Monotone interval L1,9 knots, penalty0.05 | 0.819776 | -0.004207 | [-0.009076, +0.000995] | 95.0% |

The nonlinear map has the larger observed improvement, but both conditional intervals include zero. These are development OOF comparisons among screened methods, not confirmatory intervals. The bootstrap does not repeat model fitting, calibration fitting or method selection, and therefore does not include those sources of uncertainty. The improvement fractions are descriptive bootstrap frequencies, not posterior probabilities or corrected significance tests.

Saved artifacts include `recomputed-metrics.csv`, `paired-bootstrap.parquet`, `bootstrap-summary.csv` and `verification.json`. No fitting or prediction artifacts were modified.
