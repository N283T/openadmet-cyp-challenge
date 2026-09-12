# Nonlinear calibration of saved native and global OOF

Run `20260912T113209583937Z-d94acc48` evaluated 14 fixed calibration methods on each of the native community direct mean and the previous global ensemble. Each method produces all 6,384 OOF endpoint predictions; the complete experiment contains 178,752 predictions. Source predictions are from `0c4a365a`.

## Method and evaluation

For each endpoint and held-out fold, the calibrator is fitted using predictions and labels from the other four saved OOF folds. No method or hyperparameter is selected on the evaluated fold. This is a development comparison using reused OOF, not a fully nested refit of the underlying molecular models. All fixed methods are reported separately; the lowest observed score is descriptive, not a prospectively selected winner.

The fixed methods are identity; isotonic least squares and its half-strength blend with identity; 5- and 9-knot quantile transport and their half-strength blends; monotonic piecewise-linear interval-L1 maps with 5 or 9 knots and identity penalties 0, 0.05 or 0.2; and a 5-knot point-L1 map with penalty 0.05. Knot locations are training-prediction quantiles. Predictions beyond the outer knots are clipped to the fitted endpoint values.

The interval objective minimizes average distance outside each observation's confidence interval. A linear program solves for nondecreasing knot values; the optional regularizer is the average absolute deviation of knot values from the identity map, multiplied by its fixed penalty. This uses the actual confidence intervals rather than treating every central pIC50 estimate as exact. Point-L1 uses the same solver with zero-width intervals. Quantile transport maps the training prediction marginal distribution to the training label marginal distribution; it does not use held-out labels.

## CYP2D6 results

| Fixed method | Native mean ST-RAE | Previous global ST-RAE |
|---|---:|---:|
| Identity | 0.823983 | 0.825436 |
| Isotonic L2 | 0.865352 | 0.876555 |
| Isotonic L2, half strength | 0.838963 | 0.842721 |
| Quantile 5 | 0.996868 | 1.128200 |
| Quantile 9 | 0.985051 | 1.051257 |
| Quantile 5, half strength | 0.876155 | 0.921506 |
| Quantile 9, half strength | 0.872241 | 0.893225 |
| Interval L1, 5 knots, penalty 0 | 0.820788 | 0.821403 |
| Interval L1, 5 knots, penalty .05 | 0.821968 | 0.824516 |
| Interval L1, 5 knots, penalty .2 | 0.822684 | 0.825124 |
| Interval L1, 9 knots, penalty 0 | 0.820352 | 0.820871 |
| Interval L1, 9 knots, penalty .05 | 0.819776 | 0.822711 |
| Interval L1, 9 knots, penalty .2 | 0.820699 | 0.823251 |
| Point L1, 5 knots, penalty .05 | 0.823904 | 0.824445 |

The interval-aware maps yield small D6 improvements. Isotonic least squares and marginal quantile expansion worsen the score, including half-strength variants. The apparent improvement is not a gain in molecular ranking: native D6 Spearman changes from 0.472811 to 0.466320 for the 9-knot/.05 interval map. Fold-specific monotonic maps and flat segments can alter aggregate cross-fold rank statistics.

## Four-endpoint behavior

A single fixed 5-knot, unpenalized interval map improves the native mean macro ST-RAE from 0.622574 to 0.601682. Applied to the previous global predictions, the same method moves 0.582816 to 0.581381.

| Endpoint | Native identity | Native interval map | Global identity | Global interval map |
|---|---:|---:|---:|---:|
| 1A2 | 0.695111 | 0.678298 | 0.670794 | 0.675778 |
| 2C9 | 0.528750 | 0.497859 | 0.443453 | 0.442416 |
| 2D6 | 0.823983 | 0.820788 | 0.825436 | 0.821403 |
| 3A4 | 0.442451 | 0.409781 | 0.391582 | 0.385928 |

These corrections specifically improve the interval-based objective. For example, native 3A4 point MAE increases from 0.479007 to 0.528475 while ST-RAE improves; mean prediction bias increases from +0.019329 to +0.203182. The interval map is therefore not a general improvement in central-value accuracy. It exploits the observed interval geometry, which the selected metric explicitly rewards. The original predictions and every corrected variant remain separately available.

## Validation and artifacts

Recorded audit `e43fa269` verifies all 178,752 canonical keys, labels and confidence bounds, exact identity predictions, all 40 calibration train/validation partitions, monotonicity of all 400 fitted maps, and exact reconstruction of mapped predictions from saved knots. Maximum reconstruction error is zero. No submission was performed. Full-OOF test inference is recorded separately below.

Run outputs under `artifacts/runs/regression/20260912T113209583937Z-d94acc48/outputs/`:

- `calibrated-oof.parquet`: all fixed predictions with `base_method` and `calibration` columns.
- `metrics.csv`, `macro-metrics.csv`, `fold-metrics.csv`: all candidate scores, point/interval MAE, bias and Spearman.
- `maps.json`: every fitted monotonic map's knots and values.
- `calibration-membership.parquet`: keyed training/validation membership.

Implementation: `scripts/scp_nonlinear_calibration.py`; audit: `scripts/scp_nonlinear_calibration_audit.py`; configuration: `configs/scp-calibration/nonlinear.json`. Both scripts pass Ruff. The independent calibration auditor is assessing paired uncertainty and comparing these maps with the root agent's affine alternatives.

## Fixed full-OOF test maps

At the parent's request, run `20260912T113543517968Z-5d73ca01` fits two explicitly fixed native-mean interval maps using all available OOF labels: `interval_l1_k5_reg0` and `interval_l1_k9_reg0.05`. Alongside identity, `calibrated-test.parquet` contains 3,000 rows per method (750 canonical molecules × four endpoints), 9,000 total rows. `maps.json` stores all eight endpoint/method maps and `distribution.csv` records the resulting test distributions. No test labels enter fitting.

These are exploratory local candidates requested after the development comparison, not methods independently selected before that comparison. They do not replace any existing candidate or submission. The first inference attempt `76791782` failed while concatenating Float32 identity predictions with Float64 mapped predictions; the successful new run explicitly promotes input predictions to Float64. Source values and completed runs were preserved. Inference code: `scripts/scp_nonlinear_calibration_test.py`, which also passes Ruff.
