# Focus diagnosis: CYP1A2 and CYP2D6

Run `20260912T082119640472Z-88efb7f8` compares the submitted five-component OOF ensemble, each component, four earlier controls, and the interval-trained direct neural predictor. All 6,384 rows, labels, confidence bounds and parent folds match exactly. Source: `scripts/focus_diagnosis.py`; configuration: `configs/focus-diagnosis/saved-oof-v1.json`. No model training or submission was performed.

## Main finding

The saved predictions support a representation/resolution limitation more strongly than a large global calibration error. This describes challenge OOF data; it does not establish the cause of the larger 2D6 leaderboard gap.

| Endpoint | Final ST-RAE | Original four-way ST-RAE | Prediction SD | Truth SD | Mean prediction minus truth |
|---|---:|---:|---:|---:|---:|
| 1A2 | 0.683335 | 0.735039 | 0.558 | 0.990 | +0.0615 |
| 2D6 | 0.834741 | 0.866655 | 0.366 | 0.914 | +0.0272 |

Truth-conditioned extremes show large errors: the bottom true-pIC50 quintile is overpredicted by 1.077 (1A2) and 1.017 (2D6); the top quintile is underpredicted by 0.669 and 0.850. However, this alone does **not** justify expanding the prediction range. When conditioning on prediction instead, quintile mean biases range from +0.156 to -0.005 (1A2) and +0.103 to -0.089 (2D6), while truth SD *within* a predicted quintile is still 0.55–1.10 and 0.68–1.03. Compounds with similar predictions have widely different labels. A scalar mapping cannot resolve their ordering.

The lowest predicted quintile has the largest MAE: 0.863 (1A2) and 0.704 (2D6). This is a useful stratum for checking whether new representations distinguish weak-potency compounds better. Do not interpret the per-truth-bin ST-RAE in `bins.csv` as comparable endpoint performance: each narrow bin has its own often tiny denominator. Use bin MAE, bias, and spread there.

## Model diversity and fold stability

The five final components' residual correlations are 0.904–0.943 for 1A2 and 0.932–0.975 for 2D6. Their prediction correlations are 0.838–0.905 and 0.786–0.886. Averaging helped, but these components still tend to miss the same compounds, especially 2D6.

Every final component is worse standalone than the final ensemble on both focus targets. Interval direct training is also worse: 0.735342 (1A2) and 0.859393 (2D6). Merely selecting one current model per target does not explain a large unrealized gain.

Five-fold final ST-RAE ranges are 0.652–0.763 for 1A2 and 0.806–0.868 for 2D6. 2D6 is weak across all folds rather than suffering one isolated broken fold. Neighbor scores are 0.726437 (1A2, n=78) and 0.832878 (2D6, n=68); non-neighbor scores are 0.678676 and 0.842779. The saved neighbor flag does not recreate the 2D6 LB deterioration.

## Exploratory calibration control

For each endpoint and held-out fold, fit interval absolute-error calibration on predictions from the other four OOF folds. Offset is bounded to [-0.5, 0.5]; affine additionally allows slope [0.5, 1.5]. No held-out-fold labels directly enter that calibration fit.

| Endpoint | Uncorrected | Offset | Affine | Neighbor uncorrected | Neighbor affine |
|---|---:|---:|---:|---:|---:|
| 1A2 | 0.683335 | 0.681016 | 0.679008 | 0.726437 | 0.743926 |
| 2D6 | 0.834741 | 0.834721 | 0.832580 | 0.832878 | 0.820197 |

This is an exploratory screen, **not fully nested end-to-end calibration validation**: predictions used to fit a held-out fold's calibration come from other OOF models whose training sets may include that fold. A production adoption would require genuinely inner-cross-fitted calibration. The tiny gains already do not justify prioritizing this over representation changes. Learned 2D6 affine slopes are approximately 0.897–0.910, not an expansion.

## Next controlled experiments

1. Compare target-specific training with focus-weighted multitask learning, keeping other CYP and classification auxiliaries when they help. Evaluate both focus targets separately, including low-prediction strata; do not select solely on macro improvement driven by 2C9/3A4.
2. Seek external supervision that supplies new within-target resolution. Preserve assay identity and compare rank/class/continuous objectives; the present component residual correlations argue for different information, not simply more near-identical readouts.
3. Retain the current equal ensemble as the reference. Test new model residual diversity and equal additions; a standalone gain or complementary errors are needed to justify selection.

## Saved evidence

Under `artifacts/runs/regression/20260912T082119640472Z-88efb7f8/outputs/`:

- `metrics.csv`: endpoint, fold, neighborhood, error and distribution summaries.
- `bins.csv`: truth- and prediction-quantile diagnostics.
- `correlations.csv`: prediction and residual correlations.
- `calibration-coefficients.csv`, `calibrated-oof.parquet`: all fits and held-out predictions.
- `aligned-oof.parquet`, `limitations.json`: joined source predictions and methodology qualification.
