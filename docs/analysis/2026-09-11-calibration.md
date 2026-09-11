# Regression calibration without external activity data

The first screen supports a modest interval-aware calibration gain, not a claim
that calibration closes the leaderboard gap. Fixed-OOF secondary cross-fitting
moves macro ST-RAE from **0.724704 to 0.706767**. A separate check, fitting and
scoring corrections on disjoint cluster halves inside each base-model holdout,
also improves with a half-strength correction. No new leaderboard score was
measured, and the submitted baseline remains unchanged.

## Prior PXR evidence and adaptation

Read the prior repository's `docs/deepresearch/pxr_challenge_regression_report_2026-04-19.md`,
`track1_activity/src/importance_weights.py`, calibration scripts, and the model
report's section 6. The Deep Research report discusses uncertainty calibration;
it does not by itself establish a point-prediction correction for CYP ST-RAE.
The implemented PXR method and final report provide the practical starting point:
positive affine / isotonic regression and train-to-test density-ratio weights.
PXR's reported public MAE gain of approximately 0.441 to 0.408 is source evidence,
not an expected CYP improvement. Source files are snapshotted in the screen run.

CYP adaptations:

- Fit separate maps for all four direct-inhibition pIC50 endpoints.
- Compare MSE affine, MAE affine, interval-loss affine, interval-loss shift,
  isotonic regression and monotone interpolation through five prediction bins.
  Each is tested with and without importance weights: 12 corrections plus identity.
- Interval loss is distance outside each label's confidence bounds. Minimize its
  weighted sum by linear programming; the ST-RAE denominator is constant during
  fitting. Positive affine slopes are bounded to [0.05, 3], intercept to [-10, 10].
  Shift-only fixes the slope at one. MSE uses nonnegative slope regression.
- Estimate domain probabilities with Morgan2048/radius2 and logistic regression
  C=0.1. Unlike the prior implementation's in-sample domain probabilities, use
  five connectivity-group cross-fits, correct for training domain priors, clip
  density ratios to [1/3, 3], then normalize. The domain split is not an alternative
  activity-evaluation split and never uses activity labels.
- Continue using fixed UMAP activity folds. Remove calibration-training structures
  sharing connectivity with calibration-evaluation rows.
- Do not use ChEMBL/PubChem activity data, classification targets or LB feedback
  to choose coefficients.

## Screen results

| Correction | Macro ST-RAE | Macro MAE |
|---|---:|---:|
| Identity | 0.724704 | 0.559647 |
| Interval affine | **0.706767** | 0.570936 |
| Interval affine + importance | 0.713739 | 0.576427 |
| Interval shift | 0.717349 | 0.564467 |
| MAE affine + importance | 0.721028 | 0.562537 |
| MAE affine | 0.724643 | 0.559135 |
| MSE affine + importance | 0.732610 | 0.567316 |
| MSE affine | 0.733157 | 0.562519 |
| Isotonic | 0.744232 | 0.567122 |
| Isotonic + importance | 0.754128 | 0.576316 |

The five-bin variants score 0.732693 (unweighted) and 0.731525 (weighted).
Weighted interval shift scores 0.723528. All results are in `macro.csv`.

| Endpoint | Identity ST-RAE | Interval affine ST-RAE |
|---|---:|---:|
| CYP1A2 | 0.816778 | 0.797587 |
| CYP2C9 | 0.625471 | 0.612162 |
| CYP2D6 | 0.928701 | 0.923159 |
| CYP3A4 | 0.527865 | 0.494159 |

The interval objective improves the scored metric while increasing MAE, especially
on CYP3A4. Pooled Spearman changes from 0.568038 to 0.561077 because each fold uses
a different fitted map; a single positive affine map applied to the test endpoint
preserves its ranks. The finalizer verifies this explicitly.

A shared UMAP-cluster bootstrap with 1,000 draws gives a macro interval-affine
minus identity ST-RAE delta interval of [-0.026535, -0.009965]. This is conditional
on the fixed OOF predictions and fitted corrections. It excludes fitting/seed and
method-selection uncertainty. CYP2D6's interval [-0.019584, +0.009963] crosses zero.

**Evaluation limitation:** other-fold OOF base models were trained using the
labels of the held calibration fold. Second-stage cross-fitting on these saved
predictions is therefore **not fully nested base-model CV**. This screen must not
be promoted as an unbiased deployment estimate. Full-OOF fitted-and-scored maps
are saved separately in `same_sample_fit_diagnostic.csv`; even the optimally fitted
interval-affine map has only about 0.7015 macro ST-RAE on these OOF rows. That is a
same-sample diagnostic within this bounded affine family, not a universal bound
on calibration or on hidden-test performance.

## Additional label-isolation check

For each existing base-model holdout fold, divide its UMAP clusters into two
halves. Fit calibration using one half and evaluate the other, then reverse.
Both sets of activity labels were excluded from the corresponding base model.
Repeat five predetermined cluster partitions. This avoids the screen's indirect
base-model label path without retraining the predictor.

| Correction | Mean macro ST-RAE | Range across five partitions |
|---|---:|---:|
| Identity | 0.724704 | 0.724704–0.724704 |
| Half-strength interval affine | **0.709871** | 0.708373–0.710780 |
| Full interval affine | 0.716024 | 0.712742–0.720327 |
| Interval shift | 0.720996 | 0.719327–0.723124 |
| MAE affine | 0.732716 | 0.728768–0.737748 |
| Interval affine + importance | 0.741407 | 0.726577–0.769225 |
| MSE affine | 0.745952 | 0.743285–0.749213 |
| MSE affine + importance | 0.767229 | 0.755827–0.783695 |

Half-strength means `(raw + corrected) / 2`; it reduces coefficient variance in
small calibration sets. Its MAE is 0.561853 and Spearman 0.566324. These repeated
partitions reuse labels, so their range is a stability check, not a confidence
interval. Calibration sets are much smaller than a deployment full-OOF fit;
calibration and evaluation chemistry also come from the same outer fold. This
check is not equivalent to deployment nested CV and method selection remains
exploratory. It supports further evaluation, not automatic adoption.

## Distribution shift and reviewable test predictions

Cross-fitted train/test domain AUC is 0.908–0.935. The populations are structurally
distinguishable; this does not prove conditional label shift or identify its sign.
Clipping affects 51–59% of source weights at the lower limit and 5–6% at the upper
limit. Effective sample sizes are approximately half the labeled endpoint counts.
Importance weighting did not resolve the problem in these experiments.

The full-OOF fitted interval-affine maps actually **raise test means and narrow
spreads**. They do not estimate the hypothesized lower CYP2D6 hidden population.

| Endpoint | Raw test mean | Full corrected mean | Raw SD | Full corrected SD |
|---|---:|---:|---:|---:|
| CYP1A2 | 4.9471 | 5.0657 | 0.6014 | 0.4672 |
| CYP2C9 | 4.7447 | 4.8143 | 0.5902 | 0.5334 |
| CYP2D6 | 4.5927 | 4.6621 | 0.3119 | 0.2169 |
| CYP3A4 | 4.5514 | 4.6528 | 0.8385 | 0.6816 |

CYP2D6 predictions below 4 decrease from 10 to 0 under the full map (1 under the
half-strength map). Thus the CV correction is not evidence that the test bias
problem is fixed. Applying a large negative shift improves neither ordinary OOF
nor its interval objective; whether such a shift improves hidden test cannot be
measured without new information. Preserve that distinction rather than copying
the community author's LB-tuned coefficients.

A 9-shift × 8-scale scenario grid records 288 endpoint settings, OOF scores and
hypothetical test distribution summaries. No test scores are invented. The
finalizer's `affine-sensitivity.png` uses the actual nonuniform grid coordinates;
the earlier screen PNG used uniformly spaced display cells and is superseded for
visual interpretation. Numeric scenario outputs are unchanged.

## Artifacts, validation and decision

- Failed screen `20260911T070215152437Z-5c96e736`: Float32/Float64 concatenation
  mismatch when saving test outputs; retained. Fixed by Float64 test arithmetic.
- Completed screen `20260911T070300314219Z-bb1e7e28`:
  metrics, conditional bootstrap, keyed corrected OOF/test predictions,
  domain weights, full and fold-specific maps, scenario grid and source evidence.
- Holdout check `20260911T070421577926Z-9ac84861`:
  5 repeats × 8 methods, keyed predictions, cluster-half membership and maps.
- Final audit/registration `20260911T070559006298Z-fe3022e7`:
  `review-test-variants.parquet`, `transforms.csv`, `test-distributions.csv`,
  corrected sensitivity figure, and 3,816 metric records in PostgreSQL.

Official scoring is reproduced in both experiments; predictions carry compound,
endpoint, fold and split identifiers. The final audit verifies all 750 × 4 keys
for identity/full/half-strength maps and rank preservation. Tests check the
interval-loss optimizer against an independent dense shift search and a feasible
affine construction. `pixi run check` passed; `CYP_TEST_DB=1 pixi run test`
passed all 41 tests.

Keep the existing submission. Retain full and half-strength interval-affine
predictions as reviewable prospective candidates, without claiming a measured LB
gain. Do not adopt importance weighting or isotonic calibration from this screen.
The next unresolved question is how to validate a test-population correction:
OOF-based ST-RAE optimization and LB population placement are different problems.
External-data modeling and activity classification remain deferred.
