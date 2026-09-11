# Learning the CYP interval objective

This experiment follows the distribution-summary ensemble campaign. Keep the
same outer UMAP folds and TabPFN 2.6 backbone/feature recipe, without additional
external CYP data. The complete five-fold results and audit are recorded below.

## Marginal interval models

For lower and upper bounds L ≤ U, distance outside the interval satisfies

`distance(q, [L,U]) = (|q-L| + |q-U| - (U-L)) / 2`.

Conditioning on molecular features, the expected-loss minimizer is a median of
the equal mixture of the conditional lower-bound and upper-bound distributions.
For continuous marginals its first-order condition is `F_L(q) + F_U(q) = 1`.
This does not require knowing the dependence between the two predicted marginals.
The ST-RAE normalization is constant with respect to the predictions on a given
evaluation set, so it does not change this point-decision objective.

Fit separate TabPFN readouts to training-set confidence lower bounds, upper
bounds and their midpoint. Reuse the original outer-training-only Mordred
imputer/column set and PCA256 plus nested predicted-log2fc4 features. Reload the
original fitted mean model and reproduce its OOF predictions before any refit,
checking the reconstructed input matrix. Fit new targets on the same outer
training rows. Held-out/test confidence bounds are never supplied to the models.

Use the same 2.6 checkpoint, eight estimators and existing fold seed. Preserve a
new fitted state for each target and endpoint. Request marginal quantiles at
0.01, 0.02, …, 0.99. Approximate the mixture median by the median of those 198
marginal quantile values; retain the marginal arrays for inspection. This is an
approximation to the fitted distributions, not a claim of calibrated uncertainty.
Additional point recipes are midpoint-model mean, midpoint-model median and the
average of the two marginal medians. The latter is generally not the same as
the mixture median and is explicitly kept as a different recipe.

Sixty target readout fits cover four endpoints × three targets × five folds.
These are new TabPFN fits, not just transformations of the old test outputs.

## LightGBM interval training

The same identity above also permits empirical interval-loss training by creating
two copies of each training feature row, attaching L and U as labels, and giving
each copy weight 0.5. Weighted MAE then differs from interval loss only by the
constant half-mean interval width. Molecular train/holdout separation is retained.

Use LightGBM L1 regression, learning rate 0.05, 15 leaves, minimum 20 row samples
per child, L2=1, four CPU threads, deterministic column-wise fitting and seed
42+fold. Compare:

- Actual point labels, 300 trees.
- Half-weight duplicated interval endpoints, 300 trees.
- Half-weight duplicated interval endpoints, 800 trees.

No early stopping or validation-based model fitting is used. Reuse the same
feature matrices as the TabPFN comparison, rather than the earlier Mordred/Morgan
LightGBM baseline. Row duplication changes the meaning of the minimum-row-count
constraint even though summed training weight is preserved; this is not a perfect
isolation of the loss effect. No duplicate crosses a training/holdout boundary.
Save each fitted booster and keyed OOF/test predictions. These are another sixty
models, independent of the sixty TabPFN target fits.

The reused launcher contains an inherited TabPFN-oriented purpose/scope string
in its LightGBM configurations. The explicit `backend: lightgbm`, source branch,
saved booster files and method names `lgbm_point_300`, `lgbm_interval_300`,
`lgbm_interval_800` identify the actual algorithm. Do not interpret that inherited
text as a TabPFN run.

## Structure-neighborhood alternative

Morgan 2048/radius2 Tanimoto similarity, k=10/50, similarity powers 2/8. Predict
weighted medians of either observed labels or the concatenated interval endpoints.
For every outer fold, exclude its compounds and matching connectivity blocks from
the reference pool. Average the five fold-specific point predictions on test.
Compare local-only predictions and fixed 25%/50% local contributions to the
original TabPFN mean: 24 recipes, no fitted calibration coefficients.

The best observed variant, k=50/power2/interval with 25% local weight, scores
macro ST-RAE 0.713079, MAE 0.567590, Spearman 0.569979. It improves the original
mean's ST-RAE but is worse than the three-member q0.55 ensemble (0.691191), and
its MAE is worse than the original mean. Preserve it as a different structural
signal; the combined audit also checks 10%/25% contributions to q0.55.

## Evaluation and limitations

Each point recipe is scored against original pIC50 values and original confidence
bounds, not against its modified training target. Reproduce the official scorer,
check all 6,525 compound/endpoint OOF keys and all 750 × 4 test keys, and use
explicit key ordering before converting to numeric arrays. Fixed 25%/50% blends
with the q0.55 ensemble are compared without fitting their weights.

A shared UMAP-cluster bootstrap reports conditional deltas against the original
mean and q0.55 ensemble. It excludes model refitting, seed and recipe-selection
uncertainty. Choices build on earlier results on this same dataset. No new hidden
labels or leaderboard scores are used or inferred.

Unit checks verify the mixture-median action against an independent loss grid,
the endpoint-MAE/interval-loss identity, and shuffled-key alignment. Project checks
pass and DB-enabled tests pass (44 tests; dependency deprecation warnings retained).

## Completed five-fold results

| Recipe | Macro ST-RAE | MAE | Spearman |
|---|---:|---:|---:|
| Original TabPFN mean | 0.724704 | 0.559647 | 0.568038 |
| Three-member q0.55 | 0.691191 | 0.550526 | 0.581232 |
| TabPFN interval-mixture median | 0.695389 | 0.569197 | 0.569916 |
| TabPFN midpoint mean | 0.733668 | 0.560500 | 0.569094 |
| TabPFN midpoint median | 0.718252 | 0.554284 | 0.570894 |
| Average of bound medians | 0.718516 | 0.555046 | 0.570372 |
| LightGBM point / 300 | 0.720104 | 0.560676 | 0.556580 |
| LightGBM interval / 300 | 0.719425 | 0.584612 | 0.555072 |
| LightGBM interval / 800 | 0.726148 | 0.588415 | 0.553003 |
| q0.55 75% + interval LightGBM 25% | **0.688853** | 0.554293 | 0.580535 |
| q0.55 75% + interval TabPFN 25% | 0.689186 | 0.553163 | 0.580345 |
| q0.55 50% + interval TabPFN 50% | 0.689081 | 0.557327 | 0.578405 |
| q0.55 90% + local interval median 10% | 0.689757 | 0.553805 | 0.582197 |

The interval-mixture decision improves substantially over the original mean,
whereas predicting the interval midpoint is not an effective replacement.
LightGBM interval training does not materially improve its standalone point-label
control and costs MAE. Its small ensemble contribution helps primarily CYP3A4.
The 800-tree model does not rescue standalone interval LightGBM.

The best numerical blend's macro delta interval versus q0.55 is
[-0.005289, +0.000743], crossing zero. Its delta versus the original mean is
[-0.048642, -0.024550]. The 25% interval-TabPFN blend's incremental interval is
[-0.003942, -0.000043], barely below zero before any method-selection adjustment.
Do not present 0.689 as an established improvement over the much simpler 0.691
recipe. The extra alternatives provide evidence and prospective candidates.

CYP2D6 remains unresolved on test: original mean 4.5927, q0.55 4.7392, standalone
interval-mixture 4.6830, and the 25% LightGBM blend 4.7403. All three new recipes
have zero predictions below 4. The new interval objective does not independently
identify the suspected hidden-test label shift.

## Completed runs and local artifacts

TabPFN lower/upper/midpoint runs (fold order):

- `20260911T072716421566Z-24ff1fd7`
- `20260911T073020399698Z-6f17daf6`
- `20260911T073322300653Z-61a1009a`
- `20260911T073628329393Z-94565589`
- `20260911T073931300085Z-8714b7fb`

LightGBM point/interval runs (fold order):

- `20260911T074009277562Z-be7cf1d4`
- `20260911T074131117243Z-9fd947df`
- `20260911T074301803825Z-4989e98c`
- `20260911T074417567845Z-70528e07`
- `20260911T074521296699Z-92454c6a`

Neighborhood run: `20260911T073458027748Z-50193fa7`.
TabPFN/neighbor audit: `20260911T074407984282Z-c3a37c3b`.
Combined audit: `20260911T074638497859Z-7f0af974`.
Final packaging and 3,528 metric-row registrations:
`20260911T074800083392Z-1855a7d1`.

All 120 model fits completed. The combined audit checks seven readout recipes and
fixed blends, reproduces official scoring, and preserves keyed OOF/test outputs.
The packaging run saves standalone interval and blended predictions, validates
750 official test rows per candidate, and produces distribution statistics/plots.

Local immutable alternatives:

- q0.55 + 25% interval LightGBM: `20260911T074911136468Z-1ab5c0a3`.
- q0.55 + 25% interval TabPFN: `20260911T074911361289Z-aaf5da60`.

Keep the simple three-member q0.55 recipe as the main prospective reference and
these alternatives for later evaluation. Preserve the original submitted baseline.
No new leaderboard upload occurred. Further local calibration cannot establish a
hidden population correction without additional evidence; a small CV difference
is not sufficient reason to keep enlarging this ensemble.
