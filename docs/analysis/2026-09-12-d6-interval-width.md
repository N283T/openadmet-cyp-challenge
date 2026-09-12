# CYP2D6 interval-width-conditioned point decisions

## Community placement is partly leaderboard-derived

The [current SuperCowPowers blog](https://supercowpowers.github.io/workbench/blogs/cyp_challenge/)
describes model learning and a subsequent prediction placement step separately. It
states that the ST-RAE placement was found by leaderboard sampling and identifies
learning confidence-interval widths as an unresolved offline alternative.

The linked [calibration source](https://github.com/SuperCowPowers/workbench/blob/main/ml_pipelines/OpenADMET/cyp/scripts/cyp_recalibrate.py)
confirms that `BLIND_MOMENTS` and `SOLVED_PEARSON` are recovered from three scored
submissions. Its CYP2D6 `STRAE_MOMENTS` setting (mean 3.57, standard deviation 0.90)
is based on public-board probes, not training-only cross-validation. Even its
optional OOF Pearson route multiplies by a board-derived `OOF_TO_BLIND` ratio.
These constants are not used in this experiment. The R-squared-optimal placement
and the minimum-ST-RAE placement differ; the former follows a squared-error
criterion, whereas ST-RAE permits any prediction inside the observed interval.

The fetched blog and code are preserved as `source_blog.txt` and `source_code.txt`
in the experiment outputs. This source inspection does not assert that leaderboard
placement alone explains the author's complete performance gain.

## New experiment versus prior work

Earlier work already tested affine/shift/isotonic calibration, importance weights,
TabPFN lower/upper-bound mixture medians, direct interval-L1 fitting, and structural
neighbor interval medians. This experiment instead learns a low-capacity residual
correction conditional on **predicted asymmetric interval widths**.

Run `20260912T101132109311Z-afffa016` uses frozen Jeremy600 features `c060d197` plus
four updated, outer-excluded primary-screen predictions `d272774c`, and the fixed
parent split `a3923255`. Only CYP2D6 models are trained. The 300-tree LGBM L1 base
exactly reproduces the original `1e7001de` D6 OOF predictions.

Within each outer training set, a three-way GroupKFold generates inner-held-out
base predictions and separate predictions of `log1p(y - lower)` and
`log1p(upper - y)`. Width models use 150 trees, seven leaves, minimum leaf size 50.
Base and width models are then refitted on all outer training rows. No outer
validation endpoint label or confidence bound trains any component.

Three residual models are compared:

1. Point residual: learn `y - base_prediction` from the base prediction alone.
2. Interval residual: learn the equal mixture of `lower - base_prediction` and
   `upper - base_prediction`, again from the base prediction alone.
3. Width-conditioned interval residual: use the same interval targets, adding
   predicted lower/upper log-widths, their sum and their difference as inputs.

The equal half-weight endpoint loss is algebraically equivalent to interval
absolute error up to a prediction-independent constant. Correction models use
60 trees, five leaves, minimum leaf size 80 and L2=5. Fixed 50% and 100% correction
strengths are both saved, without HPO or outer-fold weight selection. Width inputs
for correction training are inner-held-out predictions, never measured widths.
Duplicate endpoint rows change minimum-leaf effective molecule count relative to
the point residual arm; the two interval arms share exactly the same duplication.

## Results

| D6 method | ST-RAE | MAE | Spearman |
|---|---:|---:|---:|
| Matched raw LGBM | 0.842287 | 0.560907 | 0.452167 |
| Point residual, 50% | 0.843754 | 0.560543 | 0.444352 |
| Point residual, 100% | 0.851670 | 0.562817 | 0.425067 |
| Interval residual, 50% | 0.842983 | 0.560366 | 0.447974 |
| Interval residual, 100% | 0.849275 | 0.562026 | 0.434152 |
| Width-conditioned interval residual, 50% | 0.840364 | 0.560069 | 0.453231 |
| Width-conditioned interval residual, 100% | 0.845179 | 0.561908 | 0.445146 |
| Global Caruana weight-crossfit `84c28266` | 0.825436 | 0.554068 | 0.472752 |

The width-conditioned half correction gives a tiny matched gain, not a meaningful
improvement over the current ensemble. Width prediction is weak: OOF Spearman is
0.147720 for the lower-side log-width and 0.207578 for the upper-side log-width.
Thus this implementation does not reliably identify which compounds have wide
low-penalty intervals.

Independent audit `20260912T101307140418Z-e87d4248` verifies group-disjoint inner
folds, full outer-training coverage, finite predictions, and all canonical keys,
labels and confidence bounds. A 1,000-draw conditional parent-group bootstrap for
the half width-conditioned correction gives ST-RAE delta versus raw
[-0.007079, +0.002860] and versus the same-strength interval residual
[-0.007268, +0.002005]. Neither distinguishes a gain from zero. These intervals
condition on fitted models and exclude refitting and method-selection uncertainty.

The first audit `20260912T101253460579Z-ffe9028d` failed on NumPy `isfinite` applied
directly to a Polars Series; converting explicitly to NumPy fixed it. Both attempts
are retained. This did not require changing or refitting the predictions.

## Artifacts and decision

- `d6-oof.parquet`: seven methods × 1,457 D6 OOF rows.
- `canonical-oof.parquet`: eight methods × 6,384 canonical rows. Each experimental
  method changes D6 only; the other endpoints are explicitly copied from the global
  Caruana weight-crossfit reference. These are not new all-endpoint model fits.
- `width-oof.parquet`: predicted and measured log-widths for evaluation.
- `inner-memberships.parquet`: inner held-out fold membership by outer fold.
- All base, width and correction boosters are retained, as are metrics by fold and
  neighbor subset, source snapshots, configuration and verification outputs.

Retain as a small alternative candidate and negative evidence for this particular
width estimator. Do not apply the community's board-derived constants, replace the
current ensemble, or infer a hidden-test improvement. No submission was made.
