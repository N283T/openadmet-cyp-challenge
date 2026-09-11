# Regression distribution readouts and ensembles

Without additional ChEMBL/PubChem activity labels, an equal three-member TabPFN
2.6 blend using the 55th percentile improves fixed UMAP OOF macro ST-RAE from
**0.724704 to 0.691191**, MAE from 0.559647 to 0.550526, and Spearman from
0.568038 to 0.581232. This is an exploratory CV result, not a measured leaderboard
gain. Saved test predictions still concentrate CYP2D6 near high pIC50 values.

## Scope and methods

After the first calibration screen, the user explicitly requested continued
experimentation. Keep regression as the active task and retain TabPFN 2.6 weights.
Reuse the fixed UMAP folds and existing challenge-only auxiliary learning.

The initial expanded screen compares 20 methods:

- Identity, half-strength interval affine, equal top-three and equal all-five
  predictions, and four 50:50 baseline/alternative blends.
- Convex interval-loss ensembles with penalties 0.01/0.1/1.0 for weight moved away
  from the baseline; nonnegative weights sum to one and a shift is bounded to ±0.5.
- Residual correction from prediction-only context (prediction, its square, and
  hinges at 4.5/5.5), and chemical context (prediction plus molecular weight,
  logP, TPSA, donor/acceptor/rotatable-bond/ring counts, fraction sp3).
- Context features are standardized on calibration-training rows and clipped
  to ±5 standard deviations. Smooth interval loss uses temperature 0.03 and
  quadratic penalties 0.01/0.1/1.0; chemical MSE ridge is a comparison. Coefficient
  bounds for interval correction are ±1. The exact official interval score,
  rather than the smoothing surrogate, is used for evaluation.

Fitted methods use the previous five repeated cluster-half partitions inside
common base-model holdout folds, with connectivity exclusion. Both fitting and
scoring labels are unseen by that base model. This is a limited-data robustness
check, not deployment-equivalent nested CV. Fixed blends require no fitted
correction; their scores are ordinary OOF, replicated unchanged in the summary.

The five source members are baseline stacked Mordred TabPFN, compact PCA256 +
predicted-log2fc4 TabPFN, fold-specific log2fc-pretrained 2048-dimensional embedding
TabPFN, Mordred-only TabPFN, and direct CheMeleon FT. None adds the newly curated
external CYP labels. All 6,525 compound/endpoint keys, folds and labels match.

## Context and blend screen

| Method | Macro ST-RAE | MAE | Spearman |
|---|---:|---:|---:|
| Baseline | 0.724704 | 0.559647 | 0.568038 |
| Three TabPFN members, equal mean | 0.712130 | 0.554684 | 0.579980 |
| Baseline + compact, 50:50 | 0.716510 | 0.556704 | 0.574005 |
| Baseline + pretrained, 50:50 | 0.717935 | 0.556680 | 0.577230 |
| All five, equal mean | 0.722520 | 0.560657 | 0.571724 |
| Prediction-context interval correction, penalty 1 | 0.714256 | 0.560347 | 0.565004 |
| Chemical-context interval correction, penalty 1 | 0.718208 | 0.561954 | 0.559828 |
| Convex interval ensemble, penalty 0.01 | 0.719278 | 0.566656 | 0.563767 |
| Half-strength global interval affine | 0.709871 | 0.561853 | 0.566324 |

Adding Mordred-only or direct-FT predictions at 50% worsens ST-RAE to 0.734362 or
0.750013. More flexible chemical corrections overfit: interval penalty 0.01
scores 0.769056; MSE penalty 0.01 scores 0.816713. Retain these negative results.
The useful direction is simple diversity among the three stronger readouts,
rather than a fitted context model or fitted many-member weights.

## Distribution summaries from unchanged TabPFN fits

Reload all five fold models for each retained member. Remap only the copied
archive's base-weight path; preserve original states. Reconstruct the same keyed
OOF inputs and verify mean predictions against their saved values at 1e-5 absolute
and relative tolerance before producing test distributions. No model refit or
new seed is involved in this phase.

TabPFN's `predict(output_type="main")` returns mean, median, mode and requested
quantiles. Initial baseline extraction saves 12 quantiles, mean, median and mode.
The common comparison uses mean, median and q0.40/0.45/0.50/0.55/0.60/0.65/0.70,
plus the 50:50 mean/median midpoint. The duplicate median/q0.50 is preserved as an
implementation check, not independent evidence. The initial mode and extreme
quantile outputs remain saved but are outside the common comparison.

Test predictions average per-fold point summaries and, for blends, then average
across member types. **An average of fold quantiles is not the quantile of their
mixture distribution.** These are explicit point-prediction recipes.

| Recipe | Macro ST-RAE | MAE | Spearman |
|---|---:|---:|---:|
| Baseline mean | 0.724704 | 0.559647 | 0.568038 |
| Baseline median | 0.708706 | 0.552281 | 0.570013 |
| Baseline q0.55 | 0.702744 | 0.555911 | 0.570012 |
| Three-member mean | 0.712130 | 0.554684 | 0.579980 |
| Three-member median | 0.698872 | **0.547867** | 0.581126 |
| Three-member q0.55 | **0.691191** | 0.550526 | **0.581232** |
| Three-member q0.60 | 0.699709 | 0.560872 | 0.581143 |

| Endpoint | Baseline mean ST-RAE | Three-member q0.55 ST-RAE | q0.55 Spearman |
|---|---:|---:|---:|
| CYP1A2 | 0.816778 | 0.782348 | 0.514177 |
| CYP2C9 | 0.625471 | 0.595658 | 0.658473 |
| CYP2D6 | 0.928701 | 0.887998 | 0.383034 |
| CYP3A4 | 0.527865 | 0.498761 | 0.769243 |

The common UMAP-cluster bootstrap (1,000 draws) gives fixed three-member q0.55
minus baseline mean macro ST-RAE delta interval [-0.045835, -0.022301]. It is
conditional on trained models and the recipe, and does not account for choosing
members/quantiles after examining this CV or for refitting uncertainty.

## Quantile choice and correction robustness

Select a summary using only one cluster half of each base-model holdout, then
score the opposite half. Repeat the same five partitions. Three-member quantile
selection scores mean 0.698654 (range 0.697333–0.700222); blending the selected
summary 50:50 with the mean scores 0.697373 (0.696204–0.698534). This is less
optimistic than selecting q0.55 from all OOF labels and reporting its score.
Member selection still used this dataset, and calibration sets are smaller and
more same-fold representative than deployment.

Three-member median plus half-strength interval affine scores 0.695092
(0.693479–0.697751); full-strength correction worsens to 0.704482. It does not
beat the best observed fixed quantile. These partition ranges reuse labels and
are not independent confidence intervals.

## Test populations and local candidates

The three-member q0.55 test means are 5.0860 / 4.8435 / 4.7392 / 4.6815 for
1A2/2C9/2D6/3A4. CYP2D6 SD is 0.2437 and no prediction is below 4, versus the
original baseline mean 4.5927, SD 0.3119 and 10 below 4. **The CV gain does not
resolve the hypothesized low-pIC50 test shift.** No hidden test score is available.

Two local immutable candidates were validated against all 750 official test IDs,
SMILES and four regression columns:

- q0.55: `20260911T073131924313Z-6d6f5daf`.
- median alternative: `20260911T073132254302Z-da91a992`.

They are prospective artifacts, not new remote submissions. The previously
submitted baseline remains intact.

## Runs and audit trail

- Context/blend screen: `20260911T071844846604Z-595cb35a`.
- Baseline distribution inference folds 0–4:
  `20260911T071642514941Z-46dba429`, `20260911T071734625635Z-c0a1c0ca`,
  `20260911T071821705444Z-d24bb710`, `20260911T071909982887Z-282867eb`,
  `20260911T071957777491Z-48bf1b6b`.
- Compact/pretrained inference folds 0–4:
  `20260911T072106821615Z-39d5ba1d`, `20260911T072215746368Z-c51513a0`,
  `20260911T072322810219Z-9301a74e`, `20260911T072430468552Z-c902f32b`,
  `20260911T072538070064Z-d344ca0e`.
- Invalid initial comparison: `20260911T072325389052Z-0f774114`. Pivot/join output
  was not explicitly reordered before positional scoring, producing a wrong
  baseline metric. Discarded; outputs preserved. GPU extraction is unaffected
  and had independently reproduced original keyed OOF values.
- Corrected baseline comparison: `20260911T072502653437Z-d4ad44c5`.
- Combined comparison: `20260911T072717606484Z-8ec7593a`.
- Candidate packaging and metric registration: `20260911T073034292947Z-1b45d9b9`.

The corrected comparisons explicitly align/check molecule keys and reproduce
baseline mean OOF before scoring. A unit test deliberately shuffles keys and
checks missing predictions fail. Both comparisons reproduce official ST-RAE.
All source configurations, weights, fitted states, OOF/test outputs, bootstrap
results and rejected methods remain recorded. The packaging run registers both
context and distribution metrics in PostgreSQL and writes test-distribution plots.
