# CYP2D6 regression: internal primary-response tests

Three recorded experiments examine whether the challenge's own CYP2D6
single-concentration response can supply a materially different regression
signal. All comparisons use compound-excluded parent-series folds; the
supervision and primary-model screens also use the prespecified
parent-retained scaffold folds. Anchors remain in training. These are
exploratory comparisons on development folds, not independent confirmation.
No hidden test labels or leaderboard numbers enter a fit.

## 1. Measured response is informative, but its prediction loses much of that information

Run `regression/20260926T031653831897Z-65615961` fits a decreasing
measured-primary-to-direct mapping on each outer-training fold and evaluates
held-out D6 direct labels. The measured-response arm is a **nondeployable
oracle diagnostic**: the test molecules have no released primary measurements.
A separate arm fits the same one-dimensional mapping using fold-excluded
predicted primary responses for training and validation. The control is the
current submitted D6 model's unplaced OOF on the same 1,457 nonanchor rows.

| D6 parent OOF method | ST-RAE | MAE | Spearman |
| --- | ---: | ---: | ---: |
| Current submitted control | .8284 | .5549 | .4753 |
| Measured primary, linear mapping | **.6049** | **.4542** | **.8255** |
| Measured primary, isotonic mapping | .6074 | .4559 | .8113 |
| Predicted primary, linear mapping | .8871 | .5774 | .4288 |
| Predicted primary, isotonic mapping | .8494 | .5601 | .4430 |

Across the five held folds, measured primary has D6 direct Spearman magnitude
`.776–.855`, while the fold-excluded primary predictions have only
`.371–.474`. Primary prediction RMSE on the 1,457 paired direct-label rows is
`1.055–1.217` log2FC units. The oracle comparison identifies information in
the assay measurement; it does not establish an achievable blind-test gain or
that a one-dimensional model is otherwise sufficient.

The first diagnostic attempt `regression/20260926T031635816254Z-44876f06`
failed at final table concatenation because its control columns had a different
order. It produced no scientific result, was marked discarded, and was corrected
before the completed run above.

## 2. Primary-only proxy labels damage direct regression

Run `regression/20260926T031325546569Z-4f1b6599` uses the 2,883 compounds
with released D6 primary values but no direct-pIC50 label as a possible weak
supervision pool. Within each outer fold, a linear or monotone isotonic map is
fitted only on training compounds with **both** measurements. The mapped
primary-only compounds receive sample weights 0.25 or 1.0 in an otherwise
matched 647-feature LightGBM direct regressor. All held-out compounds are
excluded from both direct and proxy training. The direct-only model is fitted
with exactly the same readout and features.

| Method | Parent ST-RAE | Scaffold ST-RAE | Parent weak-label MAE |
| --- | ---: | ---: | ---: |
| Direct-only control | **.8356** | **.8272** | 1.9383 |
| Isotonic proxy, weight 0.25 | 1.0062 | 1.0036 | 1.6366 |
| Isotonic proxy, weight 1 | 1.2777 | 1.2807 | 1.3870 |
| Linear proxy, weight 0.25 | 1.3265 | 1.3111 | 1.4123 |
| Linear proxy, weight 1 | 2.0317 | 2.0315 | 1.0631 |

The primary-only cohort has mean/median log2FC `+.307/+.376`, versus
`−1.797/−1.738` among the direct-labelled cohort. The mapping is accurate
*within* held-out paired compounds (about `.45` direct-pIC50 MAE with measured
primary) but extrapolates into a region where almost no direct labels exist.
The model moves toward weak values and improves truth-defined weak-label MAE,
while harming the overall official loss on both split families. This rejects
these proxy recipes for candidate selection. It does not measure their
performance on the unlabelled cohort or on the blind test.

This selection pattern was already documented independently in the
[public-label availability audit](2026-09-11-label-selection.md): the area
at primary log2FC ≥−0.5 contains 2,823 compounds and no public D6 direct
labels. Missing direct values must not be asserted to be numerical potency
measurements.

## 3. Alternative primary predictors move the bottleneck only slightly

Run `regression/20260926T031827148427Z-1f7d39ab` holds the 647 frozen
features and outer-fold membership fixed. LightGBM, CatBoost and ExtraTrees
predict the measured D6 primary response from 4,376 labelled molecules.
Scores cover 4,313 nonanchor primary labels per split, including 1,457 with
direct D6 labels. The direct-pIC50 association is a screening readout; no
new downstream direct model was fitted from these primary predictions.

| Split | Primary model | Primary RMSE | Primary Spearman | Predicted primary vs direct Spearman |
| --- | --- | ---: | ---: | ---: |
| Parent | LightGBM | .8255 | .5851 | −.4233 |
| Parent | CatBoost | .8217 | .5908 | −.4271 |
| Parent | ExtraTrees | **.8210** | **.5951** | **−.4355** |
| Scaffold | LightGBM | .8217 | .5872 | −.4285 |
| Scaffold | CatBoost | .8206 | .5908 | −.4358 |
| Scaffold | ExtraTrees | **.8148** | **.5988** | **−.4498** |

ExtraTrees improves primary RMSE on both splits, but its direct-label ordering
gain is modest relative to the measured-response association in experiment 1.
This is not strong evidence for an expensive nested-feature and direct-regression
refit or a new bank member.

## Decision

Retain the submitted D6 vector. The internal primary assay contains substantial
D6 potency information, but simple proxy-label transfer fails badly under the
public D6 direct-label selection mechanism, and changing the primary model
family only slightly improves the fold-held signal. This combination explains
why expanding the model bank or treating missing labels as weak potency is
unlikely to close the current D6 gap. A future candidate needs a method that
handles the selection boundary explicitly and can be evaluated on relevant
held-out outcomes; this campaign does not supply such evidence.

Compact metrics and hashes are in
[`artifacts/tracked/regression-d6-primary-internal-20260926/`](../../artifacts/tracked/regression-d6-primary-internal-20260926/README.md).
Keyed OOF predictions remain in their recorded runs.
