# CYP2D6 weak-activity interventions on two regression splits

Four recorded runs test whether the weak-activity observation in
[the population diagnostic](2026-09-26-regression-d6-oof-population-diagnostic.md)
can improve held-out regression predictions. Parent-retained series is the
primary split; parent-retained scaffold is the prespecified sensitivity split.
All supervised fits exclude outer-validation compounds and retain the original
anchors in training. No hidden labels or interim leaderboard numbers enter any
fit or transformation. These are exploratory experiments on previously used
folds, not independent confirmation. OOF scores compare unplaced model
predictions; the submitted D6 test column also has a leaderboard-informed
affine placement, which is outside this experiment.

## Training-weight comparison

Run `regression/20260926T013939890389Z-249b950a` fits identical score47
CatBoost MAE readouts under five fixed training weights. The ordinary readout
is the within-experiment control. Weighting uses only outer-training labels:
either point pIC50 <4 (128 OOF cases) or confidence upper bound <4 (77 OOF
cases), at 2x or 4x weight.

| Recipe | Parent ST-RAE | Scaffold ST-RAE | Parent weak MAE | Scaffold weak MAE |
| --- | ---: | ---: | ---: | ---: |
| Uniform | **.8600** | .8568 | 1.918 | 1.901 |
| Point <4, 2x | .8674 | **.8510** | 1.854 | 1.852 |
| Point <4, 4x | .8714 | .8636 | 1.870 | 1.857 |
| Upper <4, 2x | .8704 | .8559 | 1.891 | 1.862 |
| Upper <4, 4x | .8697 | .8645 | 1.854 | 1.849 |

The 2x point-label recipe lowers weak-activity MAE, but worsens primary
ST-RAE by `.0073`. Every weighted recipe worsens parent ST-RAE. Scaffold has
one favorable point result, without agreement across splits. The simple
score47 control is also weaker than the submitted ensemble's parent OOF
`.8284`; no weighted model warrants a candidate.

## Can frozen features identify weak compounds?

Run `regression/20260926T014119877434Z-bf53d48d` fits outer-fold CatBoost
classifiers using score47 alone or score47 plus frozen medium600. The positive
label is either point pIC50 <4 or confidence upper bound <4. Both feature views
are independent of held-out challenge labels. The stronger medium view has:

| Positive definition | Prevalence | Parent AUROC / AP | Scaffold AUROC / AP | Top-decile positive rate, parent / scaffold |
| --- | ---: | ---: | ---: | ---: |
| Point <4 | 8.8% | .773 / .230 | .766 / .220 | 22.6% / 26.0% |
| Upper <4 | 5.3% | .836 / .192 | .845 / .182 | 19.2% / 19.9% |

The classifier can enrich weak activity by roughly 2.6–3.8x in its highest
probability decile. AP remains modest because weak labels are rare. This is
useful molecular signal, not evidence that a direct regression adjustment
will reduce ST-RAE.

## Fixed probability adjustment

Run `regression/20260926T014426087770Z-3eb82524` combines each classifier's
fold-held probability with the matched Jeremy anchor-all LightGBM OOF on both
splits. Within each outer-training fold, it computes weak prevalence and the
mean pIC50 gap between weak and other labels. The fixed adjustment is
`strength × (probability − train prevalence) × train mean gap`, with strengths
`.25`, `.5` and `1`. The fold-held label never enters those constants.

| Model | Parent ST-RAE | Scaffold ST-RAE |
| --- | ---: | ---: |
| Jeremy anchor-all control | **.8328** | **.8222** |
| Point classifier, 0.25 | .8365 | .8264 |
| Upper classifier, 0.25 | .8360 | .8251 |
| Point classifier, 0.5 | .8469 | .8370 |
| Upper classifier, 0.5 | .8429 | .8322 |

Strength 1 is worse still: `.8866/.8755` with point labels and
`.8660/.8573` with upper-bound labels. For point-label strength `.25`, weak
MAE falls slightly (`1.926→1.915` parent; `1.924→1.915` scaffold), while
overall ST-RAE and MAE worsen. The correction changes ranking only slightly;
parent Spearman rises `.4494→.4536` for the point recipe, but the ST-RAE
penalty dominates the primary objective.

## Nested classifier feature

Run `regression/20260926T014651429774Z-f3a4d62c` also trains a model to
combine the weak-classification probability with continuous potency. Within
each outer-training fold, a three-way molecular-group cross-fit creates the
classifier feature for **every training row**; a classifier trained on the
complete outer-training fold predicts the outer validation rows. Matched
CatBoost MAE readouts receive the same frozen score47+medium600 features and
training labels, with or without the one probability feature.

| Readout | Parent ST-RAE | Scaffold ST-RAE | Parent Spearman | Scaffold Spearman |
| --- | ---: | ---: | ---: | ---: |
| Matched control | **.8239** | **.8190** | **.4828** | **.4838** |
| Nested weak feature | .8346 | .8286 | .4582 | .4725 |

The added feature worsens all four headline values and raises parent
weak-activity MAE `1.884→1.898`; scaffold weak MAE falls slightly
`1.894→1.881`. Parent ST-RAE worsens in four of five folds, scaffold in all
five. The probability feature's training values are leakage-safe, so this
negative result cannot be dismissed as same-row target leakage.

## Decision

Discard weak-label reweighting, fixed probability correction, and this nested
classifier feature as deployment directions. Preserve the classifier as
evidence that weak activity is partly identifiable from molecular features.
The current submitted D6 vector remains the control. These experiments give
no basis for a candidate or submission. Compact metrics are preserved in
[`artifacts/tracked/regression-d6-weak-activity-20260926/`](../../artifacts/tracked/regression-d6-weak-activity-20260926/README.md).
