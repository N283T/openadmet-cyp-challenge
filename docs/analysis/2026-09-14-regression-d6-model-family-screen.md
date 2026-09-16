# CYP2D6 model-family diversity screen

## Purpose and design

Run `20260914T121807545564Z-5d24aefa` tests whether simple or structurally
different readouts supply complementary CYP2D6 molecular ordering. It is a broad
hypothesis-generation screen, not a tuning contest. All models use the retained
parent CV `20260911T104703655337Z-a3923255`; anchors remain in training and all
1,457 nonanchor CYP2D6 labels receive exactly one OOF prediction per method.

Three fixed feature views are compared:

- `scores51`: upstream public predictor scores47 plus new nested log2fc4;
- `medium604`: frozen upstream medium embedding600 plus new nested log2fc4;
- `hybrid2699`: medium600, scores47, Morgan2048 and new nested log2fc4.

Each view is read by fixed Ridge, Elastic Net, PLS16, ExtraTrees, Random Forest,
RBF-SVR, Nystroem-Ridge, kNN, MLP, XGBoost, CatBoost and compact KAN recipes.
The KAN is an explicit two-layer triangular-spline implementation on
outer-training-fitted PCA64, not the `pykan` package. Models are assessed alone
and at prespecified 10%, 25% and 50% additions to latest development OOF
`ae562a59`; moment matching uses prediction moments only and does not change
standalone ranks.

## Main result

The latest development reference scores CYP2D6 ST-RAE **0.822681**, Spearman
**0.476425**. The older simplified two-component reference scores 0.828397.

| Addition to latest | ST-RAE | Delta | Spearman | Delta |
|---|---:|---:|---:|---:|
| Hybrid CatBoost, 50%, moment matched | **0.817687** | **-0.004994** | 0.479956 | +0.003531 |
| Scores Random Forest, 25%, moment matched | 0.820409 | -0.002272 | 0.479560 | +0.003136 |
| Scores PLS16, 25%, moment matched | 0.821607 | -0.001074 | **0.481714** | **+0.005289** |
| Medium CatBoost, 25%, moment matched | 0.822096 | -0.000584 | 0.480545 | +0.004121 |

Hybrid CatBoost is also the best standalone fixed recipe: moment-matched
ST-RAE 0.823526 and Spearman 0.471292. It nearly matches the latest ensemble
alone and improves it at 25% and 50%. The improvement is not a low-correlation
new residual regime: prediction correlation with latest is 0.9404 and residual
correlation 0.9903. It is a modest but coherent tree-bias complement.

The simple 51-feature models are useful diversity probes. Random Forest improves
the blend with lower prediction correlation 0.8834. PLS produces the best tested
blend Spearman despite weaker standalone ST-RAE, illustrating why the screen
should not rank models only by standalone loss.

The high-dimensional linear recipes fail badly without stronger dimension or
regularization control. KAN also fails in the tested compact configuration:
moment-matched standalone ST-RAE is 0.966-0.977 with Spearman 0.281-0.310 across
views, and even 10% additions do not beat the reference. This rejects the fixed
PCA64/two-layer/triangular-grid recipe, not KANs generally.

## Independent audit and uncertainty

Completed audit `20260914T122524849329Z-edf03f72` verifies every method's exact
keys, labels, confidence bounds, fold coverage and row count. The first audit
attempt `20260914T122457344446Z-284cd352` failed because the baseline artifact
does not contain molecular groups; it is retained and the corrected audit uses
the groups already saved in the screen OOF.

Conditional 5,000-draw molecular-group bootstrap intervals for ST-RAE delta are:

| Candidate | 95% interval | P(improvement) | Improving folds |
|---|---:|---:|---:|
| Hybrid CatBoost 50% | [-0.01150, +0.00097] | 94.6% | 4/5 |
| Scores Random Forest 25% | [-0.00687, +0.00200] | 84.6% | 3/5 |
| Scores PLS16 25% | [-0.00555, +0.00341] | 67.9% | 3/5 |
| Medium CatBoost 25% | [-0.00409, +0.00277] | 64.0% | 3/5 |

These intervals are conditional on saved predictions and post-screen selection;
they are not multiple-comparison or selection adjusted. Hybrid CatBoost worsens
fold4 by 0.00741 and its interval crosses zero, so it is a follow-up candidate,
not an adopted replacement.

## Decision

Retain hybrid CatBoost and scores-only Random Forest/PLS as useful new D6 model-bank
entries. Confirm CatBoost across a small prespecified depth/loss/seed family and
use fold-excluded blend-weight selection before test inference. Do not promote
the current KAN recipe; any KAN follow-up should change the formulation materially,
for example an additive one-layer spline model on the compact 51-feature view,
rather than tuning this failed architecture narrowly. No test inference,
candidate, or submission was created.

## Pairwise-ranking follow-up

Completed run `20260914T123940860032Z-de101db4` tests whether learning relative
order directly supplies a different D6 signal. It keeps the same parent CV,
three feature views and 1,457-row OOF universe. Pairwise logistic regression,
linear RankSVM, and RankNet use outer-training-fitted standardization and PCA64;
their 30,000 sampled training pairs include only compounds whose confidence
intervals do not overlap. CatBoost PairLogit uses those same explicit pairs on
the original feature view. XGBoost's `rank:pairwise` objective instead orders
the complete outer-training group by its continuous label. Every learned
utility is mapped to pIC50 by a positive affine calibration fitted only on the
corresponding outer-training labels, so ranking is unchanged.

The best fixed ST-RAE result is a 10% addition of scores51 XGBoost pairwise:

| Addition to latest | ST-RAE | Delta | Spearman | Delta |
|---|---:|---:|---:|---:|
| Scores XGBoost pairwise, 10% | **0.820434** | **-0.002247** | 0.477486 | +0.001061 |
| Scores pairwise logistic, 10% | 0.821434 | -0.001247 | 0.479816 | +0.003392 |
| Hybrid CatBoost PairLogit, 10% | 0.821544 | -0.001137 | 0.479045 | +0.002620 |
| Scores RankSVM, 10% | 0.821642 | -0.001039 | 0.479857 | +0.003432 |
| Scores pairwise logistic, 25% | 0.823008 | +0.000327 | **0.481861** | **+0.005437** |

The result is useful mainly as a low-weight complement. Scores51 XGBoost has
prediction correlation 0.6738 with the latest reference, materially below the
other selected rankers, although residual correlation remains 0.9567. Direct
pairwise logistic and RankSVM give the larger rank-correlation gains. At 25%,
however, they begin to trade ST-RAE for Spearman. Standalone pairwise models are
generally weaker than the existing regression ensemble. RankNet is particularly
weak (scores51 standalone ST-RAE 1.364 and Spearman 0.292), so its apparent
residual diversity is error rather than a usable complement.

Audit `20260914T124500580752Z-2c999b86` verifies exact keys, labels, confidence
bounds and all 15 view-method combinations. Conditional 5,000-draw
molecular-group bootstrap results are:

| Candidate | 95% ST-RAE-delta interval | P(improvement) | Improving folds |
|---|---:|---:|---:|
| Scores XGBoost pairwise 10% | [-0.00534, +0.00069] | 93.1% | 4/5 |
| Scores pairwise logistic 10% | [-0.00387, +0.00137] | 82.1% | 4/5 |
| Hybrid CatBoost PairLogit 10% | [-0.00359, +0.00147] | 81.4% | 4/5 |
| Scores pairwise logistic 25% | [-0.00615, +0.00703] | 46.5% | 3/5 |

These intervals are post-screen and not selection adjusted. The XGBoost 10%
result is the strongest pairwise candidate but still crosses zero and worsens
fold4 by 0.00363. Retain the scores51 pairwise models as new model-bank entries;
next compare a small prespecified pair-margin/weight family and evaluate them
with fold-excluded blend selection. Do not promote RankNet or create inference,
a candidate, or a submission. Failed runs `20260914T123408289115Z-3b38f049`,
`20260914T123443754981Z-114af9d9`, and
`20260914T123522253242Z-f699e207` record, respectively, implicit CatBoost pair
explosion, an incremental Parquet dtype mismatch, and a final nested-record CSV
serialization error; all were corrected before the completed run.

## Ranking, kernel, and classical follow-up

Completed run `20260914T125646427763Z-1707a0e2` tests ten additional fixed
recipes under the same parent CV. The scores51 recipes are interval-plus-pairwise
linear models at pair-loss strengths 0.25 and 1.0, confidence-gap-weighted
pairwise logistic regression, five- and eight-bin ordinal logistic regression,
an additive spline GAM, and Huber regression. Three exact Tanimoto kernel-ridge
recipes use Morgan2048 with alpha 0.1, 1, and 10. Pair sampling, scaling,
ordinal thresholds, and any affine calibration are fitted within outer training.

| Addition to latest | ST-RAE | Delta | Spearman | Delta |
|---|---:|---:|---:|---:|
| Scores interval+pairwise(0.25), 25% | **0.819885** | **-0.002795** | 0.480925 | +0.004501 |
| Scores Huber, 25% | 0.820389 | -0.002292 | **0.480976** | +0.004551 |
| Scores ordinal-logit5, 10% | 0.820686 | -0.001995 | 0.479463 | +0.003039 |
| Scores weighted-pairwise logistic, 10% | 0.821202 | -0.001479 | 0.480167 | +0.003742 |
| Morgan Tanimoto KRR alpha0.1, 10% | 0.824677 | +0.001997 | 0.481026 | +0.004602 |

Audit `20260914T125740333957Z-ff207c2f` verifies complete exact-key, label,
confidence-bound, and method coverage. Its conditional 5,000-draw molecular-
group bootstrap gives:

| Candidate | 95% ST-RAE-delta interval | P(improvement) | Improving folds |
|---|---:|---:|---:|
| Scores ordinal-logit5 10% | [-0.00470, +0.00070] | **92.1%** | **4/5** |
| Scores weighted-pairwise logistic 10% | [-0.00406, +0.00109] | 86.6% | 3/5 |
| Scores Huber 25% | [-0.00709, +0.00235] | 83.5% | 4/5 |
| Scores interval+pairwise(0.25) 25% | [-0.00836, +0.00288] | 83.1% | 3/5 |
| Morgan Tanimoto KRR alpha0.1 10% | [-0.00190, +0.00604] | 15.5% | 1/5 |

The interval-plus-pairwise model has the best pooled point estimate, but its
gain is concentrated in folds0-2 and it worsens fold4 by 0.00970. Ordinal-logit5
has the strongest stability: similar uncertainty to the earlier scores51
XGBoost-pairwise 10% result, four improving folds, and only a 0.00215 fold4
regression. Huber provides a useful classical corroboration but has residual
correlation 0.9783. Tanimoto alpha0.1 is genuinely diverse (prediction/residual
correlations 0.6371/0.9150) but converts that diversity into worse ST-RAE in
four folds; do not retain it solely for its Spearman gain.

Retain ordinal-logit5 10% as the stable alternative to XGBoost-pairwise 10%,
and retain interval-plus-pairwise 25% only as the higher-variance alternative.
Any combined blend of the two must be assessed as a new prespecified recipe,
not inferred from these post-screen results. Reject this Tanimoto KRR family;
do not create test inference, a candidate, or a submission.

## Generic KERMT and MoLFormer follow-up

Runs `20260914T135108660328Z-387c0bcc` and
`20260914T135140332924Z-1d133b02` extract frozen generic KERMT-v2 and
MoLFormer-c3 representations for all 6,896 challenge compounds. These are
untouched released foundation checkpoints, not models continued on PXR or CYP
challenge labels. The finite feature matrices have 3,200 and 768 dimensions.

Screen `20260914T135408977818Z-a444a0b9` holds the parent folds fixed and uses
four fixed small-data readouts. Each encoder is tested alone and with upstream
scores47; nested log2fc4 is supplied in all cases.

| Addition to latest | ST-RAE | Delta | Spearman | Delta |
|---|---:|---:|---:|---:|
| MoLFormer ExtraTrees, 25%, moment matched | **0.819980** | **-0.002701** | 0.478349 | +0.001924 |
| Scores+MoLFormer ExtraTrees, 25%, moment matched | 0.820827 | -0.001854 | **0.480111** | **+0.003686** |
| Scores+KERMT ExtraTrees, 10%, raw | 0.822237 | -0.000444 | 0.477309 | +0.000885 |

MoLFormer also wins the standalone comparison: scores+MoLFormer ExtraTrees has
ST-RAE 0.848606 and Spearman 0.461811; MoLFormer-only ExtraTrees has
0.850825/0.441755. Corresponding KERMT results are 0.859267/0.449637 and
0.880525/0.421543. PCA64 linear/kernel and PLS recipes are substantially weaker.

Audit `20260914T135738169248Z-447140d1` gives:

| Candidate | 95% ST-RAE-delta interval | P(improvement) | Improving folds |
|---|---:|---:|---:|
| MoLFormer ExtraTrees 25% | [-0.00701, +0.00184] | 88.4% | 3/5 |
| Scores+MoLFormer ExtraTrees 25% | [-0.00584, +0.00206] | 82.4% | 4/5 |
| Scores+KERMT ExtraTrees 10% | [-0.00238, +0.00148] | 67.4% | 4/5 |

The generic MoLFormer feature is worth retaining, but the uncertainty interval
crosses zero and its residual correlation with latest is high (0.981-0.985).
Prioritize a prespecified small MoLFormer tree/readout family and combinations
with the earlier CatBoost, ordinal, and pairwise signals. KERMT is lower priority
under the tested pooling/readouts. No inference, candidate, or submission was
created.

## MoLFormer pooling and candidate combinations

Frozen feature run `20260914T140517644424Z-de551457` concatenates the CLS token
with attention-mask mean pooling, producing a complete finite 6,896-by-1,536
matrix. Screen `20260914T140653173662Z-24d529d2` compares mean and CLS+mean
pooling, with and without scores47, using fixed ExtraTrees, Random Forest,
CatBoost-MAE, and XGBoost recipes. Nested log2fc4 remains in every view.

| Addition to latest | ST-RAE | Delta | Spearman | Delta |
|---|---:|---:|---:|---:|
| Scores+CLS+mean Random Forest, 25%, matched | **0.818407** | **-0.004274** | 0.482580 | +0.006155 |
| Scores+CLS+mean Random Forest, 25%, raw | 0.818499 | -0.004182 | **0.483145** | **+0.006720** |
| Scores+CLS+mean CatBoost-MAE, 25%, matched | 0.818850 | -0.003831 | 0.479291 | +0.002867 |
| Scores+CLS+mean ExtraTrees leaf1, 25%, raw | 0.818938 | -0.003743 | 0.482515 | +0.006090 |

Audit `20260914T142335850035Z-781cf3a4` finds that both Random Forest placements
improve all five folds. The matched recipe has a conditional 5,000-draw
molecular-group interval `[-0.00859,-0.00006]` and 97.6% improvement
probability; raw has `[-0.00922,+0.00085]` and 94.4%. The upper bound for the
matched recipe is only barely negative and selection is post-screen, so this is
stronger evidence than the CLS-only screen but not confirmatory evidence.

Combination run `20260914T142603890764Z-be6e82d1` uses fixed averages of the
retained model predictions and fold-excluded level-2 fits:

| Combination | ST-RAE | Delta | Spearman | Improving folds |
|---|---:|---:|---:|---:|
| MoLFormer-RF25 + hybrid-CatBoost50 | **0.816803** | **-0.005877** | **0.482720** | 4/5 |
| All four retained candidates | 0.817780 | -0.004901 | 0.481142 | 4/5 |
| MoLFormer-RF25 + pairwise10 | 0.818485 | -0.004196 | 0.481091 | 4/5 |
| MoLFormer-RF25 + ordinal10 | 0.818703 | -0.003978 | 0.481507 | 4/5 |

Fold-excluded learned stacking fails decisively: Ridge and positive linear
stacks score 0.841-0.849 and worsen every fold, despite slightly higher
Spearman. Retain the simple fixed MoLFormer/CatBoost average only as a
prespecified follow-up candidate. Do not use learned stacking, and do not infer
or submit from this post-screen comparison.

## Cross-endpoint transfer

The D6-selected scores47+MoLFormer CLS+mean+log2fc4 feature view and three fixed
tree recipes were transferred without retuning to the other endpoints:

| Endpoint | Best fixed addition | Baseline ST-RAE | New ST-RAE | Baseline Spearman | New Spearman |
|---|---|---:|---:|---:|---:|
| 1A2 | ExtraTrees leaf1, matched 10% | 0.670809 | 0.670731 | 0.607416 | 0.607975 |
| 2C9 | Random Forest, matched 10% | 0.452634 | 0.452551 | 0.741016 | 0.741504 |
| 3A4 | ExtraTrees leaf1, matched 10% | 0.401599 | 0.401799 | 0.827518 | 0.827590 |

Runs are `20260914T143359558742Z-ec3b4fe8` (1A2),
`20260914T143751083136Z-59e6b2bd` (2C9), and
`20260914T144119889941Z-a2ea691e` (3A4). The tiny 1A2 and 2C9 changes are not
material, and every tested 3A4 addition worsens ST-RAE. This supports an
endpoint-specific D6 use rather than a global regression feature change. The
non-D6 rows were not promoted to bootstrap audit or inference.

## D6 leaderboard deployment

Full-fit run `20260914T150107799119Z-63fbf1ae` deploys the more stable standalone
MoLFormer Random Forest recipe, not the higher-gain post-screen CatBoost
combination. The model uses all 1,493 D6 training rows and 1,587 features:
scores47, frozen MoLFormer CLS+masked-mean1536, and full-fit log2fc4. Its raw
test predictions are moment-matched to the currently scored D6 vector before a
25% blend. The 1A2, 2C9, and 3A4 predictions remain bitwise unchanged.

The candidate is `20260914T150138241149Z-448b516b`; prepared submission
`20260914T150145548210Z-1ab27b1d` has SHA-256
`1fd13bfcb9f772689e1a973aa9133d4dcfe7c75954a567568530d8945ffb80db`.
Attempt `20260914T150221103444Z-46a8a9dd` was accepted as `sent` and is awaiting
leaderboard processing. No leaderboard conclusion is made from submission
acceptance alone.
Immediate D6 snapshot `20260914T150423637203Z-be8bb695` still contains the
prior 02:48 UTC N283T score (`.565520`), confirming that the new entry had not
yet propagated at that point.

Later snapshots `20260914T161853640642Z-152cf75a` (D6) and
`20260914T161857920931Z-2887788e` (overall) confirm the 15:02 UTC entry. The
blind result rejects the candidate: D6 ST-RAE worsens `.565520→.572470`,
Spearman worsens `.431248→.415591`, and D6 rank moves `21→22`. Overall
MA-ST-RAE worsens `.475974→.477712`, with overall rank unchanged at 22. This
overrides the favorable post-screen OOF evidence for deployment decisions.
