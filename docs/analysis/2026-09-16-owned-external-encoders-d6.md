# Project-owned external encoders for CYP2D6

## Scope

This experiment replaces third-party encoder dependence with random-initialized,
project-owned message-passing encoders. Evaluation is restricted to
`CYP2D6_pIC50_direct_inhibition`; the other regression targets and their banks are
unchanged. Every readout is an independent atomic prediction. No blend weights,
fractional additions, or nested members were searched.

## Data and leakage control

Pinned public OpenADMET sources were fetched in run
`20260916T072328099350Z-fa03aef6`. Curation run
`20260916T072819161332Z-a374a545` standardizes structures and excludes every
challenge connectivity before training. The resulting corpus contains 53,582
structures and 221,471 observations, including 10,086 structures and 45,135
observations newly added from ExpansionRx, PXR and the non-overlapping part of
Octant.

The overlap audit is important: 6,779 candidate observations were rejected.
Nearly the entire 1,340-row Octant inhibition set and all usable Octant reactivity
rows overlap the challenge structure universe, so they are not used as external
training labels. The AstraZeneca 120-compound SoM data are atom-level labels and
were not incorrectly converted into graph-level regression values.

## ATTNSOM representation control

Runs `2a529c61`, `df779279`, `65b4b896` and `a45372b5` fit LightGBM readouts on
the existing project-trained ATTNSOM embedding. Matched audit `ac618b17` gives:

| representation | D6 ST-RAE |
|---|---:|
| raw512 | 0.914296 |
| PCA64 | 0.939371 |
| PCA128 | 0.942972 |
| PCA256 | 0.957452 |

Raw is consistently best, but all variants are too weak for the D6 bank.

## Broad external encoder grid

Random-initialized width-256 and width-600 encoders were trained in runs
`c05c3377` and `2dda20d9` on the expanded assay corpus with the same fixed
validation connectivity split and seed 42. Audit `87e5e4d5` covers raw/PCA
LightGBM readouts; audit `b875dc14` covers selected TabPFN readouts.

| encoder/readout | D6 ST-RAE |
|---|---:|
| width600 raw, TabPFN | 0.923624 |
| width600 raw, LightGBM | 0.939020 |
| width256 raw, TabPFN | 0.947923 |
| width600 PCA128, LightGBM | 0.951189 |
| width600 PCA64, LightGBM | 0.953007 |
| width600 PCA128, TabPFN | 0.964378 |
| width600 PCA256, TabPFN | 0.970082 |
| width256 raw, LightGBM | 0.975567 |

The larger raw representation is better than the smaller representation and all
compressed versions, but remains far outside the current bank range.

## Direct-CYP2D6 focused variant

Subset run `c23f09a3` retains 31,438 leak-excluded external CYP2D6 observations
over 28,220 compounds and 1,829 assays. Width-600 training run `ba663997` uses
only those labels. Matched audit `ce97015e` gives raw LightGBM ST-RAE 0.920890
and raw TabPFN 0.938399. Focusing the labels improves LightGBM relative to the
broad raw600 encoder but does not approach the current D6 ensemble.

## Decision

Do not add any of these predictions to the CYP2D6 bank. Do not expand training
seeds: the deficit is much larger than plausible seed variation. Preserve the
runs as informative negative results and as reusable project-owned encoders.
Future work should change supervision or architecture rather than tune ensemble
weights around these candidates. Raw embeddings remain the default first test;
PCA is still evaluated when warranted, but compression was uniformly harmful in
this batch.

## Correction: actual Jeremy-replacement comparison

The initial screen above was intentionally restricted to the active D6 bank, but
that does not answer the broader question of whether the project-owned encoder
can replace Jeremy's frozen upstream encoder. Audit
`20260916T075619930058Z-e73e6670` corrects the scope: it uses the identical four
regression targets, parent-retained folds, atomic embedding-only LGBM/TabPFN
readouts, and includes Jeremy raw600 as the direct control.

| representation/readout | Macro | 1A2 | 2C9 | 2D6 | 3A4 |
|---|---:|---:|---:|---:|---:|
| Jeremy raw600 TabPFN | 0.649557 | 0.727129 | 0.533550 | 0.881106 | 0.456444 |
| Jeremy raw600 LGBM | 0.651179 | 0.728514 | 0.535216 | 0.858244 | 0.482742 |
| Owned width600 raw TabPFN | 0.762724 | 0.811255 | 0.705242 | 0.923624 | 0.610775 |
| Owned width256 raw TabPFN | 0.770749 | 0.824602 | 0.697177 | 0.947923 | 0.613295 |
| Owned width600 raw LGBM | 0.792471 | 0.837450 | 0.750695 | 0.939020 | 0.642722 |
| Owned width600 PCA256 TabPFN | 0.794830 | 0.838642 | 0.738218 | 0.970082 | 0.632379 |
| Owned width600 PCA256 LGBM | 0.812992 | 0.850947 | 0.776744 | 0.958277 | 0.666000 |
| Owned width256 raw LGBM | 0.814095 | 0.855667 | 0.768747 | 0.975567 | 0.656401 |

The correct replacement conclusion is therefore negative: the current owned
encoder loses to Jeremy on every endpoint and by 0.1132 macro ST-RAE at best.
Keep Jeremy unchanged. The D6-only work remains a valid supervision ablation,
not the replacement decision.
