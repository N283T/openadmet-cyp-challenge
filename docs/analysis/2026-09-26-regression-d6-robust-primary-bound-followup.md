# CYP2D6 robust primary-bound follow-up

**Scope.** This is a CYP2D6 direct-inhibition regression experiment. Only the
D6 column changes in the prepared submission. The 1A2/2C9/3A4 masking
experiment below is a transfer stress test because D6 has no measured direct
labels in the weak-primary cohort. It is not D6 validation and is not part of
the D6 candidate model or its ensemble.

## Why this follow-up was necessary

The [first one-sided-bound candidate](2026-09-26-regression-d6-primary-bound-candidate.md)
improved D6 parent/scaffold ST-RAE only `.0061/.0082` against its matched controls.
Both molecular-group bootstrap intervals crossed zero. Its direct-label squared
loss produced a weak standalone model; adding primary-bound constraints worsened
standalone ST-RAE further. The selected blend and 50% weight were picked on the
same OOF used to quote their gain. This is a plausible exploratory candidate,
not a demonstrated upgrade. The blind D6 pIC50 distribution also differs from
the labelled population, and no public direct values exist among the 2,823
weak-primary D6 compounds used as one-sided examples.

## Robust-loss experiment

Run `regression/20260926T085732408866Z-0167cd1b` replaced the direct-label
squared loss with pseudo-Huber loss `sqrt(delta² + residual²)` and gave weak
primary compounds a one-sided pseudo-Huber penalty only above pIC50 4. The
model used the same frozen 47 scores plus 600 medium embeddings, 300-tree
15-leaf LightGBM, parent-retained and scaffold splits, and outer-fold exclusion
as the first bound experiment. The fixed grid was `delta={.2,.5}` and
`constraint weight={0,.15,.35}`. No D6 direct labels from held molecules enter
their model fit or pseudo-bound set.

At `delta=.5`, adding weight `.35` changes standalone D6 ST-RAE
`.9104→1.0100` on parent and `.9204→1.0081` on scaffold. Spearman rises
`.4417→.4743` and `.4533→.4888` respectively. The bound again supplies
ordering information while badly shifting absolute predictions. In particular,
the robust member alone is not a replacement for the current D6 column.

## Fixed blends and selection audit

Run `regression/20260926T085959697226Z-0dbf57d4` compared weights
`.1/.25/.5`, raw and per-fold prediction-moment placement, against the current
two-member parent control and matched Jeremy anchor scaffold control. The best
common recipe is robust delta `.5`, bound weight `.35`, fold-moment matching,
and 50% blend. Moment matching uses predicted means and SDs only, not held
labels.

| Split | Control ST-RAE | Robust blend | Control Spearman | Robust blend |
| --- | ---: | ---: | ---: | ---: |
| Parent | .828397 | **.819829** | .475322 | .484868 |
| Scaffold | .822203 | **.810041** | .468652 | .492689 |

Paired molecular-group audit `regression/20260926T090048786219Z-8a2ca577`
finds improvement in all five folds of each split. Its conditional 3,000-draw
95% intervals for candidate-minus-control ST-RAE are
`[-.01727,+.00028]` on parent and `[-.02111,-.00341]` on scaffold. The parent
interval still includes no improvement. The scaffold interval excludes zero,
but these intervals are **not corrected for recipe selection** on the same OOF.

To estimate that selection effect, run
`regression/20260926T090204591664Z-07a300e1` chooses the lowest-ST-RAE recipe
from the other four folds separately for each held fold. The same recipe is
selected in all five scaffold folds and four of five parent folds. Assembled
held-fold ST-RAE is `.823924` parent and `.810041` scaffold. Thus the selected
parent point estimate `.819829` overstates the result of this fold-excluded
selection exercise by about `.0041`. This is still not an independent test:
models, feature views, split choices and the robust-loss idea followed many
earlier experiments on the same labelled population.

Direct paired comparison `regression/20260926T090326347564Z-2afbbcef`
against the previous squared-loss bound candidate gives robust-minus-previous
ST-RAE `−.00252` parent and `−.00400` scaffold. Group-bootstrap 95% intervals
are `[-.00902,+.00436]` and `[-.00969,+.00190]`. The robust version is the
better point estimate but is **not conclusively better than the first candidate**.

## Other-CYP proxy test, not D6 validation

Run `regression/20260926T090111824366Z-5a9f046c` hides measured weak-primary
direct labels from training for 1A2, 2C9 and 3A4, then evaluates those values.
It uses the selected `delta=.5`, bound weight `.35` and a matched no-bound
robust control. Weak-screen mean soft error changes as follows:

| Target | N | Parent control → bound | Scaffold control → bound |
| --- | ---: | ---: | ---: |
| 1A2 | 249 | 1.008 → .343 | 1.004 → .337 |
| 2C9 | 410 | .253 → .072 | .279 → .079 |
| 3A4 | 522 | .271 → .257 | .278 → .259 |

Strong-primary selected-like error worsens on 1A2 and 2C9. This supports a
real weak-region tradeoff and shows that the constraint is not universally
helpful. Those targets have measured hidden weak-region labels; D6 does not.
Transfer to the blind D6 distribution is an unverified assumption. The D6
decision must therefore rest on the D6 parent/scaffold OOF, its selection
audit, and ultimately a prospective D6 leaderboard result. These proxy
numbers should not be counted as additional D6 validation folds.

## Full fit and decision

Run `regression/20260926T090408836539Z-239369e6` fits on all 1,493 D6 direct
labels plus 2,823 weak-primary one-sided examples. Its raw blind D6 predictions
have mean `4.351`, SD `.568`, whereas the currently submitted D6 vector has
mean `3.570`, SD `.899`. The deployment recipe matches the raw vector's mean
and SD to the current blind vector, then blends equally. This sizeable
placement is a **material calibration dependency**; OOF improvement alone
cannot verify it. The new 750-row D6 vector has SD `.883`, Spearman `.953`
against current, mean absolute change `.132`, maximum change `.660` pIC50.
The other three endpoints and molecule identities are unchanged.

Candidate `20260926T090542937263Z-099c328f` and prepared submission
`20260926T090548323924Z-46188564` passed the 750-row dry run with
`sent=false`; payload SHA-256 is
`399cd07c66c0d4e7e1fb468b6a36a471d1572951e7f10e324e0b597d216c7b29`.
Replay `regression/20260926T090601408220Z-5b3cc170` matches the CSV to
within `8e-15` pIC50; 94 rows differ only at floating-point rounding and
the model/CSV hashes are not byte-identical across threaded LightGBM runs.

The robust candidate improves the available **D6** OOF evidence and is ready
for a prospective D6-only leaderboard probe. It remains **uncertain** as an upgrade:
the independent weak-D6 labels are absent, the head-to-head gain over the
first candidate is small, the parent fold-excluded score regresses toward the
control, and test calibration relies on matching the old submitted vector.
Do not replace the scored control on OOF evidence alone. Compact evidence,
model weights and prepared payload are in
[`artifacts/tracked/regression-d6-primary-robust-20260926/`](../../artifacts/tracked/regression-d6-primary-robust-20260926/README.md).

## Prospective result (recorded after the original analysis)

The prepared candidate was subsequently submitted at the user's direction.
The [Live result](2026-09-26-regression-d6-robust-live-result.md) confirms D6
ST-RAE `.565520→.549341` and Spearman `.431249→.441151` against the previous
scored D6 column. This is evidence on the Live scored half, not a full-test
result or proof that the untested final-SD correction is beneficial.
