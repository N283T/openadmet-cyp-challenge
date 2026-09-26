# CYP2D6 regression: weak-primary upper-bound candidate

## Objective and starting point

The goal was a submit-ready D6 direct-inhibition regression candidate, with
simple deployment even if the research model was more involved. The current
submitted D6 vector is the September 14 source-placement column embedded in
the September 24 endpoint-wise submission. Its parent-series OOF two-member
control is ST-RAE `.828397`, Spearman `.475322`. A matched Jeremy anchor
LightGBM scaffold control is `.822203`, Spearman `.468652`. The parent and
scaffold controls are different model compositions; compare **changes within
each split**, not the two absolute scores.

The [official assay description](https://openadmet.github.io/octant-cyp-inhib-blog-post/)
reports a primary screen in the TDI condition before selected compounds
received direct and TDI dose-response measurements. Public D6 direct labels
omit the weak-primary-response region almost entirely. Earlier
attempts to turn measured D6 primary responses into exact pseudo-pIC50 labels
substantially worsened validation. This experiment instead treats a weak primary
response as a **one-sided upper-bound hypothesis** for direct pIC50. This is an
assumption for model development, not an experimentally verified direct value.

## Discarded model alternatives

Run `regression/20260926T040416620014Z-fe1960e3` added 27 explicit RDKit
physicochemical and amine-proxy features to the frozen score47+medium600
CatBoost readout. The frozen control exactly reproduced the prior matched run.
Adding physchem changed D6 ST-RAE `.82393→.82477` on parent and
`.81896→.82168` on scaffold, so it was discarded. Physchem alone was much
weaker (`.96983/.95819`).

Runs `regression/20260926T040811652764Z-d970f7d9` and
`regression/20260926T041050133113Z-a4d97c86` tested a shared neural trunk
with a D6-specific output and supervision from the other three direct endpoints.
The wide 647-input model overfit badly: D6 OOF ST-RAE was at least `1.35` and
Spearman near zero on both splits. A PCA64 plus score47 input, smaller network,
stronger dropout and shorter training improved the loss but not enough: its
best parent/scaffold ST-RAE was `1.092/1.125`, again with near-zero ranking.
Neither neural recipe enters the candidate.

## One-sided upper-bound training

Run `regression/20260926T041418410211Z-076c47f6` fits a 300-tree, 15-leaf
LightGBM on score47+medium600. Direct measured labels receive squared loss.
Additional compounds whose target primary log2FC is at least `−0.5` receive
only a penalty when predicted direct pIC50 exceeds `4.0`; predictions below
that upper bound receive no pseudo-label penalty. Outer-held molecules supply
neither direct labels nor primary constraints to that fold's fit. The fixed
constraint-weight screen was `0/.05/.1/.25`.

To test whether the bound can fill a selection gap, the run also hides weak-
primary direct labels on 1A2, 2C9 and 3A4, then scores those measured held-out
labels. At weight `.25`, weak-screen mean soft error improves parent/scaffold
from `1.005/1.018→.392/.383` on 1A2 and `.243/.259→.078/.081` on 2C9.
3A4 changes little (`.263/.285→.257/.265`). This is synthetic validation of
the **method**, not a measurement of the missing weak-D6 outcomes.

On D6's labelled OOF population, the bounded model itself trades scoring for
ordering. At weight `.25`, parent ST-RAE worsens `.918→.993` but Spearman
improves `.441→.470`; scaffold ST-RAE worsens `.937→.981` while Spearman
improves `.447→.479`. The model predicts lower pIC50 on average. It is retained
only as a complementary member, not a standalone replacement.

## Fixed blend and group audit

Run `regression/20260926T041934073362Z-5e9e8cb8` screens fixed additions to
the existing parent and scaffold controls. For each held-out fold, the new
member's predicted mean and SD are matched to the control's predicted mean and
SD, using predictions only. The selected recipe is the `.25` upper-bound model
at equal weight with the control. Full-fit deployment performs the analogous
moment match on the 750 blind predictions against the current submitted D6
column; it does not use hidden D6 outcomes.

| Split | Control ST-RAE | Candidate ST-RAE | Control Spearman | Candidate Spearman | Improving folds |
| --- | ---: | ---: | ---: | ---: | ---: |
| Parent series | .828397 | **.822346** | .475322 | **.483172** | 4/5 |
| Scaffold | .822203 | **.814045** | .468652 | **.486475** | 4/5 |

Recorded audit `regression/20260926T042057198235Z-174084b3` verifies exact
compound, endpoint, fold, truth, confidence-bound and molecular-group alignment.
Its 3,000-draw paired group bootstrap gives the candidate-minus-control
ST-RAE 95% interval `[−.01496,+.00271]` on parent and
`[−.01675,+.00058]` on scaffold. Conditional improvement fractions are
91.6% and 96.6%. Both intervals cross zero. This screen was selected after
examining fixed weights and placements on the same OOF, so the intervals are
**not** selection-adjusted confirmation. The gain is modest and uncertain,
although its direction is consistent across two split geometries and eight of
ten held-out folds.

## Full-fit candidate and deployment check

Run `regression/20260926T042258195328Z-2027e07c` refits the chosen model on
all 1,493 measured D6 direct labels plus 2,823 weak-primary upper-bound
examples. The 750 blind D6 raw predictions have mean `4.363` and SD `.545`.
After moment matching and equal blending with the current submitted D6 column,
the candidate has mean `3.570`, SD `.884` (current SD `.899`). The candidate
and current D6 vectors have Spearman `.956`; their mean absolute difference is
`.125` pIC50, with maximum `.791`. The 1A2, 2C9 and 3A4 columns are exactly
unchanged, and all 750 official names and SMILES pass validation.

Candidate `20260926T042332999625Z-c6480469` was registered and prepared as
submission `20260926T042337084384Z-c3578831`, SHA-256
`adc7c5a5b600a39dd5d3f6353842e3a1ac41cfc00fc5e27aeb322e78888b08b4`.
The submission dry run reports 750 valid rows and **sent=false**. No live upload
was made. Replaying the full-fit run as
`regression/20260926T042406809391Z-6b28df55` reproduces both the candidate
CSV and 534-KB LightGBM model file byte-for-byte.

## Decision

Keep the candidate as a deliberately exploratory D6-only update. It is a
reasonable prospective submission because it uses a new selection-aware loss,
improves both established OOF geometries, has an interpretable simple deployment
recipe, and changes the blind vector moderately. The D6 weak region lacks
measured direct labels, and prior small OOF gains have reversed on the blind
leaderboard. The confidence intervals also include no improvement, so do not
describe this as a proven upgrade before a prospective leaderboard result.
The current submitted vector remains the scored control until then.

Compact evidence, the fitted model, and the candidate payload are in
[`artifacts/tracked/regression-d6-primary-bound-20260926/`](../../artifacts/tracked/regression-d6-primary-bound-20260926/README.md).
