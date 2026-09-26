# Classification positive-call rate robustness

## Question and protocol

The accepted classification vector calls 49/750 blind molecules positive for
CYP2D6 and 249/750 for CYP3A4. Neither count is imposed by the challenge.
Because the submission calls all blind molecules but the scorer uses only
confidently labelled molecules, rate evaluation must rank every molecule in an
outer fold before applying the supplied-label mask.

Run `classification/20260926T144420289486Z-14270294` sweeps eight
prespecified D6 counts on the previously reconstructed September 15 control
analog. It uses all 6,146 held-fold molecules on UMAP and scaffold, with 1,497
supplied D6 labels. The UMAP and scaffold analogs are separate refits; neither
is an exact rerun of historical cross-fitted recipe selection. For A4, the
submitted equal-seven ensemble has only labelled-molecule OOF. Run
`classification/20260926T144248722886Z-fe484e91` therefore fits a fixed
Jeremy/CheMeleon LightGBM **proxy** and predicts all 6,146 molecules on each
split (3,584 supplied A4 labels). Run
`classification/20260926T144402127052Z-20212326` sweeps seven A4 counts on
that proxy. A4 proxy behavior does not establish the submitted ensemble's
best rate.

Each rate is applied to the full held-fold pool. Projected MCC uses only
supplied labels and a prespecified scored-population prior. The paired
95% intervals resample molecular groups with 2,000 replicates; they are
conditional on these OOF predictions and do not measure blind-domain error.
The inner choice for each held fold maximizes projected MCC on the other four
folds, then its chosen rate is evaluated on the held fold. The complete
[outputs](../../artifacts/tracked/classification-rate-robustness-20260926/README.md)
include the curves and exact fold choices.

## CYP2D6 control analog

At the working 7.06% scored-population prior:

| Count per 750 | UMAP delta MCC vs 49 | Scaffold delta MCC vs 49 |
|---:|---:|---:|
| 25 | −.0087 | −.0211 |
| 35 | −.0069 | −.0176 |
| 49 | reference | reference |
| 65 | −.0051 | +.0084 |
| 80 | −.0221 | +.0195 |
| 110 | −.0324 | +.0290 |

The pooled grid maximum is 49 on UMAP and 110 on scaffold, unchanged for
5%, 7.06% and 10% prior scenarios. On UMAP, the fold-excluded selector chooses
`49,25,49,49,35`; its combined MCC delta from fixed 49 is `−.0219`, with
paired interval `[−.0505,+.0040]`. On scaffold it chooses 110 in all five
folds and has delta `+.0290`, interval `[−.0071,+.0661]`. The split disagreement
is larger than the apparent gain from a rate change. No tested alternative
improves both split families, and every fixed-rate D6 comparison interval
includes zero. The prior sensitivity does not resolve this.

## CYP3A4 static-feature proxy

At the working 29.16% scored-population prior:

| Count per 750 | UMAP delta MCC vs 249 | Scaffold delta MCC vs 249 |
|---:|---:|---:|
| 150 | −.0383 | −.0422 |
| 200 | +.0057 | −.0075 |
| 249 | reference | reference |
| 300 | +.0002 | +.0053 |
| 350 | +.0060 | +.0169 |
| 400 | +.0024 | +.0032 |
| 450 | −.0063 | −.0089 |

The 150-call proxy is worse on both splits, with paired intervals below zero.
The apparent 350-call gain is small and its intervals cross zero on both
splits. The fold-excluded selector chooses 350 on every scaffold fold, but
UMAP choices are `350,400,200,200,350`. The selected policy changes pooled
MCC by `+.0169` on scaffold (`[−.0104,+.0439]`) and `−.0352` on UMAP
(`[−.0655,−.0034]`). Across 20%, 29.16% and 40% prior scenarios, selection
remains split dependent. The proxy's UMAP AUROC/AP `.7496/.4450` also trails
the saved labelled-only current-bank ranking `.8071/.5127`; its rate curve
should not be transferred to the submitted equal-seven score.

## Decision

Retain the accepted 49 D6 and 249 A4 calls as external reference decisions.
These experiments do not show that either count is optimal, but they provide
no supported replacement count. A leaderboard upload of a count-only variant
would test a split-dependent hypothesis without a convincing local gain, so
no candidate or submission is created from this audit. The next useful A4
rate test requires complete-pool OOF for the **actual** equal-seven recipe;
D6 needs a validation design whose model/rate contrasts better explain the
historical leaderboard reversals before a new rate is selected. The scoring
mask and test-domain shift remain unobserved.

The first D6 run `d1af53a0` repeated fold-choice rows eight times because
the choice loop was inside the rate loop. Corrected run `de0a4f72` fixed that
output. Final run `14270294` uses the shared endpoint-aware audit script and
is the result cited here. The rate curves in the earlier runs were unaffected;
all runs remain preserved for provenance.
