# TDI interim leaderboard by target

## Evidence

The frozen [September 26 leaderboard snapshot](../../artifacts/tracked/leaderboard-snapshot-20260926/manifest.json)
contains Live and Interim CSVs for TDI overall, CYP2D6 and CYP3A4. Recorded
run `shared/20260926T002255140309Z-2e7b1859` compares the boards and their
N283T rows. The Interim board evaluates the September 24 submission on the
full blinded test set once; the Live board evaluates one half. N283T's
submission timestamp is September 24 06:24 UTC on both, so these are the same
submitted predictions across scoring populations.

## Result

| TDI endpoint | N283T Live | N283T Interim | Interim leader |
|---|---:|---:|---:|
| Overall MA-MCC | rank 12, .4193 | **rank 5, .3798 (Tier 1)** | nova, .3980 |
| CYP2D6 MCC | rank 31, .3086 | rank 25, .2573 | nova, .3450 |
| CYP3A4 MCC | **rank 1, .5299** | **rank 4, .5023** | EleanorRigby, .5207 |

The organizer places the top 26 overall TDI entries in Tier 1. Their score
differences are not statistically resolved by that tier test. N283T is in this
group; rank 5 alone does not support a claim that it is better or worse than
another Tier 1 entry.

## Secondary metrics for N283T

| Interim board | Accuracy (rank) | Precision (rank) | Recall (rank) | F1 (rank) |
|---|---:|---:|---:|---:|
| CYP2D6 | .9128 (29) | .2878 (22) | .3224 (66) | .3018 (23) |
| CYP3A4 | .7828 (2) | .6393 (23) | .6848 (46) | .6608 (5) |
| Macro average | .8478 (6) | .4635 (18) | .5036 (55) | .4813 (5) |

These ranks compare the corresponding metric across all 114 Interim entries;
they are separate from the official MCC ranking. F1 broadly supports the same
picture as MCC: A4 is strong, D6 has more room. Accuracy needs care because
the D6 positive class is rare; a high accuracy can coexist with poor positive
class identification. Precision and recall are more informative about that
tradeoff, but neither alone is the challenge objective.

### CYP2D6

N283T's interim precision/recall are `.2878/.3224`, versus `.4725/.3006` for
nova. N283T's F1 is `.3018` versus nova's `.3641`, and accuracy is `.9128`
versus `.9388`. N283T is `.0877` MCC behind nova on D6. This comparison suggests a
different precision/recall operating point but does not reveal the identities
or labels of individual errors. Other high-scoring D6 entries make different
tradeoffs: jeremy has precision/recall `.3232/.3887`, while newcyp has
`.2554/.5052`. Aggregate leaderboard metrics do not identify a superior
threshold for N283T's fixed score vector.

### CYP3A4

N283T's interim precision/recall are `.6393/.6848`. The endpoint leader,
EleanorRigby, has `.5633/.8623`; the second-place CYPblend has
`.6673/.6591`. N283T's F1 is `.6608` (rank 5) and accuracy is `.7828`
(rank 2), versus EleanorRigby's `.6810` and `.7503`. N283T is only `.0185`
MCC behind the endpoint leader.
Individual Precision and Recall ranks are affected by where entrants place
their decision boundaries; neither alone means the current A4 operating point
is poor. Its top-five MCC and F1 support a fairly strong balance.

The Live top-10 overlap with the Interim top 10 is only 3, and rank Spearman
among the 114 common entrants is `.756`. A4's half-test ranking is particularly
unstable, so the Live rank 1 should not be treated as the full-test result.

### Is the fixed A4 positive count appropriate?

The exact submitted vector calls **249 of 750** blind compounds positive
(33.2%). This count was retained through the September 13 and September 15
model changes; those changes swapped the identity of positive calls at the
same count. An earlier 249-call A4 vector scored Live MCC `.4307`, while a
later 378-call vector scored `.3750`; their rankings also differed. A recorded
counterfactual decomposition attributed `.0278` of the subsequent macro-MCC
recovery to restoring the A4 call rate
([analysis](2026-09-13-classification-diverse-lb-result.md)). This establishes
that call rate matters, but does not prove 249 is the full-test optimum for the
current score vector. The interim table provides
one point on its Precision/Recall curve, not the curve itself. EleanorRigby's
higher recall and lower precision could reflect a different call rate, ranking,
or both; individual predictions and labels are unavailable.

Recorded run `classification/20260926T011101091386Z-618feba5` uses the
Interim Accuracy, Precision and Recall means to infer approximate positive
call fractions **within the scored subset**. The identity is
`prevalence ≈ (1 − accuracy) / (1 − 2 × recall + recall / precision)` and
`call fraction ≈ prevalence × recall / precision`. The compact
[derived table](../../artifacts/tracked/classification-a4-interim-operating-points-20260926/operating_points.csv)
preserves the estimates.

| Entrant | A4 MCC rank | Precision | Recall | Approx. scored positive-call fraction |
|---|---:|---:|---:|---:|
| N283T | 4 | .6393 | .6848 | 33.2% |
| EleanorRigby | 1 | .5633 | .8623 | 47.4% |
| CYPblend | 2 | .6673 | .6591 | 30.6% |
| nova | 24 | .6472 | .5772 | 27.6% |

The metric columns are separately averaged bootstrap estimates, so these
fractions are approximations. They are not the other teams' exact 750-row
submission counts. Different call fractions coexist among leading A4 entries;
the comparison cannot isolate the effect of call count because their rankings
and models also differ. N283T's estimated scored fraction agrees with its
exact 33.2% full-pool call fraction, but that agreement does not establish
an optimal threshold.

The saved equal-seven OOF uses the same 249/750 rate. Its scored rows cover
only the 3,584 labelled A4 training molecules, whereas the blind decision is
made on all 750 molecules before the scoring mask is applied. A labelled-only
rate sweep would be diagnostic, not a reliable way to choose a new blind count.
The project's [classification evaluation protocol](../classification-evaluation.md)
therefore calls for complete-pool OOF scores and selecting the model and call
rate jointly outside each evaluated fold. Keep 249 as the submitted reference;
test a prespecified rate grid under that protocol before considering a new
threshold or top-K count. Do not choose a rate from the one-time Interim score.

## Macro decomposition and interpretation

Compared with the **same** interim overall leader nova, N283T trails by
`.0877` MCC on D6 but leads by `.0513` on A4, giving a macro difference of
`−.0182`. This is why a strong A4 result puts N283T in the overall leading
tier despite the weaker D6 endpoint. Across the 114 entrants common to both
boards, overall rank Spearman is `.898`, with six names shared between the
two top-10 lists.

This snapshot measures the exact September 15 classification vector resent
before the interim cutoff ([RL0160](../research/2026-09.md#rl0160--resubmit-the-exact-september-15-classification-vector)).
Prior D6 replacement submissions improved some OOF views but lost on the Live
board; the full-test interim result does not supply compound labels to resolve
that discrepancy. Keep the submitted vector as the observed control. This
analysis changes no threshold, candidate or submission.
