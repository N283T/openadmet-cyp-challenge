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
The Live top-10 overlap with the Interim top 10 is only 3, and rank Spearman
among the 114 common entrants is `.756`. A4's half-test ranking is particularly
unstable, so the Live rank 1 should not be treated as the full-test result.

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
