# CYP interim leaderboard: September 26 snapshot

## Evidence and scope

Recorded analysis run `shared/20260926T002255140309Z-2e7b1859` uses all 16
downloaded leaderboard CSVs indexed in
[`artifacts/tracked/leaderboard-snapshot-20260926/manifest.json`](../../artifacts/tracked/leaderboard-snapshot-20260926/manifest.json).
The run snapshots each CSV as an explicit input and checks its SHA-256. It compares
the September 26 Live and Interim boards for both inhibition tracks, including
their macro and every endpoint table. The Interim board scores the September 24
submission against the full test set once; Live scores one half and remains active.
This analysis contains no hidden compound-level labels or competitor predictions.

## N283T: full-test interim result

| Track / endpoint | Live rank, score | Interim rank, score | Interim leader, score |
|---|---:|---:|---:|
| Regression overall (MA-ST-RAE ↓) | 27, .4662 | **19, .4494 (Tier 5)** | preheat-to-450, .3721 |
| CYP1A2 ST-RAE ↓ | 40, .5085 | 32, .4417 | wbot, .3752 |
| CYP2C9 ST-RAE ↓ | 13, .3594 | **5, .3587** | preheat-to-450, .3130 |
| CYP2D6 ST-RAE ↓ | 31, .5655 | 26, .6104 | preheat-to-450, .4568 |
| CYP3A4 ST-RAE ↓ | 37, .4313 | 25, .3868 | preheat-to-450, .3331 |
| TDI overall (MA-MCC ↑) | 12, .4193 | **5, .3798 (Tier 1)** | nova, .3980 |
| CYP2D6 MCC ↑ | 31, .3086 | 25, .2573 | nova, .3450 |
| CYP3A4 MCC ↑ | **1, .5299** | **4, .5023** | EleanorRigby, .5207 |

N283T's submission timestamps match between the two boards: regression
September 24 05:44 UTC and classification September 24 06:24 UTC. The rank
changes therefore compare the same N283T submission across scoring populations,
although the competing population also changes. Interim has 219 regression and
114 classification entries, compared with 226 and 127 in this Live snapshot.
The blog reports rule enforcement and Live remains open; the CSVs alone do not
identify every cause of those count differences.

## Where the aggregate gaps arise

Against the **same** interim regression leader, preheat-to-450, N283T's
endpoint ST-RAE gaps are .0561 (1A2), .0457 (2C9), **.1536 (2D6)**, and
.0537 (3A4). Their mean is the macro gap .0773. CYP2D6 accounts for 49.7% of
that four-target gap. This agrees with the challenge blog's warning that D6 is
affected by a training/test activity distribution shift, but these leaderboard
values cannot establish the mechanism for N283T's errors. The D6 full-test score
is worse than the Live score (.6104 vs .5655), while the other three endpoint
scores are lower on the full test. Scores from different populations must not be
treated as paired compound-level improvements or losses.

Against interim TDI leader nova, N283T is .0877 MCC behind on D6 and .0513
ahead on A4, leaving a macro gap of .0182. On D6, N283T's interim
precision/recall are .2878/.3224; nova's are .4725/.3006. This is consistent
with N283T's D6 positive calls having lower precision at this operating point,
but without predictions or labels it does not identify which calls to change.
On A4, N283T's interim precision/recall are .6393/.6848, versus nova's
.6472/.5772. All 26 leading TDI macro entries are Tier 1 per the organizer;
rank 5 does not establish a significant difference from rank 1.

## How stable are Live ranks?

Across users present on both boards, Spearman correlation of Live and Interim
overall ranks is .986 for regression (215 common entrants) and .898 for TDI
(114 common entrants). The top-10 overlap is 8 and 6, respectively. At the
endpoint level, A4 TDI is less stable: rank correlation .756 and only 3 of the
Live top 10 appear in the Interim top 10. These are descriptive comparisons;
they are not significance tests. Live rank is a useful signal, but the Interim
result is the more direct observation of full-test performance at the cutoff.

## Decision

Keep the existing six submitted endpoint vectors as the measured reference.
Prioritize independent, fold-safe D6 regression and D6 TDI diagnostics and
candidate evaluation. Retain the successful A4 TDI vector as a control. Do not
choose a new model or threshold from these aggregate full-test scores alone:
the one-time Interim release exposes information about the final evaluation
set, and repeat tuning against it can overfit. Any replacement still needs a
recorded experiment, complete-pool evaluation where relevant, and a separate
submission decision. Prepare the required final model report before the
November 3 23:59 UTC deadline.
