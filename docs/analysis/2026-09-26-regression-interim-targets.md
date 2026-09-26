# Regression interim leaderboard by target

## Evidence

Recorded run `regression/20260926T002955923155Z-503b1907` analyzes the frozen
September 26 Live and Interim regression leaderboard CSVs for overall and all
four endpoints. Its ten explicit input snapshots and hashes trace to the
[source manifest](../../artifacts/tracked/leaderboard-snapshot-20260926/manifest.json).
The primary ranking metric is endpoint ST-RAE (lower is better). MAE and
Spearman's rho are secondary descriptive metrics. Only aggregate leaderboard
rows are available; this analysis cannot identify compound-level errors or
estimate a paired model change on the hidden test set.

## N283T among 219 interim regression entries

| Target | ST-RAE rank / score | Leader / score | Gap to leader | MAE rank | Spearman rank | Live ST-RAE rank |
|---|---:|---:|---:|---:|---:|---:|
| CYP1A2 | 32 / .4417 | wbot / .3752 | .0665 | 21 | 35 | 40 |
| CYP2C9 | **5 / .3587** | preheat-to-450 / .3130 | .0457 | **2** | 6 | 13 |
| CYP2D6 | 26 / .6104 | preheat-to-450 / .4568 | **.1536** | 36 | **89** | 31 |
| CYP3A4 | 25 / .3868 | preheat-to-450 / .3331 | .0537 | 18 | 27 | 37 |

The rank within the same interim endpoint table is the meaningful position;
raw ST-RAE values across different targets should not be interpreted as a
common activity-scale error. Relative to each target's tenth-place ST-RAE,
N283T trails by .0336 (1A2), is already top five (2C9), trails by .0796
(2D6), and by .0239 (3A4).

### CYP1A2

N283T is rank 32 by ST-RAE, while MAE is rank 21 and Spearman is rank 35.
The gap to the leader is .0665. The same submitted vector moves from Live
rank 40 to Interim rank 32, and its ST-RAE is .5085 on Live versus .4417 on
the full test. These are different scoring populations, not a paired gain.
The current submitted 1A2 vector is the September 18 foundation-v8 choice,
selected before the interim cutoff for its observed Live improvement; see
[RL0159](../research/2026-09.md#rl0159--endpoint-wise-leaderboard-best-regression-submission).

### CYP2C9

N283T is rank 5 by ST-RAE, rank 2 by MAE and rank 6 by Spearman. Its MAE
is .4280, versus .4155 for the interim ST-RAE leader. The ST-RAE gap of .0457
is real on this board, but the strong secondary metrics suggest that the
submitted 2C9 vector is already competitive. Preserve it as a control rather
than changing it to chase a few leaderboard positions.

### CYP2D6

The ST-RAE gap to the leader is .1536, substantially larger than on the other
targets. N283T's full-test ST-RAE is .6104 versus .5655 on Live. This target
alone contributes 49.7% of the macro ST-RAE gap to the same overall leader.
The rank by Spearman is 89 despite ST-RAE rank 26: N283T's rho is .4053,
compared with .5982 for the interim endpoint leader. MAE ranks 36. Thus the
available aggregate metrics point to weak **ordering** as well as point-error
room; a simple output-scale adjustment alone cannot raise Spearman. This is a
diagnostic clue, not proof that a ranking-focused model will improve ST-RAE.
The [official halfway post](https://openadmet.ghost.io/halfway-through-the-openadmet-cyp-challenge/)
explains the sampling asymmetry precisely: the **same 750 test molecules** are
evaluated for all four isoforms, but the analog expansion selected the top 25
training hits for each of 1A2, 2C9 and 3A4, then ten similar analogs per hit.
D6 hits were not a selection anchor, so the D6 test measurements skew toward
weaker activity relative to its training data. This can change the activity
range and affect a rank correlation; N283T's low D6 Spearman cannot be
attributed to its model alone from these aggregate scores. Inspect labelled
test-like analog groups before changing the model.

### CYP3A4

N283T is rank 25 by ST-RAE, rank 18 by MAE and rank 27 by Spearman. Its
ST-RAE .3868 is .0537 behind the leader but only .0239 behind tenth place.
The Live-to-Interim rank moves 37 to 25. The current submitted 3A4 vector
remains the appropriate reference: the September 18 replacement worsened its
Live ST-RAE, and the endpoint was rolled back before the interim cutoff
([RL0159](../research/2026-09.md#rl0159--endpoint-wise-leaderboard-best-regression-submission)).

## Ranking stability and next diagnostic

Among 215 common entrants, the Live-to-Interim rank Spearman correlations are
.989 (1A2), .987 (2C9), .984 (2D6), and .978 (3A4). The boards are broadly
consistent, although individual ranks and target scores move. Candidate
selection should remain based on fold-safe internal evaluation, not repeated
adjustment to this one-time full-test reveal.

For the next regression diagnosis, start with D6: compare the exact submitted
prediction vector against its keyed out-of-fold predictions on held-out,
test-like analog groups, including Spearman, ST-RAE, activity strata and
prediction spread. Then review 1A2 and 3A4. There is no evidence from these
aggregate tables to replace 2C9 or submit a new regression candidate now.
