# Classification leaderboard transfer benchmark

## Recorded evidence

Run `classification/20260924T031257156540Z-4febdf99` snapshots seven exact
accepted classification submission files and their post-scoring official D6 and
A4 leaderboard CSVs. The run validates identical 750 molecule IDs and SMILES,
the submitting team's timestamp on both endpoint boards, and one call vector
per date. Its `states.csv`, `transitions.csv`, and `changed-calls.parquet` are the
fixed reference for subsequent OOF protocol comparisons. The initial run
`0b664a0f` produced the same observations with less convenient CSV formatting;
it remains an immutable recorded result.

| Transition | D6 MCC change | D6 removed / added calls | A4 MCC change | A4 removed / added calls |
|---|---:|---:|---:|---:|
| Sep 10 → Sep 12 | +.17505 | 643 / 0 | −.05562 | 4 / 133 |
| Sep 12 → Sep 13 | +.09426 | 14 / 20 | +.05993 | 131 / 2 |
| Sep 13 → Sep 14 | +.00761 | 6 / 6 | +.01465 | 16 / 16 |
| Sep 14 → Sep 15 | +.03006 | 9 / 9 | +.08029 | 18 / 18 |
| Sep 15 → Sep 18 | −.19891 | 27 / 27 | 0 | 0 / 0 |
| Sep 18 → Sep 23 | +.12510 | 21 / 21 | 0 | 0 / 0 |

The first two D6 transitions also change the number of positive calls
substantially (686 → 43 → 49), so they are joint model/rate checks. From Sep 13
onward D6 always has 49 positive calls, making those four transitions cleaner
tests of the model's decision boundary. A4 holds 249 calls from Sep 13 onward;
the last two A4 entries are exact-vector controls, not independent scoring
successes. The seven states produce six chronological comparisons for each
endpoint, but all use the same hidden test set and cannot be treated as six
independent samples.

## Implication for an improved OOF

The primary D6 development check is whether a prespecified, genuinely refit
validation protocol recognizes the Sep 14 → Sep 15 improvement and ranks
the Sep 15 vector above the Sep 18 vector. Reserve Sep 18 → Sep 23 as a
held-back check of the protocol decision; its result is already known, so this
is a procedural guard against further tuning rather than a blind holdout. The
Sep 13 → Sep 14 gain is small relative to public score uncertainty and should
not be used as a pass/fail threshold.

Replay each historical recipe on common outer folds and rank the entire held-out
molecule pool, including unlabeled molecules, before computing MCC on labeled
ones. Save fold-safe scores with null labels for unscored rows, and evaluate
the historical submitted call counts first. Then select model and call rate
together using only other folds. The existing 1,497-row D6 OOF covers only
24.36% of the 6,146-molecule pool; its top-call boundary cannot serve as the
new benchmark. Existing UMAP OOF scores cannot be reassigned to an alternative
split. This run establishes the LB target and exact changed calls; it does not
yet establish a better OOF protocol or a new candidate.
