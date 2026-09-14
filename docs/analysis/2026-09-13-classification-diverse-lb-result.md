# Diverse classification leaderboard result analysis

## Result

The September 13 logistic-blend D6 plus Caruana-reranked A4 submission scores
macro MCC **0.352961**, rank 13, versus **0.275867**, rank 48 for the previous
submission. Both endpoints improve:

| Endpoint | Previous MCC | New MCC | Change | Previous rank | New rank |
|---|---:|---:|---:|---:|---:|
| CYP2D6 | .176698 | .270958 | +.094260 | 44 | 21 |
| CYP3A4 | .375036 | .434963 | +.059927 | 46 | 15 |
| Macro | .275867 | .352961 | +.077093 | 48 | 13 |

The new macro score exceeds the public TabICL baseline of 0.324967 by 0.027994.

## Exact gain decomposition

Endpoint MCCs can be recombined exactly because the challenge macro is their
mean. This separates the three decisions made between submissions:

| Combination | Exact macro MCC |
|---|---:|
| Previous D6 + previous A4 | .275867 |
| Previous D6 + old A4 | .303676 |
| New D6 + old A4 | .350806 |
| New D6 + reranked A4 | .352961 |

Of the total +0.077093 macro improvement:

- Restoring A4 from 378 toward 249 calls contributes +0.027809, about 36.1%.
- The new D6 logistic rank blend contributes +0.047130, about 61.1%.
- Caruana reranking of the 249 A4 calls contributes +0.002155, about 2.8%.

The result therefore validates both main hypotheses, but not equally. D6 model
diversity is the largest gain. A4's main correction is its positive-call count;
Caruana selection adds a small real improvement over the exact old A4 endpoint.

## CYP2D6 behavior

The new D6 vector retains 29 of the prior 43 positive calls, removes 14, and adds
20 for 49 total. Its leaderboard precision and recall both improve:

| Metric | Previous | New | Change |
|---|---:|---:|---:|
| Precision | .2386 | .3080 | +.0694 |
| Recall | .2313 | .3469 | +.1155 |
| F1 | .2316 | .3224 | +.0909 |
| Accuracy | .8936 | .8988 | +.0052 |

This cannot be explained by merely calling six more compounds positive: precision
rises at the same time as recall. The Morgan logistic component successfully
reordered the sparse decision boundary. Using the previously inferred 7.06%
hidden prevalence, the bootstrap-weighted true-positive share rises by about
0.00816 while false-positive share rises only 0.00291. These are aggregate
operating-point estimates, not labels recovered for the 20 additions or 14
removals.

The D6 call count is already near the public top-10 median. The remaining gap to
the endpoint leader is more plausibly a ranking-quality problem than a large rate
error. Future D6 work should extend the successful low-capacity diversity idea—
fingerprint radius/count/chirality variants, regularized linear models and local
similarity consensus—using nested selection instead of broad threshold moves.

## CYP3A4 behavior

Relative to the previous 378-call A4 vector, the new candidate retains 247 calls,
removes 131 and adds only 2. Precision improves 0.4523→0.5666 and accuracy
0.6579→0.7538; recall falls 0.8213→0.6621. This is the intended trade: the prior
submission had overcalled positives against the hidden population.

The stronger comparison is with the old 249-call A4 endpoint. The new Caruana
ranking also calls 249 compounds and swaps 26 in each direction. It improves MCC
0.430653→0.434963, precision 0.5623→0.5666 and accuracy 0.7512→0.7538, while
recall is effectively unchanged (0.66225→0.66211). Caruana therefore made a small
false-positive reduction without sacrificing sensitivity. A4 is now rank 15 and
looks substantially closer to its useful operating point than before.

## Provenance and limitation

Corrected decomposition run `20260913T120454086409Z-12597afd` uses the exact three
submitted vectors and frozen pre/post leaderboard snapshots. Run
`20260913T120416878249Z-4e1164e4` is discarded because it treated a participant's
latest row as the earlier all-negative control after that participant had replaced
their submission. The correction fixes the priors to the historical all-negative
snapshot estimates: 7.06% D6 and 29.16% A4.

Leaderboard metrics are means over bootstrap evaluations on the confidently
scored subset. Score and rank changes are official; inferred positive-call rates
and confusion-matrix shares are weighted approximations. Exact per-compound hidden
labels cannot be recovered from these summaries.

## Decision

Keep this submission as the classification reference. Do not undo either endpoint
change. For the next round, prioritize D6 rank diversity while keeping the call
rate near the current region; treat A4 as a smaller reranking opportunity rather
than another prevalence correction problem.

