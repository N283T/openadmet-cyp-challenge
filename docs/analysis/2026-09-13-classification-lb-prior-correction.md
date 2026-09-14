# Classification leaderboard-prior analysis

## Conclusion

The two available leaderboard transitions are already informative enough to reject
a symmetric correction. The current CYP2D6 contraction was directionally correct,
whereas the CYP3A4 expansion was harmful. The lowest-risk next classification
candidate is therefore the exact current CYP2D6 vector combined with the old
CYP3A4 vector. Its endpoint-score mean is known from already scored submissions:
macro MCC **0.303676**, versus **0.275867** for the current submission. This still
trails the public TabICL baseline (**0.324967**).

No new candidate was uploaded. A more exploratory CYP3A4 reranking candidate was
prepared, but its advantage exists only in OOF projections and should not replace
the known hybrid without accepting additional model-selection risk.

## What can be inferred from the public leaderboard

Run `20260913T101151963949Z-892b1e27` uses the all-negative leaderboard row and
the published accuracy, precision and recall to infer operating points. Because
the leaderboard reports means over 1,000 bootstrap samples on the confidently
scored subset, these are bootstrap-weighted estimates, not recovered raw labels.

| Endpoint | train positive rate | inferred hidden positive rate | leader call rate | top-10 median call rate |
|---|---:|---:|---:|---:|
| CYP2D6 | 21.64% | 7.06% | 2.42% | 6.54% |
| CYP3A4 | 21.32% | 29.16% | 28.97% | 28.93% |

The hidden population therefore appears shifted in opposite directions for the two
endpoints. The old/current exact prediction transitions clarify what happened:

| Endpoint | old positives | current positives | change |
|---|---:|---:|---|
| CYP2D6 | 686 | 43 | current is an exact subset; 643 calls removed |
| CYP3A4 | 249 | 378 | 133 calls added and 4 removed |

An approximate decomposition of the published metric means estimates 5.73%
precision in the removed CYP2D6 block, below the 7.06% hidden base rate. Removing
it was beneficial. The net-added CYP3A4 block is estimated at 24.93% precision,
below its 29.16% hidden base rate. Adding it was harmful. These block precisions
are derived from products of bootstrap-mean metrics and must not be interpreted as
recovered compound labels.

The observed leaderboard scores agree: CYP2D6 rose from 0.001644 to 0.176698,
while CYP3A4 fell from 0.430653 to 0.375036. Combining the already scored current
D6 and old A4 endpoint results gives `(0.176698 + 0.430653) / 2 = 0.303676`.

## Fixed-rate and reranking screen

Run `20260913T103046578961Z-129fbc5d` fixes the number of positive calls using
leaderboard-derived rates, then compares saved-model ranks under the original five
outer folds. It also evaluates MCC after replacing the source prevalence with the
inferred hidden prevalence while retaining each candidate's OOF TPR and FPR. This
expected-target MCC is a projection, not an independent validation score.

For CYP2D6, the current vector remains competitive. At 53 calls, the current rank
has projected MCC 0.08250. A current/prior-weighted rank blend gives 0.08219 at 51
calls; at the current full-data count of 43, the three-way current/Caruana/weighted
blend gives 0.07733. The extra weighted model therefore does not provide a clear
reason to replace the current D6 ranking.

For CYP3A4, old-sized predictions are more plausible than the current 378 calls.
At 249 calls, the Caruana-AUROC ranking projects to MCC 0.43472 versus 0.41166 for
the historical ranking under the same OOF projection. This candidate swaps 26 of
the old vector's 249 calls. A historical/Caruana rank mean is more conservative,
swapping 14 calls. Neither reranking has a live result, so the exact old A4 vector
remains the evidence-backed choice.

Three validated, unsubmitted vectors are retained in the run output:

- `current_d6_caruana_a4_old_count`: current D6 (43) plus Caruana-ranked A4 (249).
- `current_d6_rankblend_a4_old_count`: current D6 (43) plus conservative blended
  A4 ranking (249).
- `best_top10_rate`: blended D6 (49) plus blended A4 (217); highest extrapolation
  risk because both endpoints change.

The exact current-D6/old-A4 hybrid had already been packaged separately and is not
recreated or uploaded here.

## Prior-weighted LightGBM

Run `20260913T102143243794Z-b9ebbcba` refits six ordinary LightGBM feature views
with positive-class weights equal to target odds divided by train odds (about 0.275
for D6 and 1.52 for A4). The first attempt,
`20260913T101424160645Z-d0e02f9c`, completed all 120 fits but failed while joining
ensemble output columns; it remains recorded. Run
`20260913T103042749697Z-185307ea` records deployment call counts.

The six-view mean gives D6 AUROC 0.57245, a small increase over the current
0.56619, but projected target MCC is only 0.05473 and its learned deployment
threshold calls just 8 positives. Individual D6 views range from 1 to 452 calls,
showing severe threshold instability. For A4, the six-view mean reaches AUROC
0.77270, below the historical model's 0.79328, and calls 346 positives. Direct
prior-weighted refitting is therefore rejected; its D6 rank is retained only as a
weak blend component in the fixed-rate screen.

## Decision

1. Preserve the current CYP2D6 predictions.
2. Prefer the exact old CYP3A4 predictions if making the next low-risk submission.
3. Keep the 249-call Caruana A4 reranking as the higher-risk alternative, not as a
   demonstrated improvement.
4. Do not adopt learned thresholds from the prior-weighted LightGBMs.
5. One more leaderboard result would be most informative if it isolates the A4
   ranking change while holding D6 exactly fixed.

