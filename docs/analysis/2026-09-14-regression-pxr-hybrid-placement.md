# Strict PXR-Chemprop regression transfer and four-endpoint placement scenarios

## PXR representation transfer

The challenge-disjoint three-seed PXR Chemprop asset from classification run
`519bbb46` supplies a fixed 384-dimensional embedding for every current CYP
compound. Regression run `c33a4941` fits the established deterministic LightGBM
readout independently for all four direct-inhibition endpoints on the retained
five-fold parent CV. Audit `5d383004` verifies all 6,384 OOF keys and comparisons
with the older aligned model bank.

The PXR embedding has regression signal, but it is not a strong standalone
replacement. Run `c00f292d` reconstructs the actual simplified operational OOF
ensemble and tests fixed PXR weights 0, 10%, 20%, 25% and one third. Failed
precursor `d3f05666` completed the calculations but could not serialize a nested
list to CSV; no result was silently discarded.

| Endpoint | Baseline ST-RAE | Best full-OOF PXR weight | Best ST-RAE | Fold-selected ST-RAE | Decision |
|---|---:|---:|---:|---:|---|
| CYP1A2 | 0.687942 | 20% | 0.683686 | 0.686911 | Do not stack on the planned large placement |
| CYP2C9 | 0.463103 | 10% | 0.462422 | 0.463329 | Reject; negligible and not cross-fit stable |
| CYP2D6 | 0.828397 | 10% | 0.826175 | 0.827346 | Do not stack on the planned large placement |
| CYP3A4 | 0.398316 | 0% | 0.398316 | 0.398316 | Reject; every nonzero weight worsens |

The conditional full-OOF bootstrap intervals versus zero include zero for the
small 1A2, 2C9 and 2D6 changes. At 10%, CYP3A4 worsens to 0.400600. The result
does not support using strict PXR Chemprop merely to ensure that 2C9 and 3A4
change in the next upload.

## Hybrid source-moment scenarios

Run `46d505e2` applies the already selected full 1A2 and D6 placements to the
exact simplified test vector, then moves 2C9 and 3A4 by 25%, 50% or 100% toward
the archived source-author targets. The formula is
`raw + strength * (full affine placement - raw)`. Two precursor runs
`300c09f6` and `e6fe44b8` failed during immutable alias input-name resolution
before any transformation. All three completed CSVs pass the official 750-row
identity, schema and finite-value validator. They are scenarios, not registered
candidates or prepared submissions, and no score is assigned.

| Endpoint | Raw mean / SD | Full target mean / SD | 25% mean absolute change | 50% | 100% |
|---|---:|---:|---:|---:|---:|
| CYP2C9 | 4.922 / 0.712 | 4.830 / 0.923 | 0.045 | 0.090 | 0.180 |
| CYP3A4 | 4.768 / 0.872 | 4.880 / 1.082 | 0.050 | 0.101 | 0.202 |

Even full placement is modest relative to the planned CYP2D6 mean absolute
change of 1.119 and CYP1A2 change of 0.514. If the objective is expected-score
conservatism, retain the 25% scenario. If the objective is to make the 12-hour
leaderboard probe maximally informative and accept source-assumption risk, use
`all_full_hybrid_targets`. This is the preferred exploratory scenario under the
user's stated high-risk tolerance. It changes every endpoint while preserving
each endpoint's rank ordering.

## Live submission

The user explicitly approved the all-full exploratory scenario. Candidate
`20260914T024735834387Z-f87e9914` and prepared submission
`20260914T024743176950Z-8ef21e06` preserve SHA256
`90dda228f2e027983891579ab7b45a23a28b2ec09956a6cbf5a662417e446086` and pass
the official dry run. Before upload, the sender saved the MA and all four
endpoint leaderboard snapshots. Attempt `20260914T024810182582Z-83cfb34e`
received a confirmed `sent` response at 2026-09-14 02:48 UTC; the service said
the regression predictions were received for processing. Leaderboard reflection
and a score are not yet claimed.

## Leaderboard result

Fresh MA and endpoint snapshots at 2026-09-14 03:22 UTC confirm reflection.
Submission `8ef21e06` is recorded as accepted/scored with the MA snapshot as
evidence.

| Endpoint | Previous ST-RAE | Placed ST-RAE | Change | Previous rank | New rank |
|---|---:|---:|---:|---:|---:|
| CYP1A2 | 0.685494 | 0.526455 | -0.159038 | 69 | 38 |
| CYP2C9 | 0.369556 | 0.366939 | -0.002617 | 12 | 12 |
| CYP2D6 | 1.353707 | 0.565520 | -0.788187 | 83 | 20 |
| CYP3A4 | 0.468640 | 0.444983 | -0.023657 | 38 | 31 |
| Macro | **0.719349** | **0.475974** | **-0.243375** | **69** | **20** |

Lower ST-RAE is better. Every endpoint improves. The main effect is the intended
CYP2D6 distribution correction; its final 0.565520 is strikingly close to the
source author's reported 0.5653 at the same 3.57/0.90 placement, despite the
earlier warning that the score was not transferable. CYP1A2 also captures most
of the hypothesized gap. Full 2C9 placement is nearly neutral but positive, and
full 3A4 placement produces a useful smaller gain.

Because all transforms are strictly affine with positive scale and no floor was
activated, rank metrics are exactly unchanged: macro Spearman remains 0.724695
and Kendall remains 0.550403. The aggregate MAE improves 0.895630 to 0.657838
and R² improves 0.221760 to 0.591026. This is strong evidence that blind
population placement, rather than molecular ordering, accounted for most of the
remaining regression leaderboard gap in the previous vector.

## Recorded diagnostic and counterfactuals

Run `20260914T041058804447Z-a2287e3c` reads the immutable pre-submission and
post-reflection MA/endpoint snapshots plus the exact old and new prediction
files. It reconstructs the macro from endpoint scores to within `5.6e-17`,
verifies every endpoint improved, verifies identical ranks for all 3,000 paired
predictions and reconstructs all four affine maps with maximum error `6.3e-15`.
Failed precursors `83ecb84a`, `f33411cd` and `cca83a7e` preserve input-extension
and immutable-filename resolution errors before the analysis completed.

The additive endpoint improvements sum to four times the macro improvement.
CYP2D6 contributes **80.96%**, CYP1A2 **16.34%**, CYP3A4 **2.43%** and CYP2C9
**0.27%**. This is an accounting decomposition of the isolated affine changes,
not merely a correlation with the leaderboard movement.

Because the leaderboard macro is the exact mean of separately scored endpoints,
the endpoint substitution counterfactuals are exact:

| Reflected endpoint placements | Counterfactual macro ST-RAE | Improvement vs previous |
|---|---:|---:|
| None | 0.719349 | 0 |
| CYP1A2 only | 0.679590 | 0.039760 |
| CYP2D6 only | 0.522302 | 0.197047 |
| CYP1A2 + CYP2D6 | 0.482543 | 0.236806 |
| All four | **0.475974** | **0.243375** |

Thus changing CYP2C9 and CYP3A4 was worth an additional **0.006568 macro
ST-RAE** beyond the original 1A2+D6 candidate. It accounts for only 2.70% of the
total gain, but both endpoints moved in the correct direction and the extra
change was not wasted.

## Implication for the next regression round

Treat the four placement constants as the new fixed test-time baseline. The
single reflected comparison strongly validates the direction and nearly exactly
validates the source D6 score, but hidden labels remain unavailable and therefore
do not establish that these constants are the exact endpoint optima. Another
fine placement probe has lower expected information value than improving model
ordering.

The current leaderboard gaps to the endpoint leaders in the same snapshots are
0.122494 for D6, 0.109310 for 1A2, 0.083475 for 3A4 and 0.036372 for 2C9. D6 is
also the clear ordering bottleneck: its Spearman remains 0.431249, versus
0.779924/0.856334/0.831273 for 1A2/2C9/3A4. Prioritize a genuinely different D6
ranker, then 1A2 diversity. Preserve 2C9 and the new placement unless a new model
has strong matched evidence. The strict PXR-Chemprop regression readout tested
in this round is not that D6 ranker and remains rejected for 2C9/3A4 deployment.
