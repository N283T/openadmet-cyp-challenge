# Classification: whole intermediate-depth members in submitted TDI ensembles

## Design

The current CYP2D6 TDI score is a fold-selected core with a Jeremy PCA/raw
mean block, followed by logistic and optional tail stages. The current CYP3A4
TDI score is a seven-member equal mean. Run `fd388a59` reconstructs their
saved labelled UMAP OOF from atomic inputs (maximum prediction differences
`<1e-12` for D6 and `4.25e-9` for A4). The lower-precision A4 source accounts
for the latter difference.

CheMeleon depths 1, 3 and 6 use matched LightGBM OOF from `8017d516`.
Depth 1 additionally has PCA256 TabPFN2.6 and TabICL classifier readouts,
fit on the same outer UMAP folds in `933f232c` and `c7a759af`. The readout
fits use only labelled training molecules outside each held fold. All
candidates are complete members: a multi-depth block first averages its
depth-specific predictions, then enters the current ensemble once.

For D6, addition extends the selected core by one member before its existing
rank and tail stages; replacement substitutes for the existing Jeremy
PCA/raw block. For A4, addition creates an eight-member equal mean;
replacement substitutes for the historical CheMeleon stack in the selected
seven. No arbitrary outer blend fractions are tested.

This is a **labelled-only** OOF screen: 1,497 D6 and 3,584 A4 molecules.
At the historical labelled-fold call rates, each scenario makes 97 D6 or
1,190 A4 positive calls. Those counts are evaluation points, not fixed
requirements for a future model. Missing labels on the 6,146-molecule pool
are not turned into negatives.

## Results

| Endpoint and member operation | AUROC | AP | Target-prior MCC |
| --- | ---: | ---: | ---: |
| D6 current | .580764 | .296860 | .110109 |
| D6 add depth 1 LightGBM | .579223 | .294497 | .075258 |
| D6 add depth 1 TabPFN | .575713 | .294008 | **.115312** |
| D6 add depth 1 TabICL | .580515 | .295208 | .110109 |
| D6 replace Jeremy block with depth 3 LightGBM | .573220 | .293655 | .110109 |
| A4 current equal seven | .805784 | .522733 | **.465398** |
| A4 add depth 6 LightGBM | .804887 | **.526930** | .456205 |
| A4 add depth 1+6 block | .805753 | .524495 | .465398 |
| A4 add depth 1 TabPFN | .803750 | .519522 | .457738 |
| A4 add depth 1 TabICL | .805281 | .515506 | .456205 |

All other tested depth additions and replacements are in the tracked
`ensemble-metrics.csv`. None improves A4 target-prior MCC. D6 depth-1
TabPFN changes eleven positive calls at this labelled boundary. Its MCC
increase is `.005203`, while both AUROC and AP fall. Conditional UMAP-group
bootstrap run `acfe1ea9` (3,000 replicates, 50 groups) gives the paired D6
MCC delta interval `[-.016270,+.027342]`; 68.7% of resamples improve.
For A4 depth 1+6, the point MCC is tied and the interval is
`[-.010615,+.011102]`.

Standalone depth-1 TabPFN AUROC is `.521684/.690150` for D6/A4; TabICL
is `.540058/.707851`. Their ensemble effect is therefore a narrow call
boundary result, not a broad ranking improvement. The first TabICL attempt
`e76699f2` fit all folds but failed a final assertion hard-coded for three
views; `c7a759af` reran the same one-view experiment with the expected
count derived from configuration. The first ensemble reconstruction attempt
`8a156350` stopped on an overly strict A4 floating-point identity assertion;
the corrected tolerance was verified in `985d41ee`, then the full screen in
`fd388a59`.

## Decision and limit

Keep the submitted D6 and A4 call vectors. The only positive MCC point
estimate, D6 depth-1 TabPFN, is small and uncertain on labelled molecules,
with lower AP/AUROC. A4 depth-6 improves AP but loses at the operating
boundary. The historical D6 labelled-only OOF has previously failed to
predict blind leaderboard direction. These readouts do not cover the complete
held-out molecular pool, so they cannot support call-rate selection or a
deployment decision under the [classification protocol](../classification-evaluation.md).
No candidate or submission was prepared. Compact OOF and metric evidence is
tracked under `artifacts/tracked/classification-depth-whole-members-20260925/`.
