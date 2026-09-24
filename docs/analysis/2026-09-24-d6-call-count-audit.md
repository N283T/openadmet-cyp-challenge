# CYP2D6 TDI positive-call count is a model decision

The challenge accepts a boolean for every blind molecule and imposes no
49-positive requirement. Forty-nine is the count in the September 15
leaderboard-validated D6 package; it is a historical operating point, not an
estimated optimum. The September 24 four-swap probe also used 49 so that its
model-score change could be compared at the same call budget. That choice does
not establish that 49 should constrain future candidates.

Recorded run `7e2db311` evaluates the September 15 control analog and the
75/25 control/AP rank blend on all 6,146 held-fold training molecules. Each
model calls the same fraction corresponding to 15, 25, 35, 49, 65, 80, 110 or
150 of 750 blind compounds, then MCC is computed only where the D6 label is
supplied. The table uses the estimated 7.06% scored-population prior:

| Model | UMAP best count and projected MCC | Scaffold best count and projected MCC |
|---|---:|---:|
| September 15 control analog | 49, `.10393` | 110, `.14048` |
| 75/25 rank blend | 35, `.11404` | 80, `.15748` |

These are **in-sample maxima on the rate grid**, not unbiased performance
estimates. The best counts within each split remain the same under the 5% and
10% prior sensitivity scenarios; changing the prior does not resolve the
UMAP-versus-scaffold discrepancy. The rate chosen outside each held fold is also
unstable: for the
control, UMAP chooses 4%–6.53% and scaffold chooses 15% on every fold. The
blend chooses 4%–10% on UMAP and 10%–15% on scaffold (run `2fe73ea9`). The
full-pool OOF cannot identify a single safe replacement call count, especially
given the acquisition and label-coverage gap described in the
[analog-expansion analysis](2026-09-24-d6-analog-expansion-validation.md).

The blind score audit confirms that the saved September 15 score's top 49
exactly reproduce the accepted submission. Changing its call count to 35
removes 14 accepted calls; 65 adds 16; 110 adds 61. For the 75/25 blend, 49
exchanges four calls each way, while 80 retains all 49 and adds 31. None of
those changed compounds has a known blind class label, so these are decision
counts, not predicted LB gains.

**Decision:** evaluate model and call count jointly on complete-pool,
fold-excluded predictions; include the historical 49 count as a comparator.
Do not infer a new call count from the pooled OOF maximum. Keep the September
15 accepted vector as the LB-validated control until a better operating point
has external evidence. The first attempted run `13b3b1c9` failed because a
join reordered blind rows before the saved-score equality check; `7e2db311`
sorts by compound ID and verifies the exact submitted top-49 vector.
