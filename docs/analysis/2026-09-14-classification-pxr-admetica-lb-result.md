# ComplementNB D6 plus PXR/Admetica A4 leaderboard result

## Conclusion

The September 14 classification submission improves both endpoints. Macro MCC
rises from 0.352961 to 0.364094. The A4 revision supplies 65.8% of that gain and
the D6 revision supplies 34.2%. Both conditional CV estimates had wide intervals
crossing zero, but both predicted the correct direction and understated the
observed leaderboard gain.

This is useful evidence for two distinct modeling choices. For D6, a small
Complement Naive Bayes rank contribution improves the already successful weak
Morgan-logistic ensemble. For A4, a lower-correlation classical auxiliary bank
adds useful boundary information to Caruana. It does not establish PXR activity
as the causal source of the A4 gain because prior ablations identify Admetica as
the larger compact contribution and PXR as smaller and split-dependent.

## Leaderboard movement

| Metric | Previous | New | Change | Previous rank | New rank |
|---|---:|---:|---:|---:|---:|
| Macro MCC | 0.352961 | 0.364094 | +0.011133 | 16 | 13 |
| CYP2D6 MCC | 0.270958 | 0.278572 | +0.007614 | 23 | 23 |
| CYP3A4 MCC | 0.434963 | 0.449616 | +0.014653 | 18 | 13 |

Because macro MCC is the arithmetic mean of the two endpoint MCCs, D6 contributes
0.003807 macro points and A4 contributes 0.007326. The new macro is 0.039127 above
the public TabICL baseline of 0.324967.

The cumulative improvement across the last two submissions is substantial:

| Submission | Macro MCC | D6 MCC | A4 MCC |
|---|---:|---:|---:|
| September 12 | 0.275867 | 0.176698 | 0.375036 |
| September 13 | 0.352961 | 0.270958 | 0.434963 |
| September 14 | 0.364094 | 0.278572 | 0.449616 |

## Endpoint interpretation

### CYP2D6

The 2:1 live-D6/ComplementNB rank blend retains 43 of the 49 prior positive calls,
removes six and adds six. Public precision rises from 0.30797 to 0.31967 while
recall remains 0.34685; accuracy and F1 also rise. This supports the narrow claim
that the generative bit-frequency view cleans up the sparse decision boundary
without losing aggregate recall. It does not support broad replacement of the
current ensemble or more model-family screening.

### CYP3A4

The 80:20 live-Caruana/Morgan-PXR-Admetica rank blend retains 233 of 249 prior
calls, removes 16 and adds 16. Precision rises from 0.56662 to 0.57560 and recall
from 0.66211 to 0.67354; accuracy and F1 also rise. The fixed call count combined
with improvement in every reported endpoint metric is consistent with genuinely
better boundary ranking rather than a prevalence-only correction.

The earlier feature ablation remains important: PXR scalar predictions alone were
weak, Admetica supplied the larger compact gain, and strict PXR-Chemprop was too
correlated with Caruana to improve it. The leaderboard therefore validates the
combined lower-correlation auxiliary bank, not a standalone mechanistic PXR claim.

## CV-to-leaderboard calibration

| Endpoint | Conditional CV delta | 95% interval | LB delta | LB/CV ratio |
|---|---:|---:|---:|---:|
| CYP2D6 | +0.00503 | [-0.01905, 0.03087] | +0.00761 | 1.51 |
| CYP3A4 | +0.01074 | [-0.00363, 0.02371] | +0.01465 | 1.36 |

Both directions transferred, and both observed effects are moderately larger than
their CV point estimates. This increases confidence that fixed-rate, low-capacity,
low-correlation rank additions are useful. It should not be used as a general
1.4x correction factor: there are only two selected observations and both were
chosen after bounded screens.

## Decision

Keep the September 14 vector as the live classification reference. Do not tune
the A4 blend weight upward from one leaderboard observation, and do not add the
high-correlation strict PXR-Chemprop family. For D6, prioritize prespecified small
changes that add a distinct low-capacity count or linear view. For A4, prioritize
new genuinely low-correlation phenotype or assay blocks and require fixed-rate
cross-fit support before another submission.

Run `20260914T022414059180Z-253152e5` contains the exact endpoint metrics,
call-transition table, changed-compound list, CV/LB comparison and summary. Public
leaderboard metrics are bootstrap means over confidently scored subsets, so their
precision, recall and accuracy averages cannot be inverted into a single exact
confusion matrix or used to label individual swapped compounds as correct.
