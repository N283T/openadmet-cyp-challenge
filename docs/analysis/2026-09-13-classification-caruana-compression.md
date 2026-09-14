# Classification A4 Caruana compression

## Question

The deployed CYP3A4 TDI rank uses six effective Caruana members. Their test
predictions contain substantial redundancy, especially the two Jeremy PCA256
readouts (Spearman 0.963). Run `20260913T145411164843Z-7bc3d711` asks how many
members are needed to reproduce the existing ensemble rather than optimize a
new classifier.

Every nonempty subset of the six members was enumerated. For each subset,
nonnegative weights summing to one were fitted only to minimize OOF probability
MSE against the original Caruana prediction. The blind labels were never used.
Fidelity was then measured on the blind predictions, including changes at the
deployed top-249 call boundary.

## Best subset at each size

| Size | Members and OOF-fitted weights | Test Spearman | Top-249 swaps each way | OOF MCC |
|---:|---|---:|---:|---:|
| 1 | historical stack 1.000 | 0.92250 | 26 | 0.37950 |
| 2 | historical stack 0.622; Jeremy PCA256 TabICL 0.378 | 0.99472 | 7 | 0.40987 |
| 3 | historical stack 0.526; Jeremy PCA256 TabICL 0.364; external-class Mordred LGBM 0.111 | 0.99791 | 4 | 0.39975 |
| 4 | historical stack 0.501; Jeremy PCA256 TabICL 0.219; Jeremy PCA256 TabPFN 0.172; external-class Mordred LGBM 0.109 | 0.99929 | 2 | 0.40698 |
| 5 | preceding four plus external-class PCA256 TabPFN 0.044 | 0.99974 | 0 | 0.40843 |
| 6 | original six-member weights | 1.00000 | 0 | 0.40987 |

The five-member result shows that the 2% Morgan LGBM member can be removed
without changing any blind top-249 call. Removing both low-weight tail members
still leaves only two swaps each way. More aggressive distillation is also
possible: the two-member approximation preserves the blind ordering closely and
changes seven calls each way while reducing six model types / 30 outer-fold
readouts to two types / ten readouts.

The historical stack itself depends on frozen neural auxiliary features, so the
readout count does not represent total from-scratch training complexity.
Nevertheless, the result establishes that four of the six top-level members are
not required to reproduce the deployed rank closely.

## Decision

Treat the two-member version as the cleanest compression probe and the
four-member version as the conservative approximation. The three-member result
is not a clearly better compromise: it saves only one further type than the
four-member version, has lower descriptive OOF MCC, and the two-member version
is already simpler with comparable boundary fidelity. Do not replace or submit
the live vector solely from this analysis; matching a prediction is not evidence
that the compressed recipe will improve hidden labels.

## Unhashed-radius addition

Runs `499f49da` and `2ca7b030` test fixed late-rank additions of the A4 radius-0/1
unhashed block to the two-member compression. The previously preferred internal
block, an equal rank mean of ordinary Morgan and exact unhashed counts, lowers
AUROC, AP and fixed-rate MCC at every tested weight. It should not be added.

Pure unhashed counts behave differently when given only 1/16 of the final rank:

| Recipe | AUROC | AP | Target-prior MCC | Test Spearman to full | Top-249 swaps each way |
|---|---:|---:|---:|---:|---:|
| Full six-member Caruana | 0.80974 | 0.52512 | 0.44394 | 1.00000 | 0 |
| Two-member compression | 0.80681 | 0.52215 | 0.44394 | 0.99472 | 7 |
| 15:1 compression / pure-unhashed | 0.80561 | 0.51710 | 0.44700 | 0.99494 | 11 |

The MCC gain is only 0.00307, its grouped conditional interval versus the full
ensemble is [-0.01059, 0.01629], and ranking metrics decline. Thus this is a
plausible diversity probe, not validated improvement. If tested on the
leaderboard, use exactly the prespecified 15:1 recipe rather than tuning the
boundary further. It retains a compact three-type ensemble: historical stack,
Jeremy PCA256 TabICL and pure radius-0/1 unhashed logistic regression.
