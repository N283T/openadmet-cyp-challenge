# CYP2D6 support-aware routing

## Conclusion

Nearest-neighbor support changes which saved model looks strongest, but a
fold-excluded support gate does not improve the fixed 49-call D6 boundary. The
larger finding is that the current UMAP OOF problem is much more structurally
difficult than the blind test set. Retain the live classifier, reject the current
gate, and add a test-support-matched secondary validation design rather than
tuning the gate further on these folds.

## Design

Run `20260913T132834892734Z-5f9af66e` compares ten saved D6 ranks without fitting
new endpoint models: live and original ensembles, Tanimoto kNN, Mordred and frozen
embedding CatBoost, ComplementNB and RBF-SVM blends, radius-3 count logistic,
Chemprop-mixed logistic and radius-3/Chemprop early fusion.

For every OOF molecule, support is its maximum Morgan Tanimoto similarity to the
other four UMAP folds. Candidate gates route the lowest 25%, 50% or 75% to a
global model and the remainder to an FP/local model. Each heldout fold selects
gates using only the other four OOF folds. Test routing uses support percentile,
not the raw threshold, because full-data test models have a larger neighbor pool.

## Support strata

Below similarity 0.3, the ComplementNB blend has AUROC/AP 0.565/0.316 versus
0.558/0.298 live, while radius-3 count has the highest AUROC 0.587 but lower AP
0.286. From 0.3 to 0.4, ComplementNB has the highest AUROC 0.607 but the original
ensemble has the highest AP 0.310. Above 0.4 there are only 99 OOF molecules, and
only 18 above 0.5, so apparent kNN or RBF gains there are too sample-poor to use.
There is no stable low-support winner.

## Fold-excluded gate result

| Method | Source MCC | Projected target MCC | AUROC | AP |
|---|---:|---:|---:|---:|
| Live reference | 0.13842 | 0.09485 | 0.57672 | 0.29237 |
| Support gate, top 1 | 0.11206 | 0.07526 | 0.57872 | 0.28944 |
| Support gate, top 3 | 0.12524 | 0.08495 | 0.57842 | 0.29138 |
| Support gate, top 5 | 0.12524 | 0.08495 | 0.58009 | 0.28969 |

The gate slightly improves broad AUROC in some variants while worsening the
sparse top-49 decision boundary. Fold choices are unstable: most use live at low
support and ComplementNB or RBF at high support, but one fold instead chooses
radius-3/Chemprop fusion or the original ensemble for low support. No gate is a
submission candidate.

## Validation-distribution mismatch

| Evaluation population | Mean support | Median support | Below 0.3 | Below 0.4 |
|---|---:|---:|---:|---:|
| UMAP OOF | 0.307 | 0.299 | 51.0% | 93.4% |
| Blind test vs full D6 train | 0.464 | 0.471 | 12.5% | 42.0% |

The UMAP groups deliberately remove close structural neighborhoods, whereas most
test molecules have materially closer labeled neighbors. This explains why the
high-support OOF strata are too small to learn a reliable router and is consistent
with the live D6 leaderboard MCC substantially exceeding its projected UMAP-fold
MCC.

The actionable use of the train/test fingerprint analysis is therefore a second,
prespecified validation scheme whose validation-to-training support histogram
matches the blind-test histogram while retaining group controls. UMAP should
remain the OOD stress test. Candidate adoption should require reasonable behavior
on both, with the support-matched split used for estimating likely leaderboard
ordering. This is preferable to test-frequency bit correction or further gate
tuning on the current five folds.
