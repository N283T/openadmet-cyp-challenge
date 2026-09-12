# D6 training sanity: regularization and shuffled labels

Run `classification/20260912T130504137764Z-465afc7e` completed10 fixed model fits. This is a diagnostic comparison, not deployment selection.

All models use exactly the canonical1,497 D6 calls and original UMAP folds, Morgan radius2/2048-bit fingerprints concatenated with the frozen Mordred tables declared by `76facc71`. No observed auxiliary assay targets are features. Both recipes use balanced binary LightGBM, learning rate0.03, deterministic seed42+fold, three CPU threads, and no validation-driven stopping or tuning.

- Strong regularization:100 trees,7 leaves, minimum child80, L2 penalty10.
- Label-permutation control: original400 trees,15 leaves, minimum child20, L2 penalty2. Only outer-training labels are shuffled using a separate seed42+fold permutation. Held-out labels remain true. Each training fold's prevalence is preserved exactly.

| Fixed model | Training AUROC range against fitted labels | Pooled true OOF AUROC | OOF AP | OOF Brier |
|---|---:|---:|---:|---:|
| Original400-tree balanced model, saved `76facc71` | See prior training audit | 0.533285 | 0.261179 | 0.187446 |
| Strong regularization | 0.948464–0.967368 | 0.530677 | 0.245606 | 0.228537 |
|400 trees, shuffled training labels | 0.999991–1.000000 | 0.509864 | 0.227036 | 0.190047 |

The original unbalanced400-tree reference scores0.534943 AUROC,0.252959 AP and0.184214 Brier on exactly the same canonical OOF identities. Class prevalence is0.216433.

The shuffled-label model almost perfectly ranks its randomly assigned training targets, while its training predictions scored against the **original** labels have AUROC0.478880–0.515180. Therefore near-perfect training AUROC is achievable through memorization in this feature/model regime; it is not evidence that the learned function captures TDI biology or will transfer to other structural groups.

Reducing tree complexity lowers training AUROC slightly but does not improve held-out discrimination over the original balanced model. These two fixed settings do not resolve the weak D6 generalization. The strongly regularized balanced model also has worse Brier; balanced weighting changes the probability target, so that difference alone should not be interpreted as a pure ranking result.

The shuffled model's pooled held-out AUROC0.509864 is consistent with a weak/no-signal control. This uses only one seeded permutation per fold; it is **not** a permutation significance distribution and cannot establish that the original0.533285 is statistically distinguishable from chance.

Saved artifacts include model files, original and fitted per-molecule training labels, all fold memberships and seeds, per-fold train predictions, complete1,497-row OOF for each recipe, and train/OOF AUROC/AP/Brier metrics. True outer labels are excluded from fitting. No test predictions or submission were produced. Both scripts pass Ruff.
