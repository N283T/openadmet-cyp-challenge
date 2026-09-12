# Independent fixed classification ensemble audit

Audit `classification/20260912T121304359139Z-da4f1375` **passed** for completed ensemble run `b0a05986`. The earlier `a64253d1` failed on mixed integer/float label-column concatenation and is not the audited result.

The independent audit reconstructs each equal-weight mean using grouped long-table aggregation, rather than the source NumPy list average. When a graph output contains both binary `prediction` and continuous `probability`, the continuous probabilities are used. Labels are checked against canonical classification table `269e598d`.

Verified counts and calculations:

- Seven fixed recipes, **35,567 OOF rows** and **35,567 inner rows**: 5,081 original labels per recipe.
- **10,500 deployed test rows**: 750 compounds × two endpoints × seven recipes.
- Exact canonical compound/endpoint/fold membership and original labels across all source families. Inner validation belongs to `(outer_fold + 1) % 5`; its IDs are disjoint from the outer evaluation IDs.
- All **70 inner MCC thresholds** independently recovered using exhaustive contingency matrices and the recorded tie-breaking rule.
- All **168 pooled/fold metric rows** independently reproduced for MCC, AUROC and AP, under fixed-0.5 and inner-MCC policies.
- Each test probability reconstructed from the exact member mean and five-fold mean, within 1e-12. Deployment uses the median of the five inner thresholds; saved boolean calls match those probabilities and thresholds exactly.

## Pooled inner-threshold results

| Fixed recipe | D6 MCC | D6 AUROC | 3A4 MCC | 3A4 AUROC |
|---|---:|---:|---:|---:|
| Classical plain mean | 0.068618 | 0.559487 | 0.274063 | 0.754302 |
| All plain LGBM mean | 0.091415 | 0.565146 | 0.315875 | 0.770557 |
| Balanced LGBM mean | 0.078000 | 0.555453 | 0.338112 | 0.777900 |
| TabPFN mean | 0.069503 | 0.569830 | 0.296691 | 0.743326 |
| Graph mean | -0.016300 | 0.525223 | 0.136132 | 0.642545 |
| Classical + TabPFN mean | 0.075797 | 0.569651 | 0.295722 | 0.767674 |
| All plain + TabPFN + graph mean | 0.096228 | 0.566186 | 0.314955 | 0.768413 |

D6 probability ranking remains weak. The broad mean has the largest observed D6 MCC among these seven fixed means; balanced LGBM has the largest observed 3A4 MCC. Those are screened comparison results, not independently selected winners. The graph-only mean is materially weaker. The audit validates the calculations, not a claim that averaging necessarily improves every constituent model.

Deployment averages saved five-fold models and uses median inner thresholds; it is **not** an all-label model refit. No new model training or submission occurred during this audit.

Outputs: `verification.json`, `reconstruction.csv`, and `pooled-metrics.csv` in audit `da4f1375`. Scripts `classification_fixed_ensemble_audit.py` and `classification_fixed_ensemble_audit_launch.py` pass Ruff.
