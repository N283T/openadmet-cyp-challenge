# Independent classification LGBM reconstruction audit

Audit `classification/20260912T120959802121Z-8aeb9a25` **passed** for completed screen `76facc71`. The audit performs inference from saved models only; it does not retrain a model or change the source predictions.

The screen contains six feature views in plain/balanced variants: 12 methods, two official endpoints and five original UMAP folds. Each outer fit has a separate provisional classifier whose training excludes both the outer fold and the predefined inner fold `(outer + 1) % 5`.

## Verified

- Exactly **60,972 OOF predictions**: 5,081 canonical labels per method × 12 methods, with exact original labels, IDs and UMAP folds.
- Exactly **90,000 fold-level test probabilities**: 750 canonical test compounds × two endpoints × five folds × 12 methods.
- Training and test IDs are disjoint. Saved inner/outer membership masks agree exactly with the original UMAP assignments. Provisional training row counts exclude both held-out folds.
- All 120 provisional LightGBM boosters were loaded and used to reconstruct their inner-validation probabilities from the declared frozen/static feature inputs.
- Each threshold was independently recovered by exhaustive candidate-threshold contingency matrices. The source implementation uses cumulative sorting; this audit uses a different computation and confirms all **120 MCC-optimal thresholds**, including the fixed tie-breaking rule (nearest 0.5, then higher threshold).
- All 120 final outer boosters were loaded and used to recompute both their outer OOF probabilities and all 750 test probabilities. Every result agrees with the source within 1e-12; every saved threshold and split ID also agrees.

Feature assembly independently restores canonical ID order, original Mordred column order, external `emb_` columns, exact Morgan radius 2/2,048-bit fingerprints and float32 nonfinite-to-NaN handling. Feature dimensions are checked against each recorded fit.

## Recovered inner predictions

`8aeb9a25/outputs/inner-predictions.parquet` contains **60,972 rows** with:

`compound_id`, `endpoint`, `inner_fold`, `group`, `y_true`, `outer_fold`, `method`, `prediction`, `threshold`, `split_id`.

For subsequent fixed model ensembles, combine inner probabilities using matching **outer_fold, inner_fold, compound_id and endpoint**, then fit that ensemble's threshold on those inner labels. Do not use outer OOF labels to fit the threshold. These recovered probabilities come from the provisional models, not the final models that included the inner fold during refitting.

Other outputs are `reconstruction.csv` (120 fitted-pair checks and inner MCC) and `verification.json`. Scripts: `classification_lgbm_audit.py` and `classification_lgbm_audit_launch.py`; both pass Ruff. All original matrices, source configuration, 240 boosters, membership, probabilities and metadata were declared as immutable audit inputs. No classification submission was made.
