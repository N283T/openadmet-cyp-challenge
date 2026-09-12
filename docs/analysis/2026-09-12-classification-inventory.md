# Classification restart: labels, UMAP folds and safe feature sources

Preparation run `classification/20260912T115018743045Z-269e598d` completed without model fitting. Its canonical labels and original shared UMAP assignments are ready for the LGBM and TabPFN screens.

## Targets and official metric

The official classification targets are `CYP2D6_is_TDI` and `CYP3A4_is_TDI`. The primary metric is Matthews correlation coefficient (MCC), aggregated across these two targets. Predictions are boolean calls for all 750 test molecules; only confidently assigned hidden labels contribute to evaluation. AUROC and average precision are useful probability-ranking diagnostics, but do not replace MCC.

The pinned official challenge description explains positive shift labels, inferred positives and assigned negatives. For a measurable direct arm, a shift exceeding twofold is positive. Weak direct-arm compounds may be inferred positive when the TDI arm exceeds pIC50 4.301 or assigned negative when both arms are below 4. These conventions mean a naive sign/threshold of a newly joined potency difference is not an interchangeable label definition. This preparation keeps the original supplied boolean labels exactly and does not fill ambiguous/missing labels or derive replacements.

Official source: `references/official/cyp/challenge/20260907T135802Z/snapshot.md`, included in the preparation run. Canonical measurements: shared `a7063ccd`.

## Class balance and fixed shared split

Use the project-default original UMAP split `20260909T135130268493Z-edba6fc0`, filtering `strategy == "umap"`. It assigns all 6,146 training compounds once to five folds and preserves the shared group IDs. The 750 blind compounds are absent from assignments. Regression's later parent-retained `a3923255` split is a different design and must not be silently substituted.

| Endpoint | Labeled | Positive | Negative | Positive rate |
|---|---:|---:|---:|---:|
| CYP2D6 | 1,497 | 324 | 1,173 | 21.64% |
| CYP3A4 | 3,584 | 764 | 2,820 | 21.32% |

| Fold | D6 n | D6 positive | D6 rate | 3A4 n | 3A4 positive | 3A4 rate |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 287 | 57 | 19.86% | 736 | 131 | 17.80% |
| 1 | 325 | 70 | 21.54% | 655 | 136 | 20.76% |
| 2 | 311 | 84 | 27.01% | 761 | 196 | 25.76% |
| 3 | 327 | 63 | 19.27% | 695 | 145 | 20.86% |
| 4 | 247 | 50 | 20.24% | 737 | 156 | 21.17% |

Both classes are present in every endpoint/fold. Fold 2 has a higher positive fraction in both targets, so report fold-level probability and threshold behavior alongside pooled scores. Do not stratify the targets separately or alter these fixed assignments during the first matched screen.

## Earlier results

- Original descriptor/Morgan LGBM experiments used the same UMAP assignments (`eb0dd9c8`, metric audit `b5472b1e`). UMAP Mordred+Morgan macro MCC was about 0.1652; D6 AUROC about 0.536 and MCC about 0.094.
- The later historical shared OOF audit `243a1e23` reported D6 MCC@0.5 = 0, AUROC = 0.542860, AP = 0.237669; 3A4 MCC@0.5 = 0.258678, AUROC = 0.793282, AP = 0.492060.
- Classification threshold fit `c95bbbfc` selected D6 threshold 0.135978 (selection-only MCC 0.086555) and 3A4 threshold 0.319378 (0.383103). It predicted 686/750 D6 positives and 249/750 3A4 positives. Those same-OOF selection scores are not held-out threshold performance.
- The historical accepted submission `f2dbb957`, documented in `docs/submission-baseline.md`, had Sep10 LB D6 MCC 0.001644 (precision 0.070760, recall 0.923001) and 3A4 MCC 0.430653 (precision 0.562266, recall 0.662249). Macro MCC was 0.216149. These are the historical snapshot values, not newly fetched live ranks.

D6 was weak before and after threshold adjustment; aggressive positive calls did not repair weak probability ranking. Keep raw probability outputs and evaluate thresholds trained only on other OOF folds during the new screen.

## Feature eligibility

| Source | Initial use | Condition |
|---|---|---|
| Morgan fingerprints, Mordred, static ligand geometry/Jazzy descriptors | Eligible | Fit imputation/scaling/PCA only within each training fold where needed |
| Original frozen CheMeleon | Eligible | No challenge fine-tuning |
| Frozen public/external-pretrained embeddings and public prediction scores | Eligible | Confirm their weights use no challenge target labels; retain external source provenance |
| Frozen Admetica | Eligible | Use unchanged external weights; challenge-finetuned variants are different |
| Observed TDI pIC50, observed direct pIC50, measured delta, Emax, primary screen values or label-availability masks | Excluded as direct classifier inputs | They are unavailable at blind inference and may determine the classification label |
| Existing regression-finetuned encoders / predicted TDI or delta from `a3923255` folds | Not automatically eligible | Their held-out sets differ from classification UMAP; require verified all-head exclusion under the actual classification folds |
| Newly cross-fitted potency/TDI/primary predictions | Potential later branch | Train auxiliary models with every classification-held-out compound excluded across all heads, including any external aliases; use matching test inference |

This applies even when an input is called an embedding: supervised adaptation to a label-defining assay can leak through the representation. Predicted values require valid fold provenance, not merely a `pred_` column name.

## Prepared artifacts

Under `artifacts/runs/classification/20260912T115018743045Z-269e598d/outputs/`:

- `classification-labels.parquet`: 5,081 rows with canonical compound ID, endpoint, original value/bounds, `y_true` Int8, `strategy`, `fold`, `group`, and `split_id`. It already includes fold columns.
- `umap-assignments.parquet`: 6,146 unique training IDs with original UMAP fold/group and explicit split ID.
- `classification-compounds.parquet`: all 6,896 compounds with nullable wide classification labels; do not expose those label columns as features.
- `label-balance.csv`, `fold-balance.csv`, `summary.json`: counts, source policy and checks.

The preparation verifies binary supplied labels, unique compound/endpoint keys, exact counts, full assignment coverage, five folds and train/test disjointness. Script `scripts/classification_inventory_prepare.py` and configuration `configs/classification/umap-inventory.json` retain the procedure. No models or submissions were produced by this task.
