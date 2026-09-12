# Frozen external embeddings for TDI classification

Three frozen representation views were compared with TabPFN2.6 on the original
five UMAP folds. Only the two scored endpoints, `CYP2D6_is_TDI` and
`CYP3A4_is_TDI`, are included. No representation fitted to challenge TDI or
challenge multitask labels is reused.

## Inputs and fitted models

Prepared canonical labels/splits come from classification run `269e598d`:
1,497 D6 labels (324 positive) and3,584 3A4 labels (764 positive). Each endpoint's
test set has750 compounds. Its original split identifier is retained in all
probability outputs.

| View | Frozen source | Original dimensions |
|---|---|---:|
| External classification encoder | `a3ea2ec3` | 600 |
| External mixed-objective encoder | `2c0a0953` | 600 |
| Original CheMeleon | `e9576045` | 2,048 |

For each outer fold, StandardScaler and PCA256 are fitted only on that fold's
training rows. TabPFN uses the saved local2.6 classifier checkpoint, four
estimators and seed42+fold. No hyperparameter search is performed. Run
`20260912T115242271820Z-51023dac` completed30 outer fits, saved each fitted
classifier state, PCA/scaler, membership and all probabilities. The installed
library's deprecation warning mentions its default8 estimators, but source
inspection confirms an explicitly supplied4 is never overridden.

## Threshold protocols

Fixed0.5 is the primary untuned classification threshold. The initial run also
records exploratory thresholds selected from other folds' OOF probabilities;
because those base models have overlapping training folds, this is not fully
nested threshold selection.

Separate run `20260912T115749418593Z-397ddee7` fits30 additional classifiers for
an independent threshold comparison. For outer fold `f`, validation fold
`(f+1)%5` selects the threshold, and both folds are excluded from that inner
classifier's training and PCA. It chooses maximum inner MCC over the fixed
0.01–0.99 grid, breaking ties toward0.5, then toward the smaller threshold.
The chosen threshold is applied to the unchanged saved outer probabilities.

All-fold OOF thresholds saved for test use are deployment choices, not additional
OOF accuracy estimates. Test probabilities are the arithmetic mean of five fold
classifiers. Nothing has been submitted.

## Outer probability comparison

| View | D6 AUROC | D6 MCC at0.5 | 3A4 AUROC | 3A4 MCC at0.5 |
|---|---:|---:|---:|---:|
| External classification | 0.5398 | 0.0000 | 0.6944 | 0.1490 |
| External mixed | 0.5607 | 0.0000 | 0.6849 | 0.1142 |
| Original CheMeleon | 0.5395 | 0.0000 | 0.7120 | 0.1713 |

All D6 predictions fall below0.5, so the fixed threshold produces only negative
classes. Its AUROC also remains weak; changing the threshold cannot itself fix
the ordering. CheMeleon is the strongest of these three views for3A4 ranking.
The independent inner-threshold stage completed all30 fits. Its results are:

| View | D6 AP | D6 inner-threshold MCC | 3A4 AP | 3A4 inner-threshold MCC |
|---|---:|---:|---:|---:|
| External classification | 0.2395 | 0.0082 | 0.3619 | 0.1862 |
| External mixed | 0.2588 | 0.0371 | 0.3541 | 0.2133 |
| Original CheMeleon | 0.2486 | 0.0348 | 0.3955 | 0.2241 |

Positive prevalence is0.2164 forD6 and0.2132 for3A4, giving context for average
precision. D6 remains weak after honest threshold selection. For3A4, CheMeleon's
MCC improves from0.1713 at0.5 to0.2241 with independent inner threshold selection.
The exploratory other-OOF-fold threshold result0.2709 was higher; it should not
replace the independent result when comparing against the LGBM inner protocol.

## Completed-output audit

Audit `20260912T120154640294Z-401ba682` passed:15,243 total OOF probability rows,
4,500 mean test probability rows and22,500 fold-test rows, covering all three
views and two endpoints. It verifies canonical labels/IDs/folds, the masks of
all30 outer and30 inner fits, and recomputes all30 selected thresholds from
saved inner validation predictions. The threshold stage leaves every outer
probability unchanged. Recomputed five-fold test means have maximum error zero.

The audit's `outputs/metrics.csv` contains the final AP/AUROC and fixed/inner MCC.
Outer run `51023dac` retains `oof-probabilities.parquet`,
`test-fold-probabilities.parquet` and `test-probabilities.parquet`. Inner run
`397ddee7` retains `oof-inner-threshold.parquet`, inner validation probabilities,
thresholds and `test-fold-inner-threshold.parquet`. Every trained model has a
saved fit state, preprocessing asset and membership table.

The scripts are `tdi_tabpfn_frozen.py`, `tdi_tabpfn_inner.py`, their recorded
launcher, and the probability/membership audit under `scripts/`. Configurations
are under `configs/classification-tabpfn/`. Model and data artifacts remain
outside Git.
