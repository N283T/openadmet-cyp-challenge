# Honest inner predictions for the historical TabPFN stack

The historical strong classification candidate is the original `stack_mordred`
recipe: PCA256 of frozen CheMeleon, four nested predicted log2fc features and
Mordred, read by TabPFN2.6 with eight estimators. Its original OOF artifact is
`243a1e23`; five-fold test predictions and their mean are in `343eb5ec`.

These existing predictions can remain unchanged while a new inner layer supplies
the probabilities needed for Caruana selection. Simply using the old saved
features for that layer would be invalid: their PCA includes the designated
inner fold, and some auxiliary log2fc encoders also saw that fold's log2fc labels.

For each outer fold `o`, the designated inner fold is `i=(o+1)%5`. Readout training
uses the three other folds. A training row in fold `r` receives log2fc predictions
from an encoder excluding **all of o,i,r**. Inner-fold rows receive the mean of
the three corresponding encoders, each excluding o and i. Ten unique triples
cover all five outer/inner combinations. The original pair-exclusion features
do not provide these quantities, so ten auxiliary encoders are newly trained.

The auxiliary trainer, initialization, target normalization, early-stopping rule
and hyperparameters are the original recipe from `540ba590`, with only the
excluded-fold set expanded. In particular, random stopping-validation rows are
drawn only from the two eligible auxiliary folds. No outer or designated inner
log2fc label enters auxiliary training, stopping or target scaling.

Each classifier then fits unscaled PCA256 only on its three eligible training
folds, using the historical seed42+outer. Mordred's constant/empty-column removal
and median imputation also use only those training rows. TabPFN settings remain
eight estimators and seed42+outer. A single pair-excluded encoder used in-sample
on readout training would alter the original stacking recipe; that shortcut is
not used here.

The run exports normalized `historical_stack_mordred` tables for original outer
OOF, new inner probabilities and original fold-test probabilities. It verifies
canonical labels/folds for all5,081 old classification OOF rows and reconstructs
the old1,500 test means from7,500 fold predictions. Memberships, preprocessing
assets and classifier fit states are retained. This is a nested reconstruction
of the same recipe with less training data; it does not replace or refit the
historical outer/test model predictions.

Scripts are `classification_stack_aux_launch.py`, `classification_stack_inner.py`
and `classification_stack_inner_launch.py`; configurations and the completed
auxiliary-run registry are under `configs/classification-stack-inner/`. Each
training and export stage declares and snapshots its inputs through `cyp run`.
## Completed artifacts

All ten triple-exclusion auxiliaries completed, with run suffixes `3c5d2ee2`,
`2988ffaf`, `f0c8012a`, `59bed70f`, `24f7c55a`, `863ac716`, `57064499`,
`5e5d9dd2`, `bf7515d4` and `7ef863e3`. The registry records their exact excluded
folds and full run IDs.

Classification run **`20260912T132823500863Z-130ac6af`** completed all ten
inner readouts and normalized the historical predictions. Its outputs are:

| File | Rows | Meaning |
|---|---:|---|
| `outer-oof.parquet` | 5,081 | Unchanged historical OOF probabilities |
| `inner-probabilities.parquet` | 5,081 | New independently excluded inner predictions |
| `fold-test-probabilities.parquet` | 7,500 | Unchanged historical five-fold test predictions |
| `test-probabilities.parquet` | 1,500 | Unchanged historical test means |

All tables use method `historical_stack_mordred` and numeric `probability`.
Inner tables explicitly carry `outer_fold` and physical `inner_fold`; the
inherited `fold` column is the physical inner fold. The Caruana reader should
index these selection predictions by `outer_fold`, as for the other candidates.

The original test mean reconstruction passed at tolerance1e-12, and canonical
outer labels/folds were matched exactly before fitting. Expected NumPy warnings
about all-missing Mordred columns occurred while identifying columns to drop;
the actual fitted feature matrices passed finite-value assertions. Classifier
fit states, preprocessing assets, membership tables, auxiliary provenance and
per-inner-fold metrics are preserved. No old probability value or submission
was changed by this reconstruction.

Independent audit **`68fe5a1c`** passed all ten auxiliary exclusion masks,
training-only target scalers, selected epochs and checkpoint hashes; all ten
readout memberships, PCA fits and descriptor medians; all thirty feature-lineage
entries; and all 5,081 inner labels and recomputed AUROCs. It also confirmed
that historical outer and test probabilities remained exactly unchanged.
See [the independent provenance report](2026-09-12-classification-historical-stack-provenance.md).
