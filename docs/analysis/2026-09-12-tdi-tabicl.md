# TDI classification: matched TabICL readouts

This comparison replaces the readout while retaining the original five UMAP
folds, endpoint labels and frozen representations. The endpoints are
`CYP2D6_is_TDI` (1,497 labels) and `CYP3A4_is_TDI` (3,584 labels).

Run `7884042f` evaluates external classification-pretrained embeddings
(`a3ea2ec3`), mixed external embeddings (`2c0a0953`) and the original frozen
CheMeleon embeddings (`e9576045`). Each endpoint and outer fold fits its own
StandardScaler and PCA256 on training rows only, exactly as the previous frozen
TabPFN comparison (`51023dac`, inner stage `397ddee7`). PCA uses seed42;
classifiers use seed42+outer. The new readout is TabICL2.2.0 with eight
estimators and the classifier-v2-20260212 checkpoint, explicitly snapshotted
from recorded download `8c08c1de`; automatic downloads are disabled during fits.
The previous TabPFN comparison used four estimators, so this first comparison
includes that estimator-count difference.

For each outer fold, an additional inner model excludes both the outer fold
and designated next fold `(outer+1)%5`; the next fold alone selects the MCC
threshold on the same 0.01–0.99 grid as TabPFN, with ties resolved toward 0.5,
then the lower threshold. That threshold is applied to the separately trained
outer model's validation probabilities. Both fixed-0.5 and honest inner-selected
MCC are reported, together with AUROC and average precision. These inner
probabilities also support subsequent Caruana selection.

A separate matched eight-estimator TabICL versus TabPFN2.6 comparison uses the
actual frozen Jeremy600 representation (`c060d197`). Both backends use identical
train-only StandardScaler/PCA256, outer/inner folds, estimator count and seeds.
No challenge-finetuned representation is reused here.

Every fitted estimator retains its training context; TabICL assets reference
the snapshotted checkpoint without duplicating its weights. Per-fit preprocessing
assets and molecular membership tables are saved. All probability tables carry
compound ID, endpoint, method and split/fold metadata. Inner `fold` is the
physical validation fold; `outer_fold` identifies the selection task.

Scripts: `tdi_tabicl_frozen.py`, `tdi_tabicl_launch.py`,
`tdi_jeremy_readout.py`, and `tdi_jeremy_launch.py`. All fits run through the
experiment recorder with explicit configuration and input snapshots. No
submission is made by this comparison.

## Frozen-view results

Run `7884042f` completed all 60 outer/inner fits. The following pooled OOF
values use the fixed outer folds. `MCC inner` applies each fold's cutoff chosen
only from its designated inner fold.

| Representation | Endpoint | AUROC | AP | MCC 0.5 | MCC inner |
|---|---|---:|---:|---:|---:|
| External classification | D6 | 0.548595 | 0.253519 | 0.050660 | -0.000398 |
| External mixed | D6 | 0.575784 | 0.276615 | 0.054295 | 0.020119 |
| CheMeleon | D6 | 0.545617 | 0.252647 | 0.028624 | 0.036826 |
| External classification | A4 | 0.695725 | 0.362137 | 0.114134 | 0.225658 |
| External mixed | A4 | 0.676621 | 0.339710 | 0.085954 | 0.195888 |
| CheMeleon | A4 | 0.738778 | 0.427603 | 0.222757 | 0.240609 |

The earlier same-view TabPFN run used four estimators. TabICL improves A4
CheMeleon AUROC from 0.711984 to 0.738778 and inner MCC from 0.224095 to
0.240609. D6 changes remain small. This comparison does not isolate the learner
from estimator count.

Independent audit `e71978cd` passed all 60 saved fits. It verifies exact equality
of all 6,146 assignments with canonical UMAP run `edba6fc0`, outer and inner
memberships, train-only scaler/PCA statistics, saved TabICL training arrays and
labels, positive-class probability ranges, thresholds, and the fivefold test
means. Failed audit attempts `51608d5d` and `b2c15e4a` used exact or near-exact
float assertions below the saved float32 precision; both remain recorded. The
completed audit uses tolerances no larger than 1e-6 for reconstructed raw-feature
moments and 1e-10 for encoded training arrays.

## Matched Jeremy600 comparison

Run `0e4bcc8d` uses the same frozen 600-dimensional representation, train-only
StandardScaler/PCA256 matrices, folds, seeds and eight estimators for both
learners.

| Learner | Endpoint | AUROC | AP | MCC 0.5 | MCC inner |
|---|---|---:|---:|---:|---:|
| TabICL8 | D6 | 0.598758 | 0.299654 | 0.064714 | 0.022169 |
| TabPFN8 | D6 | 0.580581 | 0.277636 | -0.013588 | 0.068954 |
| TabICL8 | A4 | 0.768924 | 0.462047 | 0.285636 | 0.263135 |
| TabPFN8 | A4 | 0.766157 | 0.464641 | 0.271691 | 0.312712 |

TabICL has higher ranking performance on both endpoints, especially D6. The
honest transferred thresholds favor TabPFN, so the result does not support a
universal learner ordering. Jeremy600 is the strongest tested TabICL
representation on both endpoints by AUROC.

## Historical stack readout replacement

Feature preparation run `665df66e` freezes the original CheMeleon PCA, nested
predicted log2FC and Mordred matrices. Run `181ca91b` changes only the final
eight-estimator classifier to TabICL and saves 20 new inner/outer fits and test
probabilities. Its OOF results are D6 AUROC 0.516294, inner MCC 0.064554, and A4
AUROC 0.763943, inner MCC 0.269580. The old historical TabPFN remains stronger
for A4 at inner MCC 0.365149. TabICL is therefore useful through Jeremy600 and
the expanded ensemble rather than as a direct replacement for that stack.

The completed runs feed the 24-candidate Caruana comparison `97df4c71`; its
independent replay `d8ca5f7b` passed. Target-wise packaging run `d30a0365`
retains the scored D6 vector and uses the A4 AUROC-Caruana deployment prediction.
Candidate `f0deac96` and prepared submission `82d53182` contain 43 D6 and 361 A4
positive calls. Dry-run validation passed with SHA256
`c6cbfed75f43dcce6a333c6e0f614caa8a9f12f19956f1ac494c9c8d9f1e596f`.
It was not uploaded.
