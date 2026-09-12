# Latest all-endpoint deployment inventory

Recorded metadata-only inventory: `20260912T114524817386Z-0456dff2`. No training, inference, candidate replacement or upload was performed.

**The latest strong all-endpoint combination is not yet a completed test artifact.** The inventory covers 41 endpoint-weight entries, 30 unique globally selected candidates, and 36 source/recursive child runs. Most required fitted assets already exist, but the global ensemble's test inference and calibrated assembly have not been implemented. This is a multi-family deployment task, not a ready concatenation or a justified substitution of the earlier submitted five-component mean.

## Candidate evidence and readiness

Scores below are observed development OOF results. Selecting the smallest observed value among screened methods does not make it an independently validated optimum.

| Endpoint | Strong current OOF option | ST-RAE | Corresponding test status |
|---|---|---:|---|
| 1A2 | Global126 + interval-loss offset (`4a9dd6ab`) | 0.668194 | Global126 test aggregate absent; requires source-model inference and offset application |
| 2C9 | Global126 + interval-loss offset (`4a9dd6ab`) | 0.442042 | Same missing global test aggregate; old global107 uncalibrated is 0.442189 |
| 2D6 | Native community mean + 9-knot interval-L1 penalty .05 (`d94acc48`) | 0.819776 | Ready: `5d73ca01/calibrated-test.parquet`, matching calibration, 750 rows |
| 3A4 | Global126 + unpenalized 9-knot interval-L1 (`d94acc48`) | 0.383986 | Global test aggregate absent; full-OOF global map also needs fitting/application |

Global174 (`7f98ba4d`) leaves the other three endpoint predictions unchanged within 1e-12 and slightly worsens weight-crossfit D6 to 0.825952. Global126 (`84c28266`) has 1A2/2C9/2D6/3A4 weight-crossfit scores 0.670794/0.443453/0.825436/0.391582. Full-OOF weight-fit scores are a different fitting result and must not be ranked against these cross-fit scores as though equivalent.

The best already-completed **native** calibration test files cover all endpoints (`5d73ca01`), but using those for 1A2/2C9/3A4 would knowingly replace stronger current global OOF models with weaker native members. That is not the user's requested strongest-current all-endpoint deployment.

## Exact selected support

The authoritative complete list is `0456dff2/outputs/selected-support.csv`: every endpoint, mID, full-OOF fitting weight, method name, source run, source script, source model artifact count, test artifact count and nested child run IDs. `support-manifest.json` records original input paths/hashes, model/preprocessor filenames and selection files for all 36 runs.

| Required family / examples | Existing assets | Remaining deployment work |
|---|---|---|
| `m0118` joint geometry/substrate, `4ec73c9f` | Per-endpoint/fold LGBM models and preprocessing joblibs; shared descriptors cover all 6,896 compounds | Exact feature assembly and five-fold inference |
| `m0109` geometry mean/std, `ed3e4d03` | Per-endpoint/fold TabPFN fit states and preprocessors; all-compound frozen/geometry features | Load fit states, reconstruct fold-specific auxiliary inputs, predict test |
| `m0096` fine-tuned primary features, `036b9595` | Per-endpoint/fold LGBM text models; retained embedding/parameter inputs | Exact fold embedding and auxiliary routing, then inference |
| External embedding / adapted embedding LGBM, TabPFN and TabICL | Source result inventories contain saved LGBM, TabPFN fit-state or TabICL joblib assets and relevant PCA preprocessors | Family-specific loading and exact feature-column order; no matching global test aggregate exists |
| Neural direct members `m0042`, `m0071`, `m0077` | Fold checkpoints and exported molecular representations retained | Restore matching direct head and produce canonical test predictions |
| Observation members `m0111/m0112` | Saved neural checkpoints, latent outputs and LGBM readouts | Restore neural feature transformations; generate direct/latent test values and corresponding LGBM inference |
| Nested Caruana `m0032/m0123` | Per-fold selected base LGBM text models and `selected.json`/candidate configurations | Reconstruct each base feature view, predict and apply exact saved inner selections |
| Fixed means `m0079/m0080` | Explicit five-child component recipes recursively inventoried | Assemble exact child predictions; do not replace them with the old submission's different full-data implementation |
| Audit-derived `m0002` parent-retained full LGBM | Upstream five-fold models **are saved**; audit run itself has only OOF | Resolve through source chain below and load exact parent-retained full models |

The initial shallow m0002 inspection saw only audit OOF and raised a possible missing-model refit. The recursive inventory resolves this: `e9f6c457 → 5167f1a0 → e998787e/562f006a/a178d271/66728975/8ec1d5bb`, one run per fold. Each upstream run has the `parent_retained-full.txt` model plus columns/imputer/PCA artifacts. No missing model refit is established by the completed inventory.

## Preconditions and execution mode

1. Reconstruct test auxiliary features from saved pair-excluded predictions. Those source predictions cover 6,896 compounds, while the assembled nested-log2fc training table has only 6,146 compounds × five outer folds (30,730 rows). The current training table cannot simply be joined to test IDs.
2. Restore exact model-family preprocessors, feature ordering, fold identity and neural heads; assemble nested and fixed-mean components using their own saved recipes.
3. Apply the chosen global deployment weights. The inventory enumerates **full-OOF fitted weight support**; its fitting score and the weight-crossfit estimate remain distinct.
4. Fit requested global calibrators on the appropriate OOF predictions and apply to the newly produced global test values. Audit all 750 molecules × four endpoints against canonical IDs and native original SMILES.
5. Establish whether the deployment uses saved five-fold models or new all-label final refits. Most inventoried model assets are CV fits, not all-label refits. The user's earlier full-training preference is not fulfilled merely by averaging these models. No source-family-wide full-data assets matching the latest global recipe were established here.

This bounded inventory supports implementing the missing inference and assembly rather than blindly retraining everything. It does **not** establish that loading every saved estimator will succeed without compatibility work, nor that the result is ready for submission. The parent task owns the authorized submission and recorded cooldown handling.
