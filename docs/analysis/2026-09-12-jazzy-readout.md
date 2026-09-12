# Controlled Jazzy descriptor comparison

## Question and fixed design

Can DG-generated Jazzy molecular descriptors add information beyond the existing
external-pretrained representation and the improved primary-screen predictions?
The primary targets of interest are CYP1A2 and CYP2D6; all four direct-inhibition
endpoints remain in the comparison.

The control uses the frozen Jeremy medium 600-dimensional representation
(`20260912T055819206227Z-c060d197`) plus the four improved nested log2fc predictions
(`20260912T082554377951Z-d272774c`). It reproduces the fixed LightGBM readout
`20260912T083010595992Z-1e7001de`: L1, 300 trees, 15 leaves, minimum 20 rows per leaf,
learning rate 0.05, L2 regularization 1, and seed 42 plus outer-fold index.
The retained parent split is `20260911T104703655337Z-a3923255`.

Three prespecified arms share every split, training label, base feature and model
parameter: control; control plus Jazzy descriptor means; control plus Jazzy means
and standard deviations across conformer seeds. No molecules are dropped when
descriptor extraction fails. Jazzy finite/nonconstant feature selection and median
imputation use only the outer training partition, with missingness indicators
appended for retained descriptor columns. The standard deviation describes
seed sensitivity, not a Boltzmann-weighted conformer ensemble. These descriptors
do not claim to model an assay-specific protonation ensemble.

Each arm produces 6,384 keyed OOF predictions and 20 saved fitted models. The run
asserts that the control predictions are bitwise identical to the previous control.
It saves train/validation memberships and each fitted descriptor preprocessor.

The independent audit checks all keys, labels, bounds and memberships, reports
per-endpoint/per-fold and neighbor scores, and compares prespecified 25% and 50%
additions to both the latest unsubmitted candidate
(`20260912T083627802028Z-ae562a59`, 50% combined variant) and the submitted five-way
ensemble (`20260912T072935321335Z-1e8d2643`). The same blends with the control are
included so any gain from model averaging can be distinguished from a Jazzy gain.
These are development comparisons, not selection-adjusted independent estimates.

## Execution

Feature extraction `20260912T085529851955Z-9df14c0c` used Jazzy 0.1.4,
RDKit 2024.3.1 and Kallisto 1.0.10 in the separately locked Jazzy environment.
The standard Jazzy path is **ETDG**, without force-field minimization or a 2D
fallback; seeds were 11, 42 and 89. The six molecular descriptors are `sdc`, `sdx`,
`sa`, `dga`, `dgp` and `dgtot`. Exactly 6,884 of 6,896 molecules succeeded for all
three seeds. The remaining 12 have nonzero net charge and retain missing values
for all three seeds because the standard Jazzy API assumes total charge zero
in its Kallisto calculation; we therefore leave net-charged molecules missing
instead of applying that assumption to them. There were no embedding failures among the
supported molecules. `feature-report.json` preserves coverage and versions.
The complete 6,896-molecule, three-seed extraction took approximately 114 seconds
with eight CPU workers, including run setup and serialization.

The first launch `a1908f3f` failed because the relative interpreter path was not
available in the source snapshot. It remains recorded; the successful extraction
used the absolute project-environment interpreter path.

Implementation: `scripts/jazzy_readout.py`, `scripts/audit_jazzy_readout.py`, and
`scripts/launch_jazzy_readout.py`. Extraction and readout runs are recorded
separately, with all feature and reference inputs snapshotted. No submission is
part of this experiment.

## LightGBM result

Paired fit `20260912T085800716037Z-db5687c1` and independent audit
`20260912T085911904925Z-b8962a03` completed. The control is bitwise identical to
the earlier LGBM OOF on all 6,384 labels. All 60 models and fitted preprocessors
are saved; keyed membership and official-score checks pass.

| OOF ST-RAE | Control | Add means | Add means and standard deviations |
|---|---:|---:|---:|
| 1A2 | 0.691998 | 0.686434 | 0.693365 |
| 2C9 | 0.494360 | 0.492430 | 0.481813 |
| 2D6 | 0.842287 | 0.842709 | 0.834582 |
| 3A4 | 0.442934 | 0.441020 | 0.442689 |
| Macro | 0.617895 | 0.615648 | 0.613112 |
| Neighbor macro | 0.729222 | 0.721243 | 0.717054 |

There is a small, inconsistent endpoint-specific signal: means help 1A2, whereas
mean-plus-variability helps 2D6 and 2C9. Adding the latter at 25% to the latest
candidate worsens macro OOF from 0.586931 to 0.589044. Its 2D6 score improves only
from 0.822681 to 0.821798; the 1A2 mean-only 25% addition improves 0.670809 to
0.670192. These tiny, differently selected endpoint gains do not justify replacing
the candidate. The old submitted reference improves when this readout is added,
but most of that improvement also occurs with the Jazzy-free control and must not
be attributed to Jazzy.

The mixed result motivated a second matched readout using TabPFN 2.6 with eight
estimators, the same 600+4 base features, and the same three descriptor arms.
It uses control `20260912T083032513346Z-568bc526` and a prespecified maximum
prediction discrepancy tolerance of 1e-5 for repeated GPU inference. Exact
observed discrepancy is saved rather than assuming bitwise GPU determinism.

## TabPFN result and decision

Paired TabPFN run `20260912T090040918188Z-695fb6ee` and independent audit
`20260912T090905337416Z-03395b00` completed all 60 fits. The reproduced control
was in fact bitwise identical, with maximum absolute prediction difference 0.
All 6,384 keys, labels, bounds and memberships pass the same verification.

| OOF ST-RAE | Control | Add means | Add means and standard deviations |
|---|---:|---:|---:|
| 1A2 | 0.705017 | 0.703028 | 0.705406 |
| 2C9 | 0.495867 | 0.499524 | 0.494894 |
| 2D6 | 0.874374 | 0.874372 | 0.874545 |
| 3A4 | 0.439395 | 0.435385 | 0.433501 |
| Macro | 0.628663 | 0.628077 | 0.627087 |
| Neighbor macro | 0.714844 | 0.721833 | 0.710916 |

Means again have a small 1A2 benefit, but no meaningful 2D6 gain appears in this
readout. Neither descriptor arm improves the latest candidate when added at 25%
or 50%; the two 25% macro scores are 0.591396 and 0.591504 against 0.586931.

Retain the inexpensive extracted feature asset and the full negative/mixed
comparison. Do not replace the latest candidate or submit a new model on this
evidence. This recipe adds little on the two focus endpoints and does not provide
the larger improvement sought. It does not establish that every geometry-based
feature or every way of using Jazzy is ineffective. No additional tuning was
performed to chase these small development differences. All launched fits are
complete, and both model families reproduce their historical controls exactly.
