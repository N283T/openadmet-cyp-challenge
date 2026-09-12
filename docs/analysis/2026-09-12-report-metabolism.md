# External substrate supervision as CYP regression features

## Question and implementation

This evaluates report idea 4 using observed substrate/non-substrate labels, rather than heuristic reactivity or inferred atom-level sites. The TDC Carbon-Mangels datasets cover CYP2C9, CYP2D6 and CYP3A4. There is no CYP1A2 substrate head in this source. All three predicted probabilities can nevertheless serve as regression auxiliaries for all four challenge endpoints.

Source snapshot: `20260912T092811159257Z-d5082036`. The official TDC metadata maps the three tasks to Harvard Dataverse file IDs 4259584, 4259578 and 4259581. The downloaded Dataverse dataset metadata explicitly includes all three files and identifies the license as CC0 1.0. Raw TSVs preserve Drug_ID, Drug and binary Y. TDC attributes the task to Carbon-Mangels & Hutter (2011); this experiment uses the benchmark labels without relabeling inhibition data as substrate data.

Primary references:
- [TDC metadata](https://github.com/mims-harvard/TDC/blob/main/tdc/metadata.py)
- [TDC deposited data and license](https://doi.org/10.7910/DVN/21LKWG)
- [ATTNSOM reference implementation](https://github.com/dmis-lab/ATTNSOM): a distinct, untested atom-level option. This experiment does not train or claim to reproduce ATTNSOM.

## External model preparation

Feature run: `20260912T092919346028Z-bc6e15fc`.

Structures use the project ChEMBL parent standardization. Every connectivity block found among all 6,896 challenge molecules is excluded before external model training. Conflicting labels within an isoform/connectivity group are quarantined and agreeing duplicates collapsed. Raw rows, selection records, curated identities, external CV membership, fitted classifiers and features are retained.

| Isoform | Raw records | Challenge-overlap records removed | Retained | Positive |
|---|---:|---:|---:|---:|
| CYP2C9 | 669 | 147 | 521 | 110 |
| CYP2D6 | 667 | 146 | 520 | 131 |
| CYP3A4 | 670 | 147 | 522 | 271 |

Morgan radius 2 / 2,048 bits feed two fixed classifiers: logistic regression (C=1) and random forest (500 trees, minimum leaf 3). Three-fold scaffold-grouped external CV is diagnostic only; no external or challenge validation score selects hyperparameters. Models are then fitted to all retained external records and predict six columns for every challenge molecule.

| Isoform | Logistic ROC-AUC | Forest ROC-AUC | Logistic AP | Forest AP |
|---|---:|---:|---:|---:|
| CYP2C9 | 0.6762 | 0.6477 | 0.3379 | 0.3039 |
| CYP2D6 | 0.7364 | 0.7276 | 0.4851 | 0.5305 |
| CYP3A4 | 0.6441 | 0.6680 | 0.6454 | 0.6563 |

These small datasets contain a learnable substrate signal. That alone does not demonstrate usefulness for pIC50 regression.

## Matched challenge experiment

Readout run: `20260912T092952805447Z-28db8483`.

Control is frozen Jeremy medium 600-d embedding + four new nested log2fc predictions, fixed L1 LightGBM (300 trees, 15 leaves, minimum leaf 20, learning rate 0.05, L2=1). It must reproduce `1e7001de` exactly. Two prespecified additions are the three logistic probabilities and all six probabilities. All arms use the same parent CV, anchors and 6,384 validation labels. No challenge molecules are removed. Each descriptor has a missing-value flag; flags are constant zero here.

Audit run: `20260912T093208821534Z-e9701f47`. Control reproduction is exact and all keys, confidence intervals and fold membership checks passed.

| Endpoint | Control | +3 logistic substrate scores | +6 substrate scores |
|---|---:|---:|---:|
| CYP1A2 | 0.691998 | 0.690803 | 0.689724 |
| CYP2C9 | 0.494360 | 0.491880 | 0.489097 |
| CYP2D6 | 0.842287 | 0.842690 | 0.842759 |
| CYP3A4 | 0.442934 | 0.438763 | 0.442862 |
| Macro | 0.617895 | 0.616034 | 0.616111 |

Metric: official ST-RAE, lower is better. Small gains on some endpoints do not improve the focus ensemble: adding the six-score arm at a prespecified 25% weight changes the ae562a59 macro from 0.586931 to 0.589362 (worse). Its 1A2 changes 0.670809 → 0.670938 and 2D6 0.822681 → 0.823078. The 50% addition is also worse. This audit compares against ae562a59 and the submitted 1e8d2643, not a reoptimized global Caruana ensemble.

Decision: retain as a measured small/negative feature experiment, not a promoted replacement. The external substrate task is learned, but these six probabilities do not resolve the 1A2/2D6 errors. This does not test richer substrate encoder pretraining or atom-level SoM supervision. The OOFs are available for the shared candidate pool; no all-data regression refit or submission was made.

## Reproduction and artifacts

- `scripts/report_metabolism_fetch.py`: immutable source snapshot.
- `scripts/report_metabolism_features.py`: curation, external diagnostics, full external classifiers and all-compound scores.
- `scripts/report_metabolism_launch.py`: declared-input paired readout and audit launches.
- `scripts/report_metabolism_readout.py`, `scripts/report_metabolism_audit.py`: fixed controls, saved models and keyed official scoring.
- `configs/report_metabolism/`: resolved experiment configurations.
- Feature file: `artifacts/runs/shared/20260912T092919346028Z-bc6e15fc/outputs/substrate-features.parquet`.
- Regression OOF: `artifacts/runs/regression/20260912T092952805447Z-28db8483/outputs/substrate-oof.parquet`, method `substrate_logistic` and `substrate_six`.
- Audit metrics: `artifacts/runs/regression/20260912T093208821534Z-e9701f47/outputs/metrics.csv`.

All jobs completed; project Ruff checks for the five new scripts passed. No dependency or infrastructure changes were required.
