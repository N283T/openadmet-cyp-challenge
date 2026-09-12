# Continuous TDI auxiliaries for direct pIC50 regression

## Design and source interpretation

This implements the paired-assay suggestion in the data-scarcity research report.
The official [underway snapshot](../../references/official/cyp/underway/20260907T135802Z/snapshot.md)
identifies the direct arm as the NADPH-free preincubation comparator and the TDI
arm as active preincubation with NADPH. These are separate labels, not pooled
potency observations. `TDI_minus_direct` is the measured difference in pIC50;
it is not inferred for unpaired compounds, and no default TDI booleans are used.

Preparation run: `20260912T092802562244Z-dbb38755`.
Raw `TRAIN_TDI.csv` rows are matched by both molecule name and original SMILES to
the immutable standardized aliases. Every mapping succeeds. The direct-arm
values present in this file exactly reproduce the canonical direct measurements
within absolute tolerance 1e-10. No duplicate compound/head keys survive.

| CYP | TDI labels | Paired differences | TDI-only labels | Direct/TDI Pearson |
|---|---:|---:|---:|---:|
| 1A2 | 1,413 | 1,412 | 1 | 0.9863 |
| 2C9 | 1,285 | 1,285 | 0 | 0.9623 |
| 2D6 | 1,497 | 1,493 | 4 | 0.8964 |
| 3A4 | 3,583 | 2,334 | 1,249 | 0.9471 |

The additional measurements chiefly add condition information for the focus
endpoints, rather than many previously unlabeled focus compounds. There are
7,778 TDI values and 6,524 measured differences in total.

## Matched training protocol

All three arms start from the same Jeremy medium encoder bridge `7604ffeb`:

- Control: four direct pIC50 and four primary log2fc heads.
- TDI arm: control plus four continuous TDI pIC50 heads.
- TDI arm plus delta: control plus TDI pIC50 and four measured-difference heads.

The exact existing interval-loss control `20260912T072355846776Z-2e678154` is
reused. Its saved training source and the treatment training source have identical
SHA256 `8ce782bc32f5acc4c91510ff7e9cb343f1a6cecf41c64e15c74de1d6506b3ecc`.
Parameters are identical except added outputs and their losses: maximum 160 epochs,
patience 20, batch 64, seed 42 plus fold, encoder/head learning rates 1e-5/1e-4,
weight decay 1e-4, dropout 0.1, direct-head weights 2, primary weights 1, and new
auxiliary weights 1. Inner-training normalization is per head. Direct pIC50 uses
the official interval distance, auxiliary outputs use MSE. Checkpoints minimize
the mean four-endpoint direct-pIC50 ST-RAE on grouped inner validation.

All labels of outer-holdout compounds, including TDI and derived differences,
are excluded from fitting and checkpoint selection. The fixed `a3923255`
parent assignments and always-training anchors are retained. Adding outputs
changes head initialization and can expand eligible training compounds, as
expected for this auxiliary-task intervention; this is not a paired-random-head
ablation. The standardized differences give equal loss weight even when their
raw variance is small.

Each selected encoder produces fold-specific 600-dimensional embeddings. All
readouts use the same fixed LGBM recipe and the improved nested primary features
`d272774c`. The control's new readout is `20260912T092824781073Z-fef7788e`.
Direct neural predictions and LGBM readouts both retain 6,384 keyed OOF rows.

## Results

| LGBM representation | 1A2 | 2C9 | 2D6 | 3A4 | Macro |
|---|---:|---:|---:|---:|---:|
| Control | 0.692959 | 0.470690 | 0.837797 | 0.430793 | 0.608059 |
| Add TDI arm | 0.690035 | 0.482349 | 0.835413 | 0.432441 | 0.610060 |
| Add TDI arm and paired delta | 0.690031 | 0.483367 | 0.835201 | 0.434798 | 0.610849 |

Direct neural predictions show the same small focus improvement without an
overall gain: control/arm/arm+delta macro 0.629855/0.630811/0.631020;
1A2 0.735342/0.732796/0.731923; 2D6 0.859393/0.855388/0.857389.

For fixed 25% additions to the existing six-component plus focused-2D6 candidate
`ae562a59`, TDI-arm LGBM changes 1A2 from 0.670809 to 0.670078 and 2D6 from
0.822681 to 0.821164. Arm+delta LGBM gives 1A2 0.670890 and 2D6 0.820755.
However, all-endpoint macro worsens to 0.587466 and 0.587856, respectively,
versus 0.586931. The matched control's 25% LGBM addition scores 0.587337.
Neural direct additions slightly improve macro but the original control is better
than either new auxiliary treatment (0.584411 vs 0.585400/0.585176). This is
existing component diversity, not evidence for a TDI-task improvement.

The TDI auxiliary tasks therefore provide a small focus-endpoint effect and no
clear overall improvement under this protocol. Retain the four new OOF candidates
for the expanded ensemble inventory; do not replace the main model automatically.
No calibration, full-data fitting, or submission was performed.

## Run inventory and validation

| Stage | Run |
|---|---|
| TDI preparation | `20260912T092802562244Z-dbb38755` |
| Reused exact control fit | `20260912T072355846776Z-2e678154` |
| Control LGBM with new nested primary | `20260912T092824781073Z-fef7788e` |
| TDI arm fit | `20260912T092823518869Z-8e57e8e7` |
| TDI arm LGBM | `20260912T093415808252Z-708fd656` |
| TDI arm + delta fit | `20260912T092945291275Z-4fe962be` |
| TDI arm + delta LGBM | `20260912T093615055236Z-2a25a497` |
| Successful independent comparison | `20260912T093755925162Z-ad74f080` |

The comparison verifies all six 6,384-row OOF files against the existing candidate
by compound, endpoint, fold, split, truth, and confidence bounds. All 15 training
membership files pass outer-compound exclusion, anchor, and inner-group checks.
Fold-specific best/last checkpoints, histories, embeddings, memberships, neural
predictions and all 60 LGBM models across control/treatments remain saved.
`ad74f080` contains `metrics.csv`, `macro.csv`, `audit.json`, and all prespecified
blends. No performance-selected blend was registered as a production candidate.

The first comparison run `20260912T093728861173Z-390622f6` failed while concatenating
prediction outputs because neural/readout column ordering differed. Model training
and keyed score checks had succeeded. The corrected audit explicitly selects a
common output schema and completed as `ad74f080`; the failed run is retained.
