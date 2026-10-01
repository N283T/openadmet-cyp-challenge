# CYP frozen Boltz trunk pooling experiment

The user authorized a thorough exploratory pooling comparison on October 1,
2026 (JST). Expected improvement is deliberately uncertain. The objective is
to measure whether frozen trunk features add useful information, including
complementarity to ordinary molecular features, rather than promote a model.

## Inputs and feature blocks

The canonical post-r1 manifest contains 10,779 complexes across CYP1A2,
CYP2C9, CYP2D6 and CYP3A4, including 750 blind molecules per target. The frozen
Boltz-2 structure trunk includes the full CYP sequence, HEM and the variable
ligand. It retains `s` (384 channels) and `z` (128 channels), without diffusion
coordinates or affinity predictions. Source tensors remain preserved locally.

Pooling uses chain identities and canonical residue indices, with core and
extended pocket definitions fixed before target-label fitting. HEM atom order
is checked against the archived CCD and the extraction parser/tokenizer source.
HEM has 43 heavy-atom tokens; Fe is heavy-atom index 42 (zero-based).

Reusable blocks contain protein, ligand, core/extended pocket and HEM `s`
mean/std/max; both directions of protein/pocket/HEM/Fe/non-Fe-HEM to ligand
`z`; mixed pocket+HEM `z`; and ligand-internal `z`. Pair statistics are
mean/max/std/q10/q90. Reductions use float32 even though source tensors are
float16. Pair embeddings are latent features, not measured distances.

## Fixed comparison

Twelve feature views are declared:

1. Ligand `s` mean and ligand-internal `z` mean/max.
2. PXR basic pooling with whole-protein to ligand `z`.
3. PXR basic pooling with core-pocket to ligand `z`.
4. PXR basic pooling with extended-pocket to ligand `z`.
5. Core pocket and HEM pooled together.
6. Extended pocket and HEM pooled together.
7. Core pocket and HEM as separate blocks.
8. Extended pocket and HEM as separate blocks.
9. Core pocket, Fe and remaining HEM as separate blocks.
10. Core pocket and HEM separate, retaining both `z` directions.
11. Core pocket and HEM separate, adding HEM `s` mean.
12. Core pocket and HEM separate, with mean/std/q10/q90 and ligand `s` std.

PXR basic pooling is protein `s` mean + ligand `s` mean + selected
protein-region-to-ligand `z` mean/max. Omitting HEM from readout does not remove
HEM conditioning from the already-computed trunk.

Each view is fitted alone and concatenated with frozen CheMeleon2048 plus
Mordred1613. The molecular-only input is a matched control for each readout.
Regression uses fixed LGBM L1, LGBM L2, and TabPFN 2.6 with eight ensemble
members. TDI uses LGBM binary and the corresponding TabPFN classifier.
TabPFN uses training-fold-only LGBM gain selection of at most 500 columns.
Imputation medians are learned only from training-fold molecules. All effective
settings, model assets, selected columns and training identities are retained.

All six endpoints use the original five UMAP folds from shared run
`20260909T135130268493Z-edba6fc0`; parent-series and scaffold scores are not
mixed with these OOF results. There are 1,500 regression fits and 500 TDI fits
(including molecular controls), plus fold-local feature-selector fits. These
are fixed recipes, not an Optuna search.

## Evaluation and limitations

Regression reports pooled and per-fold ST-RAE, MAE and Spearman, error
correlation with the matched molecular control, and fixed 10/25/50% prediction
blends. Paired group-bootstrap error differences are conditional on the fitted
models and use the original OOF reference denominator. No multiplicity
correction, nested recipe selection, or confirmation on an independent set is
claimed. Best-of-grid rankings are exploratory, not selected-model estimates.

TDI reports AP, AUROC and MCC at 0.5 as diagnostics. Its extraction covers the
target-labelled pool rather than all 6,146 training molecules. These outputs do
not support deployment call-rate selection. Missing labels remain missing;
unlabelled compounds are not negatives. No candidate or upload is authorized
by this experiment.

## Recorded execution

Every substantive pooling, fitting and reporting command runs through
`pixi run cyp run`, with explicit input aliases. The 840 GB raw inputs use
explicit hardlink snapshots on the same filesystem, record SHA-256 hashes,
verify against the canonical manifest before pooling, and recheck input hashes
at run completion. Hardlinks preserve retained data after a path is unlinked but
share an inode: source tensors must remain immutable. They provide no protection
against in-place edits or disk failure.

The durable suite supervisor is `scripts/launch_cyp_boltz_suite.py`. Its plan,
launch logs, pool/run identities and failures are saved under
`artifacts/boltz-trunk/experiments/20260930T232717Z/`. Actual model console logs
and scientific outputs are in their individual recorded run directories.
Completed runs and informative failures are retained. The suite produces a
recorded English report after all screens finish.
