# Rustuna tuning of the primary-screen auxiliary model

## Question and design

Test whether auxiliary-model hyperparameters improve predicted log2fc and the
downstream pIC50 readout. Rustuna 0.1.0 is installed through pixi; the lockfile adds
only Rustuna and leaves existing model dependencies unchanged.

The search uses 16 trials including the original recipe as trial zero. Rustuna TPE
uses seed 42 and six startup trials. The criterion is the mean squared error of
four standardized log2fc targets on the fixed 261-compound internal validation
set; 2,357 compounds train the model. The fold-pair 0/1 exclusion is retained:
neither fold participates in preprocessing, fitting, stopping or HPO selection.
All trials retain their configuration, membership, scaler, epoch history, best
and final weights inside the recorded campaign run. No outer metrics are computed
by the tuning path.

The candidate recipe is selected once and refitted on pairs 0/1, 0/2, 0/3, 0/4.
These four models feed a fixed full-feature LightGBM readout on untouched fold 0.
The comparison uses identical parent-retained validation molecules, labels and
readout parameters. The primary-data evaluation and pred_log2fc-to-pIC50
correlations are also measured on fold 0. This is an independent outer-fold
evaluation of the tuning decision, not five-fold nested HPO.

## Search space

| Parameter | Search |
|---|---|
| Head peak learning rate | Log-uniform 2e-5 to 1e-3 |
| Encoder/head learning-rate ratio | Log-uniform 0.03 to 1 |
| Weight decay | Log-uniform 1e-6 to 0.05 |
| Head dropout | Uniform 0 to 0.5 |
| Pooled BatchNorm | On/off |
| Schedule | Noam-style warmup/exponential decay or ReduceLROnPlateau |

New trials use AdamW, 80 maximum epochs, early-stopping patience 15, head width
256 and batch size 64. Noam uses three warmup epochs, 0.1 times peak initially,
and 0.01 times peak finally. Plateau starts at the sampled peak and halves LR
after four non-improving epochs (PyTorch patience=3), with a 0.01-times-peak floor.
Encoder and head preserve their learning-rate ratio.

Trial zero retains Adam, 150 maximum epochs, original learning rates and all other
original settings. Its training and validation losses reproduce the saved original
run, including best validation MSE 0.586266040802002 at zero-based epoch 7.
The new controls are optional; original configurations retain their behavior.

## Execution records

- Initial run `20260911T120309477536Z-74878cd7` failed after two training fits:
  incremental summary writing used the append-only JSON writer on an existing
  path. All artifacts remain. The corrected driver appends JSONL records and
  writes the final JSON summary once.
- Corrected search: `20260911T120540987924Z-3f5823d1`.

Source: [Rustuna](https://github.com/optuna/rustuna). Rustuna accelerates the search
implementation; GPU model fitting dominates this experiment's runtime.

## Results and decision

Trial 14 wins: validation MSE 0.539817 versus original 0.586266 (7.92% lower).
Its best checkpoint is epoch 15 (one-based), stopping at epoch 30. Parameters:
head peak LR 0.000479166; encoder/head ratio 0.194756; AdamW weight decay
0.00591952; dropout 0.332311; no BatchNorm; Noam-style schedule. Thus the winning
change is a combination of settings, not evidence that plateau scheduling alone
was necessary. No single-parameter causal attribution is supported by this search.

Refit pair 0/1, 0/2, 0/3, 0/4 runs, respectively:

- `20260911T123506505890Z-15c68a1d`
- `20260911T123637922123Z-692ffa81`
- `20260911T123754944170Z-4d39cd56`
- `20260911T123923896514Z-aa20fce2`

Full-feature LightGBM readout: `20260911T124052504276Z-c22b0df7`.
Paired audit: `20260911T124152993630Z-ee516e03` (48 registered readout metrics).
Every validation key, label, confidence interval and neighbor flag matches the
original fold-zero readout. The figure `tuning.png` was visually inspected.

On all 874 primary-screen compounds in untouched fold 0:

| Enzyme | Original RMSE | Tuned RMSE | Original Spearman | Tuned Spearman |
|---|---:|---:|---:|---:|
| 1A2 | 0.579382 | 0.575349 | 0.589413 | 0.610052 |
| 2C9 | 0.330328 | 0.313795 | 0.714117 | 0.738454 |
| 2D6 | 0.893576 | 0.875488 | 0.523392 | 0.536036 |
| 3A4 | 0.934014 | 0.913086 | 0.733537 | 0.744814 |

Predicted log2fc versus direct pIC50, all nonanchor labeled fold-zero compounds
(all available pIC50 labels, including those without observed primary):

| Enzyme | n | Original Spearman | Tuned Spearman |
|---|---:|---:|---:|
| 1A2 | 283 | -0.407574 | -0.432707 |
| 2C9 | 244 | -0.644339 | -0.660996 |
| 2D6 | 297 | -0.340389 | -0.332895 |
| 3A4 | 513 | -0.725683 | -0.740283 |

Downstream LightGBM results:

| Subset | Original ST-RAE | Tuned ST-RAE | Original MAE | Tuned MAE |
|---|---:|---:|---:|---:|
| All nonanchors, fold 0 | 0.744666 | 0.740726 | 0.524782 | 0.523594 |
| Hit neighbors, fold 0 | 0.929463 | 0.989661 | 0.659707 | 0.678479 |

Per-endpoint all-nonanchor ST-RAE: 1A2 0.905047 → 0.885725; 2C9 0.622268 →
0.626637; 2D6 0.906370 → 0.914864; 3A4 0.544977 → 0.535678. Neighbor endpoint
counts are 21/9/6/30, so the neighbor result is retained with its actual support.

The auxiliary objective improves on untouched primary data for all four enzymes.
The downstream gain is small overall, absent in 2D6, and reversed on this fold's
neighbor subset. Preserve the tuned candidate, but do not replace the five-fold
baseline based on this experiment. No test inference or submission was made.
The experiment supports improving auxiliary fit; it does not resolve the 2D6
pIC50 limitation. Checks and 44 DB-enabled tests passed.
