# CheMeleon-initialized log2fc pretraining

## Objective and scope

Adapt the CheMeleon message-passing encoder to CYP single-concentration log2fc,
then preserve frozen representations for a later TabPFN readout. This stage does
not train on pIC50 or TDI, tune a downstream model, or submit predictions.
The scratch-initialized comparison is supported by the implementation but deferred.
No claim of improvement over scratch is possible from this experiment alone.

The reference is PXR's `track1_activity/scripts/run_chemprop_pretrain.py` and the
later raw-2048 CheMeleon correction. The PXR recipe contributes target-standardized
MSE, auxiliary-first training, a small regression head, warmup/decay and early
stopping. The network dimensions and pooling follow the supplied CheMeleon weights.

## Auxiliary labels

Preparation run: `20260910T100337258849Z-927eac1c`.
Input: `cyp-challenge-single-concentration-TRAIN.csv`, mapped to the immutable
standardized compound IDs using the dataset's auxiliary aliases.
There are 4,376 compounds, each with four labels: CYP1A2, CYP2C9, CYP2D6 and CYP3A4,
all at `4.95049505e-5 M` (approximately 49.5 µM). There are no repeated
compound/enzyme/concentration groups in this snapshot. Original rows are preserved;
the declared policy for such repeated groups is to average log2fc estimates.
No filtering by FDR, uncertainty weighting, sign reversal or target clipping.

Four heads are weighted equally after per-head standardization. Unlike PXR's
concentration-specific 1:0.5 weighting, there is no basis here to favor one CYP.

## Fold and selection boundaries

Use the existing UMAP assignment `20260909T135130268493Z-edba6fc0` only.
Train one encoder per outer fold, excluding every auxiliary label for compounds
in that fold across all four heads. Auxiliary-only compounds follow the same
assignment. Test structures do not enter training or stopping.

Within the remaining auxiliary compounds, a seeded 90/10 molecule split provides
inner validation for early stopping, as in PXR. This is not an additional random
CV benchmark. Target means and standard deviations use **inner training only**;
neither inner validation nor outer holdout contributes to scaling or gradients.
Constant-target standard deviations fall back to one.

The best inner-validation checkpoint is used for outer prediction and embeddings.
Outer scores are computed only after selection and do not select hyperparameters.
The PXR script scaled before splitting and also exported final-epoch weights; those
behaviors are not carried forward. Explicit train/validation/holdout membership,
scalers, best epoch and checkpoint hash are retained.

## Network and optimization

- CheMeleon checkpoint SHA-256:
  `c376624d3407204e780a0ed13a9ac097cc9bb1c13ef89cdbc633c1715c183651`.
  Its recorded [provenance limitation](chemeleon.md) still applies: local PXR asset,
  pinned by SHA-256, with official Zenodo re-verification unavailable.
- Pretrained BondMessagePassing: width 2,048, depth 6, ReLU, dropout 0; mean pooling.
  All encoder weights are trainable during this auxiliary stage.
- Trained BatchNorm on the pooled representation, then Linear(2048,256), ReLU,
  dropout 0.1 and Linear(256,4). The exported fingerprint is the 2,048-dimensional
  representation after BatchNorm; it does not pass through the predictor head.
- Adam, batch size approximately 64 (balanced batches avoid singleton BatchNorm).
- Initial LR `1.364559692954765e-5`, peak LR `1.364559692954765e-4`, three warmup
  epochs, exponential decay toward `1.364559692954765e-6` over the configured horizon.
  These rates are 10× lower than the prior scratch recipe to adapt pretrained weights.
- Maximum 150 epochs; patience 15; no sweep. Seed 42, float32, deterministic Torch,
  TF32 disabled. Actual per-epoch learning rates and losses are saved.

## Artifacts and use

Launch the fixed five-fold campaign with `bash scripts/launch_log2fc.sh` from the
repository root. Each fold is an independent `pixi run cyp run` with its own config,
source snapshot, declared input copies and immutable completed artifacts.

- `membership.parquet`: auxiliary compound IDs and their optimization/selection roles.
- `metadata.json`: configuration, head order, scalers, architecture and weight provenance.
- `history.jsonl`: append-only epoch losses, learning rates, timing and best epoch.
  Saved epoch indices are zero-based (index 5 means the sixth epoch).
- `best.pt`: selected weights and metadata. `last.pt`: final training/optimizer/RNG
  state for recovery, distinct from the selected best weights. These files can be
  updated while the run is active; completed runs are never overwritten.
- `outer_predictions.parquet`: original-scale log2fc predictions, compound IDs,
  auxiliary endpoints, outer fold and split ID. No in-sample predicted log2fc is
  mislabeled as OOF.
- `embeddings-*.parquet`: selected model's representations for all 6,896 compounds,
  with compound ID and outer fold. Use the encoder associated with the readout's
  outer fold; do not average encoders trained on that fold's held-out compounds.
- `summary.json`: selected epoch, outer auxiliary scores and feature/checkpoint metadata.

The vectors are emitted by float32 inference and serialized as Parquet list values.
These are fold-specific supervised features, not the globally reusable frozen
foundation fingerprints registered in `compound_embeddings`.

PostgreSQL `pretraining_sets`, `pretraining_membership`, `auxiliary_predictions`
and `metrics` store structured provenance, membership and scores. Model/embedding
files remain in the recorded run folders. `scripts/register_pretraining.py` checks
hashes, split lineage, holdout exclusion and best-checkpoint selection before loading.
`scripts/audit_pretraining.py` checks full OOF label coverage and computes pooled
auxiliary metrics; each encoder is evaluated exactly once per held-out compound.

## Deferred downstream work

A TabPFN readout still needs a recorded experiment. The auxiliary OOF scores do not
establish pIC50/TDI performance. Predicted log2fc used on readout training rows will
need inner cross-fitting within that readout's outer-training data; simply joining
these global OOF predictions can leak the readout outer fold through other encoders.
No scratch training, inner cross-fitting for TabPFN, all-data production refit or
competition submission is included in this stage.

## First campaign — 2026-09-10 UTC

| Outer fold | Run | Selected epoch (1-based) | Epochs trained |
|---|---|---:|---:|
| 0 | `20260910T100719601351Z-d8f745ca` | 6 | 21 |
| 1 | `20260910T100846602294Z-f7a00ae9` | 6 | 21 |
| 2 | `20260910T101011080180Z-d3607e4f` | 6 | 21 |
| 3 | `20260910T101136302254Z-e2a36516` | 6 | 21 |
| 4 | `20260910T101301049877Z-5d886b77` | 5 | 20 |

All five runs completed. Selected checkpoints and 6,896 compound embeddings per
fold are in each run's `outputs/`. The rapid decrease in training loss continued
after inner-validation loss stopped improving; retain selected best checkpoints,
not final-epoch weights, for downstream work.

Pooled OOF audit: `20260910T101445700767Z-388556ae`. All 17,504 compound/head pairs
are present exactly once with source labels matching within float32 tolerance.
All fold memberships and best-checkpoint hashes were verified during DB registration.

| log2fc head | n | RMSE | MAE | Spearman |
|---|---:|---:|---:|---:|
| CYP1A2 | 4,376 | 0.6090 | 0.4416 | 0.5797 |
| CYP2C9 | 4,376 | 0.3449 | 0.2718 | 0.7093 |
| CYP2D6 | 4,376 | 0.9228 | 0.7123 | 0.5323 |
| CYP3A4 | 4,376 | 1.0048 | 0.7897 | 0.7053 |

These are auxiliary regression results only. There is no scratch comparison or
pIC50/TDI readout result yet. See [RL-0005](research/2026-09.md#rl-0005--chemeleon-initialized-log2fc-pretraining).
