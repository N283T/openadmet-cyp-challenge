# TabPFN readout on log2fc-pretrained embeddings

Recorded 2026-09-10 UTC. The first embedding-only readout improves pooled UMAP
regression ST-RAE over the earlier LightGBM pipelines, but lowers classification
MCC at the fixed 0.5 threshold. This compares complete pipelines; it does not
isolate the benefit of auxiliary pretraining.

## Method

- Fixed UMAP five-fold assignment: `20260909T135130268493Z-edba6fc0`.
- Dataset: `20260909T132734119233Z-a7063ccd`.
- Each fold uses its own selected [log2fc-pretrained encoder](log2fc-pretraining.md),
  which excluded that outer fold across all auxiliary heads. The same frozen
  encoder transforms both readout training and validation compounds.
- Inputs are the 2,048-dimensional post-BatchNorm representations, before the
  auxiliary predictor. No predicted-log2fc columns, descriptors, fingerprints,
  external feature selection or dimensionality reduction are added.
- Fit one TabPFN per endpoint and fold: 30 fits in total. TabPFN package 8.5.0,
  explicitly selected model version 2.6, eight estimators, seed 42 + fold, CUDA,
  `fit_preprocessors`, one preprocessing worker, automatic memory saving.
  TabPFN's internal preprocessing remains enabled. `ignore_pretraining_limits=True`
  permits the 2,048 feature inputs. No TabPFN weight fine-tuning or parameter search.
- Fit on endpoint-labeled compounds outside the holdout. Classification uses
  positive-class probabilities and a predeclared 0.5 threshold, without tuning.

## Runs and results

| Fold | Readout run |
|---|---|
| 0 | `20260910T111808494092Z-67df22cf` |
| 1 | `20260910T111905846314Z-8e3ab14b` |
| 2 | `20260910T111959350302Z-5529e407` |
| 3 | `20260910T112052970926Z-dbe145e9` |
| 4 | `20260910T112147200339Z-2bbfe757` |

Audit `20260910T112409569324Z-25cecb67` verifies all 11,606 unique labeled
compound/endpoint pairs occur exactly once in OOF predictions and match source
labels. Regression scores agree with the snapshotted official scorer at 1e-12
tolerance; fold scores reproduce the saved readout scores.

Pooled OOF metrics (not averages of fold scores):

| Regression endpoint | n | ST-RAE ↓ | MAE ↓ | Spearman ↑ |
|---|---:|---:|---:|---:|
| CYP1A2 | 1,412 | 0.833545 | 0.655902 | 0.494348 |
| CYP2C9 | 1,285 | 0.636717 | 0.454392 | 0.641175 |
| CYP2D6 | 1,493 | 0.923544 | 0.597411 | 0.387547 |
| CYP3A4 | 2,335 | 0.530652 | 0.543378 | 0.758748 |

| TDI endpoint | n | MCC ↑ | AUROC ↑ | AP ↑ | Predicted positive rate |
|---|---:|---:|---:|---:|---:|
| CYP2D6 | 1,497 | 0.000000 | 0.566133 | 0.265695 | 0.00% |
| CYP3A4 | 3,584 | 0.225348 | 0.756017 | 0.465194 | 4.44% |

The observed positive prevalence is 21.64% and 21.32%, respectively. CYP2D6
predicts no positives at 0.5; its weak AUROC also limits any threshold-only remedy.
Threshold selection would require training-side validation rather than optimizing
and reporting on these same OOF labels.

| Pipeline, same UMAP split | Macro ST-RAE ↓ | Macro MAE ↓ | Macro MCC ↑ | Macro AUROC ↑ |
|---|---:|---:|---:|---:|
| LightGBM + Mordred | 0.765844 | 0.586536 | 0.138333 | 0.636499 |
| LightGBM + Mordred + Morgan | 0.753365 | 0.580243 | 0.165203 | 0.647892 |
| Pretrained embedding + TabPFN | 0.731115 | 0.562771 | 0.112674 | 0.661075 |

Earlier comparison source: audit `20260909T140626284108Z-b5472b1e`,
[baseline report](analysis/2026-09-09-split-baselines.md). Macros equally average
the four regression or two classification endpoint metrics. This is one fixed
split and seed configuration, without uncertainty estimates or LB confirmation.
A frozen foundation-encoder + TabPFN control is needed before attributing changes
to log2fc pretraining; scratch initialization remains deferred.

## Reproducibility and storage

`scripts/launch_tabpfn_readout.sh` launches the five recorded runs using
`configs/readout/tabpfn-embedding-fold*.json`. Every fold snapshots its encoder
shards, membership, summary, dataset labels/aliases, split assignment and local
TabPFN checkpoint files. Local checkpoint SHA-256:

- Regressor: `c855169116d79d2c957a29107ffc04464555f39aec0eef90c90390f2c2b2008b`
- Classifier: `0578fa56f97e11024e31735aaec2c4e7332584b7730242fbaf6c0bbd0299206a`

Runtime metadata records torch 2.14.0, CUDA 13.0, RTX 5080, installed versions and
effective model constructor parameters. The auxiliary foundation checkpoint
retains its previously documented local-weight provenance limitation.

Each run preserves per-endpoint and combined keyed OOF Parquet, metrics, runtime
metadata and six `.tabpfn_fit` states. These states require the separately saved
base checkpoints; preserve run input paths, or remap them on relocation. A saved
state reload round trip has not been tested. GPU bitwise replay is not guaranteed.

Migration 008 and `scripts/register_readouts.py` register readout lineage, 11,606
predictions and fold metrics after verifying source hashes and fold membership.
The loader's durable record is `20260910T112410751916Z-0da6be2c` under
`artifacts/readout-loads/`. All five fits and the audit completed without failures.
No blinded-test predictions, candidates or submissions were produced.
