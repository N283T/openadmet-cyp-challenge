# Native SuperCowPowers Chemprop reproduction

## Scope and implementation

The four published model scopes were trained using the **original Workbench `train_chemprop_fold` function**, pinned to historical 18-head commit `4cf07a8cb4366eb29111fb42396336b20b69ada0`. This replaces the earlier custom 600-dimensional transfer and LGBM experiments as the implementation used for this reproduction. Source snapshot: shared run `97b9c270`; source split implementation: `9388c902`; exact author data: `ae7e37fc`.

`scripts/scp_repro_train.py` loads the archived source module and provides a small read-only Polars-backed frame adapter. The original code constructs the native Chemprop model, standardizes targets using the training partition, creates its loss and optimizer, trains with Lightning, restores the best validation checkpoint, and predicts in original target units. The adapter adds durable membership tables, epoch histories and native prediction arrays. It does not replace the training objective or optimizer.

Verified source settings are scratch initialization, message-passing width 700 and depth 5, norm aggregation, batch normalization, two FFN layers of width 300, dropout 0.1, native MAE, batch size 64, maximum 400 epochs and patience 40. The native Adam/Noam schedule warms up for two epochs from 1e-4 to 1e-3, then decays to 1e-4. Training uses mixed precision; prediction uses source FP32 behavior. Seeds are 42 + fold. CPU threads are limited to two per process and dataloader workers to zero; these are resource settings. No rows, epochs or batches are artificially capped.

The multitask loss uses the source's weighted error numerator and **unweighted observed-label count denominator**. Union18 and challenge8 use inverse-count weights for the four direct endpoints, normalized to mean one, with every auxiliary head weighted 0.3. Counts are computed globally before CV, as in the author scripts. D6isoform6 uses weights `[1, .3, .3, .3, .3, .3]`. D6single1 uses the native unweighted criterion.

## Data and validation

| Recipe | Rows | Heads |
|---|---:|---|
| union18 | 31,670 | Direct4 + log2fc4 + ChEMBL5 + max-response5 |
| challenge8 | 4,905 | Direct4 + log2fc4 |
| D6isoform6 | 23,264 | D6 direct, log2fc, TDI pIC50, **direct-condition Emax**, ChEMBL, max-response |
| D6single1 | 1,493 | D6 direct |

The union has exactly 161,400 observed labels, matching the historical published entry. Author SMILES and label values are preserved. One publicly measured Veith compound shares blind-test connectivity; the public response remains included, matching the author data and the authorized external-data scope. No hidden challenge label is available or used.

The main comparison deliberately keeps challenge parent folds `a3923255` rather than claiming to reproduce the author's scaffold OOF. Whole rows, including public labels, are withheld at held-out challenge connectivity; anchors remain in training. Public-only molecules retain the author's scaffold assignment algorithm: first-appearance scaffold IDs, shuffled group order with seed 42, then greedy least-loaded fold assignment. Thus validation includes public tasks as well as challenge tasks. Checkpoint selection uses native held-fold `val_loss`, as in the author pipeline; these are development OOF predictions, not the previous nested selection design.

A separate D6single1 run uses source-exact scaffold CV over all 1,493 rows. It has a distinct split identifier and can be compared on the same 1,457 nonanchor canonical labels as the parent-CV run.

The source's "100%" prediction policy is the five-fold ensemble over the available corpus. Each test prediction is the mean of five native model predictions. No separate all-row refit or LGBM readout is inserted.

## Run ledger and export correction

| Recipe | Native training run | Canonical prediction export |
|---|---|---|
| D6single1 parent | `fec30715` | `0f300d9d` |
| challenge8 parent | `423d2caf` | `dcf54932` |
| union18 parent | `288cae8b` | `d5c85419` |
| D6isoform6 parent | `e6e93fbf` | `e6fbb180` |
| D6single1 source scaffold | `03b86209` | `ea687d29` |

Initial adapter startup `3d49be88` failed because declared inputs have recorder-generated numeric filename prefixes. It remains recorded. The corrected launcher uses these prefixes.

The first launched configurations used bare target names when joining canonical evaluation labels. Training and all native validation/test arrays were unaffected, but the initial long-form OOF exports were empty. This was detected from the first completed summary. The separate recorded `scp_repro_export.py` reconstructs canonical keyed predictions from the immutable native arrays, checks hard expected counts (6,384 for four endpoints; 1,457 for D6 parent CV; 1,493 for source CV), and verifies parent fold membership. All reported metrics must use the recovered exports, not the original empty long tables. Future launcher configurations and trainer assertions are corrected.

D6single1 parent training completed with best checkpoint epochs 10, 10, 14, 13 and 12 (zero-based). It produced 1,457 keyed OOF predictions and 750 test predictions after export recovery. The public union processes 394 full training batches per fold-0 epoch; D6isoform6 processes 293, confirming actual full-corpus passes.

The published local comparison script uses a plain absolute-error baseline denominator, while this project's previously used scorer soft-thresholds the baseline denominator. The independent metric audit and root evaluation report distinguish those definitions; training always uses the native MAE objective.

Independent data audit `08a9df69` confirms exact source cohort/head counts and zero cross-fold connectivity leakage in all four parent-adapted recipes. The source-scaffold single control preserves the source algorithm's two cross-fold connectivity groups; this diagnostic does not have the parent design's connectivity isolation. Export audits `b1080b38`, `edea60d1` and `6b3569b5` verify parent single, source single and challenge8 respectively, with exact named-head native-array correspondence and exact five-model means.

D6isoform6 completed all five folds: best checkpoint epochs 12, 11, 17, 14 and 13 (zero-based). Recorded per-fold training/prediction callback durations were approximately 300, 203, 196, 199 and 182 seconds. Score run `a1cb2ca5` gives D6 ST-RAE 0.869874 and Spearman 0.415397, compared with 0.825436 and 0.472752 for the fixed previous global ensemble. It improves over the reproduced challenge-only single/auxiliary models, but does not beat the previous global model. This is a member-level result; the four-model equal prediction average is evaluated separately in the root reproduction report.

## Completed native training

All 25 planned fits finished: four recipes with five parent-adapted folds each, plus five source-scaffold D6single folds. Every fold has a best native checkpoint, saved model, membership table, epoch history, all-head validation predictions and all-head test predictions. All five recipe/control exports are complete. Nothing was submitted.

| Recipe | Best zero-based epochs by fold | Per-fold recorded seconds, rounded |
|---|---|---|
| union18 | 48, 44, 27, 93, 28 | 521, 408, 273, 428, 232 |
| challenge8 | 18, 17, 10, 17, 33 | 65, 78, 63, 67, 78 |
| D6isoform6 | 12, 11, 17, 14, 13 | 300, 203, 196, 199, 182 |
| D6single1 parent | 10, 10, 14, 13, 12 | 11, 10, 13, 15, 26 |
| D6single1 source | 29, 14, 13, 8, 22 | 32, 25, 25, 20, 27 |

Each selected checkpoint was restored after patience-based stopping; the selected epoch is not the final executed epoch. Durations include fit and prediction operations instrumented inside each fold and reflect concurrent execution, not a standardized speed benchmark.

Union18 score run `bd89c19e` reports:

| Endpoint | Native ST-RAE | Previous global ST-RAE | Native Spearman | Previous global Spearman |
|---|---:|---:|---:|---:|
| 1A2 | 0.736395 | 0.670794 | 0.576032 | 0.611303 |
| 2C9 | 0.572785 | 0.443453 | 0.680005 | 0.742306 |
| 2D6 | 0.920374 | 0.825436 | 0.431863 | 0.472752 |
| 3A4 | 0.483195 | 0.391582 | 0.787438 | 0.825636 |

The faithfully implemented native union member does not improve the fixed prior global ensemble on these parent-CV endpoints. This establishes an actual source-training comparison; it is not evidence that the earlier custom transfer/LGBM pipeline was the published model. The final fixed four-member direct average, fold comparisons, distribution analysis and test candidate are handled in the main reproduction evaluation report.

Final independent export audits also passed: D6isoform6 `e5654aa2` and union18 `08b26a21`. Across all members the audits verify canonical keys, labels, confidence bounds, fold membership, exact native named-head prediction mapping and exact five-fold test means. See [the independent export audit](2026-09-12-scp-export-audit.md).
