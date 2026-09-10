# Mordred and direct fine-tuning screen

Recorded 2026-09-10 UTC. Compare three requested alternatives under the existing
UMAP five-fold assignment. No routine random/scaffold evaluation, feature search,
threshold tuning or blinded-test predictions.

## Prespecified configurations

1. Mordred alone → TabPFN-2.6, eight estimators, seed 42 + fold.
2. Saved fold-safe PCA256 + nested predicted log2fc4, concatenated with Mordred
   → the same TabPFN. Reuse the exact 260-dimensional matrices from RL-0007,
   preserving its PCA and auxiliary cross-fitting protocol.
3. Direct CheMeleon fine-tuning: a four-head pIC50 regressor and a separate two-head
   TDI classifier per outer fold. Start from the original foundation weights,
   with new BatchNorm and 256-unit prediction heads. Train the encoder too.

Mordred source: `20260909T132918829974Z-bc81af58`. Use all 1,613 available 2D
columns before train-only removal of entirely missing and constant columns.
Convert nonfinite values to missing; fit median imputation on training rows.
Persist retained column names, imputer states and fitted TabPFN states. There is
no supervised feature selection or external PCA on Mordred. TabPFN's native
internal preprocessing is enabled. Base model checkpoints are copied into each run.

Direct FT uses Adam at a fixed 1e-4 learning rate, batch size 64, maximum 100 epochs,
patience 15 and a seeded 90/10 inner partition for checkpoint selection. Outer
holdout compounds never participate in updates, scaling or early stopping.
Regression targets are standardized on inner training only. The objective equally
averages observed per-head MSE (regression) or BCE-with-logits (TDI); missing labels
contribute no gradient. There are no class weights. Classification predictions use
sigmoid and the unchanged 0.5 threshold. Selected weights are not refitted on the
inner validation rows. This is one initial FT recipe, not tuned FT performance.

The checkpoint architecture remains width 2048, depth 6 and mean pooling. Preserve
best/last weights, optimizer state, inner/outer membership, target scalers, learning
curves and keyed OOF outputs. The last weights are recovery artifacts rather than
exact resume checkpoints (no RNG snapshot). Foundation checksum provenance remains
as documented in RL-0004. Each run includes the installed environment inventory;
TabPFN effective constructor arguments are also saved per endpoint.

`scripts/launch_model_screen.py VARIANT` creates configurations under
`configs/screen/` and invokes recorded runs. Variants are `mordred`, `stack_mordred`
and `direct_ft`. Existing configurations are never overwritten. The inherited
`features` metadata field is unused; `variant` controls the actual feature/model
branch. `scripts/finalize_model_screen.py VARIANT` runs the official scoring audit
and verifies lineage/fold membership before inserting structured predictions.
All attempts, including failures, remain recorded. Audit registration is atomic.

Interpret results against the same pooled OOF metrics from RL-0006 and RL-0007.
Mordred-only versus Mordred-augmented estimates the effect of the combined 260
features; it does not isolate log2fc from CheMeleon PCA. Direct FT changes the
learner and representation jointly. Reusing the same CV splits across exploratory
choices creates selection uncertainty; no significance claim is warranted.

## Completed Mordred-only evaluation

Audit `20260910T123026818816Z-7968c1c0` passed: 11,606 unique OOF pairs and
agreement with the official regression scorer. Macro ST-RAE 0.753071,
MAE 0.573650, MCC 0.123266 and AUROC 0.652575. This does not improve on the
previous compact 260-feature model's ST-RAE 0.729505 or MCC 0.124720.
Database load: `20260910T123028057829Z-bdefcd30`.

Mordred-only fold runs, in fold order:

- `20260910T122243830616Z-ba02a70e`
- `20260910T122343854890Z-3a03f908`
- `20260910T122509144654Z-bd47d2b2`
- `20260910T122650970011Z-303d797c`
- `20260910T122823042449Z-ddd5cc8e`

## Completed stacked Mordred evaluation

Audit `20260910T123856868493Z-243a1e23` passed with the same 11,606 OOF pairs.
Macro ST-RAE 0.724704, MAE 0.559647, MCC 0.129339 and AUROC 0.668071. Adding
Mordred to the saved 260 features improves ST-RAE by about 0.0048 and MCC by about
0.0046 versus RL-0007, a small gain on this fixed split. CYP2D6 TDI remains
all-negative at 0.5 (AUROC 0.542860); CYP3A4 TDI MCC is 0.258678 and AUROC 0.793282.
Database load: `20260910T123857997147Z-b7077239`.

Stacked Mordred fold runs, in fold order:

- `20260910T123025501770Z-492cf5f6`
- `20260910T123205649940Z-20f4f93f`
- `20260910T123336884511Z-63b24af3`
- `20260910T123510230339Z-c92321ca`
- `20260910T123640049734Z-b387c093`

## Completed direct FT evaluation and decision

Audit `20260910T124213825407Z-df00102c` passed for all 11,606 OOF pairs.
Database load: `20260910T124214870821Z-1fcb7dab`. Macro ST-RAE 0.824623,
MAE 0.610799, MCC 0.115048 and AUROC 0.648803. Regression is worse than the
previous readouts for all four endpoints. CYP2D6 TDI has some positives but still
weak pooled MCC (0.040064); CYP3A4 TDI MCC is 0.190033. Reject this FT recipe as
a replacement for the current baselines, while preserving all artifacts. This
single fixed learning-rate trial does not establish that direct FT cannot work.

Direct FT fold runs, in fold order:

- `20260910T122404575916Z-3524a751`
- `20260910T122802955829Z-867dbeeb`
- `20260910T123155955776Z-91d5fc70`
- `20260910T123546583878Z-34e31c6d`
- `20260910T123920632611Z-29f2ba24`

| Pipeline | Macro ST-RAE ↓ | MAE ↓ | MCC ↑ | AUROC ↑ |
|---|---:|---:|---:|---:|
| Earlier LightGBM + Mordred + Morgan | 0.753365 | 0.580243 | 0.165203 | 0.647892 |
| Earlier fine-tuned embedding + TabPFN | 0.731115 | 0.562771 | 0.112674 | 0.661075 |
| Earlier PCA256 + log2fc4 + TabPFN | 0.729505 | 0.563181 | 0.124720 | 0.645634 |
| Mordred + TabPFN | 0.753071 | 0.573650 | 0.123266 | 0.652575 |
| PCA256 + log2fc4 + Mordred + TabPFN | 0.724704 | 0.559647 | 0.129339 | 0.668071 |
| Direct CheMeleon FT | 0.824623 | 0.610799 | 0.115048 | 0.648803 |

The augmented TabPFN is the lowest observed regression ST-RAE in this screen;
the earlier LightGBM still has the highest classification MCC in this comparison.
Keep regression augmentation as a small observed improvement, not a statistically
established gain. No automated selection of per-endpoint winners or ensemble was
performed. All 15 fold runs completed with no failed experiments. Nonfinite-column
NumPy warnings during Mordred filtering are retained in logs; training-empty
columns were removed before imputation. Checks and 32 tests passed, including
masked-loss gradient checks for missing labels and empty-head handling.
