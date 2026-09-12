# Full-data medium encoder for the five-component submission

The selected Jeremy-medium challenge fine-tuning recipe was refitted on every
available non-test pIC50 and actual primary-screen log2fc label. This is a
submission fit, not another OOF estimate.

- Fit run: `20260912T074319275987Z-483dcce3`.
- Independent verification: `20260912T074601765511Z-538aa567` (passed).
- Adopted OOF recipe: `20260912T060538153709Z-8d434eee`.
- Starting encoder-only bridge: `20260912T060426566774Z-7604ffeb`.
- Output: `outputs/full-model.pt`, `external-embeddings.parquet`,
  `direct-predictions.parquet`, `training-labels.parquet`, `membership.parquet`,
  `history.jsonl`, `effective-config.json`, and `summary.json`.

The five selected one-based OOF epochs were 94, 112, 56, 55, and 103. The full fit
used their median, **94 epochs**, with no reserved validation subset and no early
stopping on full-training error. All other adopted settings were retained:
600-dimensional bond-message-passing encoder, depth 4, norm aggregation / 100,
fresh 600→256→8 predictor with ReLU and 0.1 dropout; AdamW encoder learning rate
1e-5, predictor learning rate 1e-4, weight decay 1e-4, batch size 64, seed 42,
masked standardized MSE, pIC50 head weight 2, log2fc head weight 1, gradient norm
clip 5. Per-head means and standard deviations were refitted using all available
training labels. The model was initialized from the upstream encoder bridge,
not from one of the already challenge-trained fold models.

All 6,525 pIC50 labels were used: CYP1A2 1,412; CYP2C9 1,285; CYP2D6 1,493;
CYP3A4 2,335. All 17,504 actual log2fc labels were used: 4,376 per CYP. The union
contains 4,906 labeled compounds out of 6,146 non-test compounds. All 750 test
compounds were inference-only. No TDI targets or additional Emax heads were used.
Embeddings and direct predictions were exported for all 6,896 compounds with
compound IDs, `split_id=full_training`, and `fold=-1`.

The separate audit reconstructed the original label table exactly, checked every
exported ID, confirmed all training-only target scalers exactly, loaded the saved
model strictly, and recomputed 16 evenly spaced embedding/prediction rows on CPU.
Maximum absolute CPU/export differences were 5.96e-8 for embeddings and 7.72e-7
for predictions, within the declared numerical tolerance. All values were finite.
The adopted downstream tuned LightGBM component uses the exported embeddings,
not this predictor's direct pIC50 predictions.
