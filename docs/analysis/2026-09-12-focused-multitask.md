# Focused interval multitask experiment

The submitted ensemble improved CYP2C9 and CYP3A4 but retained weak CYP1A2 and
CYP2D6 performance. This bounded comparison changes the training objective and
checkpoint selection, retaining all four potency and four measured primary
screen log2FC heads. Other targets remain auxiliary supervision.

## Prespecified comparison

Both new arms initialize the same upstream Jeremy medium encoder-only bridge
`20260912T060426566774Z-7604ffeb`, use the retained parent CV, and keep all outer
nonanchor compounds out of every training and checkpoint-selection head.

- pIC50 weights: CYP1A2 2, CYP2C9 0.5, CYP2D6 2, CYP3A4 0.5.
- Interval-distance pIC50 loss normalized by inner-training mean-predictor
  interval error; log2FC uses standardized MSE.
- Checkpoint selection: arithmetic mean inner ST-RAE for CYP1A2 and CYP2D6.
- Arm A: each log2FC head weight 0.25.
- Arm B: each log2FC head weight 1.0.
- Shared settings: five outer folds, inner molecular-group holdout 10%, anchors
  retained in training, seed42, maximum160 epochs, patience20, batch64,
  encoder LR1e-5, head LR1e-4, weight decay1e-4, dropout0.1.
- Produce all 6,384 matched OOF labels and fold embeddings; evaluate direct
  predictions and fixed LGBM on embedding600 plus nested log2FC4.

`focus_multitask.py` is a separate copy of the existing interval trainer with
only endpoint loss weights and selection-target filtering changed. The original
trainer remains untouched. The earlier equal-target interval run
`20260912T072355846776Z-2e678154` is the multitask control.

No test labels, LB labels, or tuned OOF blend weights enter this comparison.

## Arm A completed

- Training: `20260912T082146155304Z-f5b1cef2`.
- Direct audit: `20260912T082505957126Z-5b9497e4`.
- LGBM+log2FC readout: `20260912T082508554737Z-6457cb39`.
- Readout audit: `20260912T082540622249Z-42d23a38`.
- Independent membership/history verification:
  `20260912T082601915661Z-b7105514` (passed).

| Model | CYP1A2 | CYP2C9 | CYP2D6 | CYP3A4 | Macro |
|---|---:|---:|---:|---:|---:|
| Equal-pIC-weight interval FT + LGBM/log2FC control | 0.711416 | 0.502763 | 0.846144 | 0.448488 | 0.627203 |
| Focused aux0.25, direct | 0.740913 | 0.536862 | 0.866215 | 0.461705 | 0.651424 |
| Focused aux0.25, LGBM/log2FC | 0.711941 | 0.512063 | 0.856750 | 0.447974 | 0.632182 |

This arm does not improve the two focus endpoints relative to the comparable
existing LGBM readout. The matched verification checks all 6,384 OOF rows,
outer compound exclusions, inner group disjointness, and that every epoch's
selection score equals the mean of the two designated inner endpoint scores.

## Arm B completed

- Training: `20260912T082603039086Z-4dcb4559`.
- Direct audit: `20260912T082845486897Z-0342f49c`.
- LGBM+log2FC readout: `20260912T082847782759Z-b1e19c7b`.
- Readout audit: `20260912T082916732176Z-e435c1e4`.

| Model | CYP1A2 | CYP2C9 | CYP2D6 | CYP3A4 | Macro |
|---|---:|---:|---:|---:|---:|
| Focused aux1, direct | 0.735000 | 0.524647 | 0.858658 | 0.449637 | 0.641985 |
| Focused aux1, LGBM/log2FC | 0.724084 | 0.512435 | 0.857508 | 0.449007 | 0.635758 |

Increasing the auxiliary weight improves direct predictions relative to arm A,
but the embedding readout becomes weaker. Neither focused objective improves
CYP1A2 or CYP2D6 relative to the existing equal-weight interval embedding
readout. These two prespecified tests do not support adopting target loss
reweighting/checkpoint focusing as the next model improvement. Both retain the
full four-endpoint OOF key set; no submission or ensemble replacement was made.

Independent arm B verification `20260912T082944602413Z-9f7909a5` also passed
all membership, finite-loss, checkpoint-selection and 6,384-row OOF checks.
