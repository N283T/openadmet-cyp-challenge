# Actual TDI Chemprop fine-tuning

Run: `classification/20260912T120135164747Z-1da37202` (completed successfully).

Two fixed recipes initialize their message-passing encoder from the external-data-only classification checkpoint `a3ea2ec3`. Both use a600-dimensional, depth5 mean-aggregation Chemprop encoder followed by a256-unit ReLU/dropout0.1 prediction head. Encoder/head AdamW learning rates are1e-5/1e-4, weight decay1e-4, batch size64, seed42+outer fold. The imported architecture implementation is the existing `finetune_external_multitask.py`; this run trains both the encoder and new classification head.

- `tdi_only`: two official boolean TDI heads, CYP2D6 and CYP3A4. Loss is mean masked per-head BCE.
- `tdi_aux`: the same BCE plus0.1 times the mean masked per-head standardized MSE of ten auxiliaries: four primary log2fc, four direct pIC50, and two TDI-condition pIC50. The0.1 coefficient weights the total auxiliary mean, not each head independently relative to total classification loss. Head means and standard deviations are computed only from each phase's training partition.

Classification labels are the supplied canonical confident calls:1,497 D6 and3,584 3A4, totaling5,081. Ambiguous missing calls are never converted to negatives. Auxiliary observation table `271e3020` supplies actual log2fc and TDI-condition potency; direct potency and official classification labels come from canonical measurements `a7063ccd`. These auxiliary values are targets, not inference features.

The outer split is the original UMAP assignment `edba6fc0`, with no regression-parent fold substitution. Inner fold is `(outer+1)%5`. Each inner fit excludes both outer and inner compounds across **all twelve heads**, selects the best inner two-head BCE within80 epochs with patience10, and obtains separate endpoint MCC thresholds on a fixed0.05–0.95 grid with0.01 increments. The encoder, new head, optimizer and RNG then reset to the same initialization; all outer-training compounds are used to refit for the selected epoch count. Outer labels never select epochs, scaling or thresholds.

Artifacts include per-phase membership tables covering all6,896 compounds, training-only scaler/count JSON, complete epoch histories, selected inner and final refit model states, all observed inner probabilities with outer-fold keys, outer OOF probabilities and calls, and750 test probabilities per endpoint per fold. Test results average the five probabilities and five inner thresholds; no test labels are accessed. No submission is made.

## Results and independent verification

| Recipe | D6 MCC | D6 ROC-AUC | 3A4 MCC | 3A4 ROC-AUC |
|---|---:|---:|---:|---:|
| TDI BCE only | 0.049559 | 0.520171 | 0.111170 | 0.634316 |
| TDI BCE + auxiliary regression | 0.025323 | 0.529741 | 0.121850 | 0.648071 |

The auxiliary heads modestly improve3A4. D6 remains weak: its slightly higher AUC does not translate into better MCC with the independently chosen thresholds. These outcomes do not establish this model family as an improvement over frozen-feature alternatives; the common classification comparison evaluates that separately.

Selected epochs across folds0–4 are19/2/12/36/61 for BCE-only and25/2/13/33/45 for auxiliary training. Every refit starts fresh from the external checkpoint and uses exactly its selected epoch count. No fit selected the80-epoch cap.

Independent recorded audit `classification/20260912T120615145022Z-7e38be86` checks all20 phase memberships and training-only scalers/counts, selected inner epochs against refit counts, all20 thresholds reconstructed from saved inner probabilities, exactly5,081 canonical OOF and1,500 averaged test rows per recipe, and exact test means. The audit passed. See its separate report for details. No submission was made.

This is a classification experiment distinct from the preceding regression community reproduction.
