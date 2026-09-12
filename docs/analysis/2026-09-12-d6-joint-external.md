# Joint CYP2D6 graph training with public efficacy labels

## Concrete distinction from frozen efficacy features

The preceding `50533a8f` experiment froze external tree predictions and added those columns to a challenge readout. This experiment instead backpropagates public assay losses through the same600-dimensional molecular message-passing encoder that predicts challenge CYP2D6 pIC50. It tests the joint-learning distinction highlighted by the community report; it is an adaptation, not a full reproduction of its26-head union model or leaderboard placement.

The encoder starts from the existing external600 bridge `7604ffeb`. The four outputs are challenge2D6 direct pIC50, AID1851 maximum response, AID1851 activity at57.14µM and the active/inactive call. Actual newly used external data is the16,244-compound efficacy table `616de5d8`, whose connectivity blocks are disjoint from every challenge molecule. Inconclusive binary outcomes remain missing. External numeric responses are clipped to[-150,50] and kept in separate standardized heads; there is no conversion to challenge pIC50.

## Matched controls and objective

Public per-observed-label weights are0,0.05 and0.30. The0-weight control has the identical encoder, four-output head, optimizer, direct target, seed and epoch budget. Inactive uses BCE logits; the three continuous heads use standardized squared error. Every scaler uses only permitted training labels, with the binary head left at0/1. Inner CYP2D6 ST-RAE selects epochs; outer labels never select a checkpoint.

Let D be the number of training direct labels plus w times the number of observed public labels across the three public heads. The objective is

`L = (sum_direct MSE + w * sum_public observed-head loss) / D`.

Each epoch processes every training direct molecule and a uniform subset of4,096 external compounds without replacement. Public losses are multiplied by16,244/4,096=3.9658203125. Splitting this epoch into S mixed batches and multiplying each batch sum by S/D makes its averaged gradient an unbiased estimator of the stated full-label objective at a fixed parameter state. Ordinary successive optimizer updates still make the optimization path stochastic. Missing head labels contribute zero, rather than being replaced by a zero target. There is no public oversampling masquerading as a0.30 aggregate loss fraction.

For example, inner fold0 has1,095 direct labels and16,244/13,024/13,169 public labels. At public weight0.05 the public label-weight mass is65.96%; at0.30 it is much larger. Counts, exact weighted masses and fractions are serialized for every fold/arm.

The network reuses the Chemprop600 message-passing architecture and256-hidden head. AdamW uses encoder learning rate1e-5, head learning rate1e-4 and weight decay1e-4, dropout0.1 and seed42+fold. Inner selection caps at80epochs with patience10. No new hyperparameter search is conducted.

## Actual encoder gradient and convergence evidence

The first batch of every arm records the encoder gradient norm from direct and public losses separately. For inner fold0:

| Public weight | Direct encoder gradient norm | Public encoder gradient norm |
|---|---:|---:|
| 0 | 0.019187 | 0 |
| 0.05 | 0.006531 | 0.004331 |
| 0.30 | 0.001520 | 0.006046 |

Thus public labels change the shared encoder, and increasing their weight materially changes its gradient composition. This is not just a frozen feature transformation.

All15 inner-selection fits stopped by patience after21–49epochs, below the80epoch cap. Best epoch budgets by outer fold0–4 are:

- No public:13,14,11,13,14.
- Public0.05:14,39,13,11,15.
- Public0.30:23,22,23,21,23.

## Inner-checkpoint variant, retained as a separate experiment

Run `20260912T101831048657Z-478527aa` restores the selected inner-training checkpoint. It does not refit the encoder on its inner-validation molecules. That original variant remains recorded and is not overwritten.

| CYP2D6 ST-RAE | No public | Public0.05 | Public0.30 |
|---|---:|---:|---:|
| Direct graph output | 0.889216 | 0.892561 | 0.880635 |
| Learned600 + nested log2fc → LGBM | 0.843085 | 0.833024 | 0.833210 |

The latent readout improves by about0.010 relative to the matched no-public control. Its separate LGBM trains on all permitted outer-training pIC50 labels. This is developmental evidence for representation improvement, not yet an ensemble or leaderboard improvement.

Independent audit `20260912T102255549289Z-18175a32` passed exact1,457-molecule OOF key/label/bound matching, inner group separation, outer exclusion, external connectivity exclusion, scaler recomputation, per-label weight-mass calculations and encoder gradient checks.

## Full outer refits

Run `20260912T102254368746Z-570cbe07` resets each encoder/head/optimizer to the same initialization and trains all allowed outer-training challenge labels for exactly the selected epoch count above. It uses a snapshot of the inner run's epoch table. Every target scaler is recomputed on the full outer-training membership; all public rows remain challenge-disjoint. There is no outer validation checkpoint selection. These six OOF methods have a `_refit` suffix and are the final comparison for promotion.

Full outer refits completed. Independent audit `20260912T102527546486Z-a81e812b` passed all checks, including exactly the selected epoch budgets and inclusion of every allowed outer-training pIC50 label.

## Retained artifacts

`joint-oof.parquet` contains direct and LGBM readout predictions for each public weight. Each fold/arm directory preserves encoder checkpoint, all6896-compound600-dimensional embeddings, scaler and label-mass metadata, training history, and LGBM weights. `memberships.parquet`, `external-membership.parquet`, `encoder-gradient-audit.csv` and `progress.json` support independent review. Scripts are under `d6_joint_external*`; no packages or infrastructure were changed and no submission was made.

| Final full-outer CYP2D6 ST-RAE | No public | Public0.05 | Public0.30 |
|---|---:|---:|---:|
| Direct graph output | 0.893099 | 0.900522 | 0.884174 |
| Learned600 + nested log2fc → LGBM | 0.838061 | 0.841291 | 0.839558 |

Public0.30 improves the direct graph output by0.008925 versus its matched control. However, the latent-readout gain seen in the inner-checkpoint variant does not survive the full outer refit: the no-public LGBM control is now best. Thus the evidence supports a modest direct-output benefit under this high public weight, not a clear improvement to the current strong regression ensemble. The public label-weight mass in the full-outer0.30 arm is about91%, rather than30%.

Final comparison uses the `_refit` methods in `570cbe07/outputs/joint-oof.parquet`. The earlier unsuffixed methods remain in478527aa solely as their separately recorded variants. No result is replaced or silently selected between those two training protocols. All five new joint scripts pass Ruff; both independent audits completed successfully.
