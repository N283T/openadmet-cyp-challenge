# Neighbor ablations and the primary-response bottleneck

The user authorized local-parent trials but questioned whether they address the
whole problem. Test local labels under corrected parent-retained CV and diagnose
which information the model is failing to recover.

## Local predictions

Use only endpoint-labeled outer-training molecules (including retained anchors).
Morgan radius2/2048 Tanimoto, k=1/5/20, similarity-to-the-fourth weighted means,
fixed blend weights .25/.5/1. Compare all-training pools and anchor-only pools;
for weight .5 also gate by clip((maximum similarity−.3)/.3,0,1). No validation
labels enter local predictions. Evaluate all nonanchors and the same neighbor set.

| Method | Neighbor ST-RAE | Neighbor MAE | All-nonanchor ST-RAE |
|---|---:|---:|---:|
| Current full-feature LightGBM | .836776 | .547848 | .687660 |
| All-training k20, weight .25 | .785259 | .535877 | .689139 |
| All-training k20, gated weight .5 | .789704 | .538116 | .687484 |
| Anchor-only k20, gated weight .5 | .826443 | .561957 | .688401 |

Best fixed local recipe within this exploratory screen improves neighbor scores,
but not the full population. The similarity-gated recipe is more balanced by
ST-RAE. No recipe has an independent selection-adjusted evaluation. Parent-only
information is not the strongest local signal; broader labeled neighbors help more.
K20/.25 endpoint neighbor ST-RAE is .813376/.818694/.916450/.592516 in
1A2/2C9/2D6/3A4 order. The single-nearest-label approach is not generally competitive.

## Observed-primary diagnostic

Fit a decreasing isotonic primary-log2fc→direct-pIC50 mapping on each outer
training set, then evaluate using **measured** validation primary response. This
is diagnostic, never a deployable candidate. Compare baseline on exactly the same
primary-available rows (3A4 has only56 of107 neighbors with public primary data).

| Endpoint | Matched n | Current ST-RAE | Observed-primary diagnostic ST-RAE |
|---|---:|---:|---:|
| 1A2 |78| .866457 | .337052 |
| 2C9 |46| .901990 | .456398 |
| 2D6 |68| .982010 | .648940 |
| 3A4 |56| .524975 | .260761 |

Macro .818858→.425788 on this matched population. This is evidence of available
primary-response signal, not an achievable hidden-test score, a formal upper bound,
or proof that retraining the auxiliary learner can recover the whole difference.

## Check actual predicted-primary quality

Reconstruct the four-model outer-held-group average from the ten newly trained
auxiliary models and join by compound/outer-fold/enzyme to the primary measurements.
Every prediction has exactly four contributing excluded-pair models.

| Endpoint | Neighbor primary Spearman | log2fc RMSE | Primary/pIC50 residual Spearman |
|---|---:|---:|---:|
| 1A2 | .4206 | .7083 | −.7286 |
| 2C9 | .5917 | .3443 | −.7021 |
| 2D6 | .2790 | 1.3527 | −.5149 |
| 3A4 | .7472 | 1.0879 | −.6801 |

The negative association is consistent with inhibition's sign: predicting log2fc
insufficiently negative tends to accompany underpredicted pIC50. It is descriptive,
not causal. Primary-response prediction is a concrete bottleneck candidate,
especially for2D6 and1A2; nearest-parent placement alone is insufficient.

Structural error strata are saved separately. High similarity is not uniformly
low error, and high-similarity counts are sparse (only3 neighbors above.6 for1A2).
Do not equate nearest-label disagreement with proven experimental activity cliffs.

## Decisions and lineage

Retain the local/gated blends as exploratory alternatives, not a blanket replacement
for the full-population baseline. Prioritize comparison of molecular representations
and auxiliary learners, plus isoform-specific readouts, over further global shifts.

- Local and observed-primary audit: `20260911T110644999918Z-f0e04121`.
  Outputs predictions.parquet, metrics.csv, macro.csv, diagnostics.parquet,
  error-strata.csv, summary.json.
- Predicted-primary residual audit: `20260911T110751082997Z-dc7a1851`.
  Outputs aligned-residuals.parquet, metrics.csv, summary.json.

Both runs are registered with complete config/input snapshots. No new remote
submission. Source checks pass.
