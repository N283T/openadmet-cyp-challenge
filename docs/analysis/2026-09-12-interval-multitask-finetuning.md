# Official interval-error multitask fine-tuning

The official scoring source `artifacts/metric-reference/custom_scoring_functions.py`
defines soft absolute error as `max(lower - prediction, 0) + max(prediction - upper, 0)`.
Its relative-error denominator is the same error for the evaluation subset's mean
observed target. Predictions inside the confidence interval incur zero error.

The optional `training_loss=pic50_interval` uses that interval error for the four
pIC50 heads. Each head is normalized by its **inner-training** mean-predictor
interval error, so neither outer targets nor inner-validation targets determine
the training denominator. Bounds are joined by compound ID and endpoint, checked
for finite ordered values wherever a target is present, and saved for the inner
training compounds. All other heads use explicitly configured `auxiliary_loss`
(`mse` by default, or `mae`). Existing point MSE/MAE behavior remains unchanged.

The source scorer is declared as a run input. For every fold/head, an independent
constant-prediction calculation is asserted to agree with the official scorer
within 1e-12. Missing labels are masked. Every training step checks finite loss
and gradients. Full outer-compound exclusion applies to all eight heads;
checkpoint selection uses group-separated inner validation only.

## Smoke

Regression run `20260912T072034139020Z-92ad0df3` completed all five folds with two
epochs per fold, produced all 6,384 keyed OOF rows, and exported 6,896 embeddings
per fold. All confidence-bound, group-disjointness, outer-exclusion, and finite
gradient assertions passed. This short run is an execution check, not a model
quality comparison.

## Full trial

Regression run `20260912T072355846776Z-2e678154` starts from Jeremy's upstream
600-dimensional, depth-4 encoder bridge `20260912T060426566774Z-7604ffeb`, with
fresh prediction heads. The fixed configuration uses encoder LR 1e-5, head LR
1e-4, pIC50 weight 2, log2fc weight 1, maximum 160 epochs, patience 20, and
`checkpoint_selection=pic50_st_rae`. The launcher schedules direct-prediction
evaluation and a fold-specific LightGBM readout using embeddings plus the safe
nested predicted log2fc package `20260912T052549384296Z-bac5fb3a`.

The upstream encoder's public provenance and incomplete molecule manifest retain
the context documented in `2026-09-12-jeremy-medium-frozen-control.md`.

## Results

All scheduled work completed successfully. Both audits matched all 6,384 OOF
keys without weight search.

| Prediction | OOF ST-RAE | Neighbor ST-RAE | Audit run |
| --- | ---: | ---: | --- |
| Direct multitask pIC50 | 0.629854683 | 0.718056211 | `20260912T072937188761Z-26684ea2` |
| FT embeddings + nested log2fc, LightGBM | 0.627202554 | 0.732228620 | `20260912T073000950973Z-ce100a46` |

The LightGBM readout is run `20260912T072938900055Z-5d9f4a03`; its keyed output is
`outputs/external-oof.parquet`. It slightly improves on the earlier corresponding
FT-log2fc LightGBM component (0.628743095), but does not by itself supersede the
separately established five-component ensemble (0.5976166864). No further
training or post-hoc mixture search was performed for this trial.
