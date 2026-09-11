# Complete downstream replacement after auxiliary adoption

The user requested an actual model replacement. The earlier LGBM-only update
did not complete the downstream TabICL/TabPFN work. The earlier trial-stage
preference for LightGBM does not restrict this replacement.

This campaign fits both TabICL 2.2.0 and TabPFN 2.6 on all five parent-retained
folds, four direct-inhibition pIC50 endpoints each: 20 fitted models per backend.
Each uses the adopted Rustuna auxiliary models and the full feature recipe:
frozen CheMeleon PCA256, four nested predicted log2fc values and imputed Mordred.
This campaign updates the full models; it does not rerun the compact ablation.

The shared readout reconstruction validates every fold against the saved updated
LightGBM predictions with strict tolerance before fitting. PCA and descriptor
imputation remain training-only. Recipient-fold auxiliary features use the
appropriate excluded pair; validation and test average the four auxiliary models
containing the outer fold. All anchor compounds stay in readout training.

Both backends use eight estimators, seed 42+fold, and mean predictions. TabPFN
explicitly uses the 2.6 factory and local 2.6 checkpoint, not the installed package's
latest model default. TabICL saves fitted joblib models; TabPFN saves fitted states.
Their complete preprocessing and fixed checkpoint inputs remain in immutable runs.

For each backend, outputs include all 6,384 nonanchor labeled OOF rows, 15,000
fold-level test predictions and 3,000 five-fold-averaged endpoint predictions,
also exported as a 750-molecule regression table with names and SMILES.
The final audit checks keys, labels, confidence intervals, folds, endpoints and
complete test coverage. Official ST-RAE is compared on matched rows against the
previous full TabICL and updated LightGBM. Fold 0 remains the only HPO-excluded
evaluation fold; pooled results are post-selection evaluation of the adopted recipe.

## Completed runs and results

| Fold | TabICL | TabPFN 2.6 |
|---|---|---|
| 0 | `20260911T130541636074Z-98075594` | `20260911T131313006781Z-a8508a11` |
| 1 | `20260911T130703156203Z-8d3c8859` | `20260911T131412685950Z-800631af` |
| 2 | `20260911T130823926945Z-a08d12eb` | `20260911T131512444869Z-ee4cdffc` |
| 3 | `20260911T130944812832Z-48eff1e0` | `20260911T131612138159Z-011ed838` |
| 4 | `20260911T131106047746Z-2ea2c432` | `20260911T131712066480Z-7bc3ef15` |

Final audit `20260911T131826604125Z-ebb033cc` registers 576 metrics and verifies
6,384 OOF rows and 750 complete test compounds for each backend. All 40 fitted
model assets exist. No failed training runs occurred in this replacement campaign.

| Model | All ST-RAE | All MAE | Neighbor ST-RAE | Neighbor MAE |
|---|---:|---:|---:|---:|
| Previous TabICL full | 0.702055 | 0.543191 | 0.833343 | 0.547310 |
| Updated TabICL full | 0.695341 | 0.540527 | 0.834264 | 0.552283 |
| Updated TabPFN 2.6 full | 0.688696 | 0.536118 | 0.854329 | 0.556130 |
| Updated LightGBM full | 0.684192 | 0.537255 | 0.850116 | 0.557958 |

TabICL improves overall versus its matched previous model; neighbor ST-RAE is
essentially unchanged. TabPFN 2.6 is now actually fitted on the adopted parent-CV
features; no matched previous parent-CV TabPFN result existed. Its earlier UMAP
scores are not used as a before/after control here.

`configs/regression/current-foundation-models.json` now points to actual fitted
models and the audited OOF/test outputs for both backends. The two
`*-regression.parquet` outputs contain molecule names, SMILES and all four direct
pIC50 predictions. Check and 44 DB-enabled tests pass. No remote submission.
