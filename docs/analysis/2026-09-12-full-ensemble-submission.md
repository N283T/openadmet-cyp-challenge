# Full-data five-component regression submission

The user authorized submission of the uncalibrated five-component recipe whose
development OOF macro official ST-RAE is **0.5964141092282933**. Every challenge
readout was refitted on all available endpoint labels before test inference.
The official API returned **sent** for attempt `20260912T074820067880Z-3d2e5c50`.
Leaderboard scoring is pending; the OOF score is not a leaderboard result.

## Exact artifact and receipt

- Package run: `20260912T074740645032Z-45d23980`.
- Candidate: `20260912T074758708521Z-da569eed`.
- Prepared submission: `20260912T074802231782Z-eaafa860`.
- Send attempt: `20260912T074820067880Z-3d2e5c50`.
- File SHA256: `4e70a727c532cad2e19aff04dba4f888c874b7e654bda9bbafbf084ce10287e6`.
- Uploaded file: `artifacts/submissions/regression/20260912T074802231782Z-eaafa860/predictions.parquet`.
- Receipt: `artifacts/submission-attempts/regression/20260912T074820067880Z-3d2e5c50/result.json`.

The package retains Parquet/CSV copies, long keyed predictions, exact component
hashes and weights, full training membership checks, test distribution summaries
and component correlations. Dry-run validation passed before the single upload.
The existing private profile was reused. No public repository/report release or
classification submission was made.

## Full-data recipe

All components have weight **0.2**; no offset, slope, clipping or fitted ensemble
weight was applied. Original pIC50 labels remain in their original units.

| Component | Full-data run | Deployment details |
|---|---|---|
| Caruana LGBM | `20260912T074545177210Z-0d561e17` | Average five saved inner-selected weights per endpoint; refit 77 positively weighted base models across four endpoints, fixed seed42 |
| Frozen medium600 → TabPFN | `20260912T074409483212Z-f3024998` | Original upstream frozen features; TabPFN2.6, eight estimators, all endpoint labels |
| Upstream scores47 + log2fc4 → TabPFN | `20260912T074557646349Z-da857069` | Original 47 public predictor outputs, eight-estimator TabPFN2.6 |
| Adapted medium600 + log2fc4 → tuned LGBM | `20260912T074624192288Z-bd58b762` | Refit each of five CV-selected parameter sets on all labels; equal mean within this component |
| Adapted-public CheMeleon2048 → PCA256 + log2fc4 → TabICL | `20260912T074558733332Z-98ab8ceb` | Endpoint training-only embedding scaler/PCA; append log2fc after PCA; eight estimators |

The encoder fine-tune `20260912T074319275987Z-483dcce3` starts from the same
upstream medium encoder as the adopted OOF run. It uses every **6,525 pIC50 +
17,504 primary log2fc labels** across 4,906 compounds. Fixed **94 epochs** is the
median of the five saved one-based checkpoint epochs (94,112,56,55,103); no
validation subset is withheld in this final fit.

The auxiliary CheMeleon model `20260912T074339967571Z-ff194fcf` fits all **4,376
measured compounds** for 13 epochs, the median of earlier selected epoch counts.
All test auxiliary predictions come from this full model. Readout training rows
retain the recorded recipient-excluded cross-fitted auxiliary predictions,
averaging four pair-excluded models. Thus auxiliary train features remain
cross-fitted while the deployed auxiliary predictor uses all measured labels.

The public upstream checkpoint corpus details remain as documented in
[the OOF result](2026-09-12-oof-goal-result.md). This refit does not change their
provenance or retrain frozen external features.

## Validation

Every component membership exactly matches all original endpoint labels:
1A2 **1,412**, 2C9 **1,285**, 2D6 **1,493**, 3A4 **2,335**. Each predicts the same
750 test compounds and four endpoints. Test compounds are absent from supervised
challenge training. Package joins use compound IDs and endpoints; original
Molecule_Name/SMILES pairs match the reference exactly. All predictions are finite
and pass the official local submission contract.

- Encoder audit `20260912T074601765511Z-538aa567`: exact label/scaler coverage,
  checkpoint epochs, 6,896 exported identities, and independent CPU checkpoint
  inference on 16 molecules; maximum embedding difference 5.96e-8.
- Auxiliary audit `20260912T074551656789Z-beecf03c`: training cross-fit features
  exactly match their source; saved full model reproduces all test log2fc values
  with maximum absolute difference 0.
- Caruana audit `20260912T074646849858Z-1a37f7b6`: exact original labels and
  saved-inner-selection weight transfer, no test/training overlap.
- Package `45d23980`: independent alignment/membership checks for all five
  components, 3,000 predictions per component, equal mean and original identifiers.
- `pixi run check` passed; `pixi run test`: 37 passed, 7 skipped.

Final test prediction mean / standard deviation:

| CYP | Mean pIC50 | Standard deviation |
|---|---:|---:|
| 1A2 | 4.934455 | 0.690488 |
| 2C9 | 4.901427 | 0.683965 |
| 2D6 | 4.638156 | 0.405163 |
| 3A4 | 4.742966 | 0.902667 |

These are diagnostics, not calibration targets. No decision was changed based
on the test distribution. All-data deployment uses fixed CV-derived epochs,
parameters and weights; its future LB score must be read from the scoring receipt.
