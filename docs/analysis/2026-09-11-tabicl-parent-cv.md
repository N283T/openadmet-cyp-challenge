# TabICL baseline on corrected parent-retained CV

The previous TabICL comparison used UMAP folds. This experiment runs it on the
corrected parent-retained map, using the already rebuilt nested auxiliary models.

## Fixed comparison

Two TabICL feature recipes, each fit for four endpoints and five outer folds:
- Compact: frozen CheMeleon PCA256 plus four nested predicted log2fc features.
- Full: compact plus Mordred.

TabICL2.2.0, checkpoint tabicl-regressor-v2-20260212.ckpt (same snapshotted weights
as the earlier UMAP test), eight estimators, ensemble batch1, CUDA, automatic
AMP/FA3/offloading, seed42+fold, mean readout. No calibration or quantile selection.
The compact construction follows the previously inspected 256+4 baseline recipe;
local nested cross-fitting and CV are the project's explicit evaluation protocol.

Use the exact saved parent-retained PCA, imputer and descriptor column list from
the LightGBM run. Reconstruct input matrices and reload its booster before each
endpoint's TabICL fits. Assert predictions match saved LightGBM OOF within1e-10.
Verify all four auxiliary parent membership exclusions again. Anchors always train;
validation remains identical across estimators and feature recipes.

This is40 TabICL fits. No auxiliary refitting is required because these auxiliary
models were already trained for this partition. No external labels or test scores.
Compare all-nonanchor and neighbor subsets separately using the official scorer,
reporting ST-RAE, MAE and Spearman. LightGBM comparator is the unblended full-feature
parent-retained model, not its subsequently selected local blend.

## Results

| Model | All-nonanchor ST-RAE | All-nonanchor MAE | Neighbor ST-RAE | Neighbor MAE |
|---|---:|---:|---:|---:|
| LightGBM full | .687660 | .537846 | .836776 | .547848 |
| TabICL compact256+4 | .724348 | .549602 | .883880 | .572975 |
| TabICL full | .702055 | .543191 | .833343 | .547310 |

| Endpoint | n neighbors | LightGBM full | TabICL compact | TabICL full |
|---|---:|---:|---:|---:|
| 1A2 |78| .866457 | .837975 | .844744 |
| 2C9 |46| .901990 | .891155 | .915235 |
| 2D6 |68| .982010 | 1.143111 | 1.000136 |
| 3A4 |107| .596648 | .663281 | .573256 |

The full TabICL recipe is numerically close to LightGBM on neighbors, but worse on
the full nonanchor population. Compact TabICL does not improve aggregate scores.
Neighbor Spearman is .510554 LightGBM, .484350 compact TabICL, .494127 full
TabICL. No broad replacement is supported by this fixed baseline comparison.
No claim is made about a public leaderboard score; these are local corrected-CV
results on the stated populations.

## Reproduction check and retained failed attempt

Fold4 first attempt `20260911T111644447819Z-ad17198d` failed the pre-fit LightGBM
reconstruction check on one CYP1A2 validation molecule: maximum prediction delta
.011962676. Recreating the original PCA fit/transform sequence using the same
training rows and seed yielded PCA parameters matching the saved state and a
maximum transformed-feature difference of2.384185791015625e-7 from the reloaded
transform. With that reconstruction the LightGBM predictions passed the original
1e-10 check. The tiny feature difference changed a tree branch. The failed run is
retained, and the repaired fold runs from a new config and source snapshot. No
error threshold was relaxed. Exact parameter and preprocessor files are inputs.

Successful folds0–4:
- `20260911T111242726192Z-f4e821b7`
- `20260911T111346031387Z-2da109d7`
- `20260911T111445421875Z-125ca4de`
- `20260911T111544943864Z-0c11f380`
- `20260911T111844883202Z-14521360`

Audit `20260911T111952347762Z-e9f6c457` checks matching keys/labels and registers
432 endpoint/fold metrics. Outputs oof.parquet, metrics.csv, macro.csv and
summary.json. Forty TabICL fits completed. Checks and44 DB-enabled tests passed
(16 existing dependency deprecation warnings). No new test predictions/submission.
