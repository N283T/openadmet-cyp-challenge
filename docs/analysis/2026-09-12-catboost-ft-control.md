# CatBoost on fine-tuned Jeremy embeddings plus nested log2fc

The fixed CatBoost control did not improve the retained three-model ensemble. Its standalone macro ST-RAE is **0.638134**, compared with **0.628743** for LightGBM on the same feature recipe. Adding CatBoost as an equal fourth component changes the retained three-model score from **0.605039** to **0.607041**. Keep the three-model mean as the stronger observed result from this comparison.

## Recorded runs

- Encoder fine-tuning source: `20260912T060538153709Z-8d434eee`.
- CatBoost readout: `20260912T062520343225Z-49647baf`.
- Standard matched audit: `20260912T063041498362Z-415261e9`.
- Fixed three/four-component audit: `20260912T071605395541Z-e57a40d1`.

The CatBoost fit and standard audit completed before the conversation resumed; they were inspected and reused without rerunning training. The final comparison is a new recorded audit. All earlier artifacts remain intact.

## Configuration and verification

Inputs are the five fold-specific 600-dimensional embeddings from the completed Jeremy-medium multitask FT, plus four safe nested log2fc features from shared run `20260912T052549384296Z-bac5fb3a`. All 20 endpoint/fold models use 604 features and the unchanged parent-retained split.

CatBoost parameters were fixed in advance: 500 iterations, depth6, learning rate0.03, L2 leaf regularization5, MAE objective, seed42+fold, three CPU threads, no early stopping, and no hyperparameter search. Every fitted model was saved. No GPU training occurred. The backend addition leaves LGBM and TabPFN recipes and embedding-only PCA behavior intact.

Both audits verified 6,384 unique compound/endpoint/fold keys and identical labels, confidence bounds, neighbor markers, and split IDs. The four-component blend uses fixed 0.25 weights. The comparison's original three components are original Caruana24 (`96fce32b`), frozen Jeremy-medium TabPFN (`2f0fbfbd`), and fine-tuned Jeremy-medium plus nested-log2fc LightGBM (`730c8bea`). No weights were selected on these OOF results.

The upstream checkpoint provenance and available audit evidence are recorded in [the Jeremy medium control report](2026-09-12-jeremy-medium-frozen-control.md); this experiment inherits that context.

## Results

| Model | All OOF ST-RAE | Neighbor ST-RAE |
|---|---:|---:|
| FT + log2fc LightGBM | 0.628743 | 0.747805 |
| FT + log2fc CatBoost | 0.638134 | 0.774141 |
| Existing three-component mean | **0.605039** | **0.696275** |
| Four-component mean including CatBoost | 0.607041 | 0.706221 |

| CYP | LGBM | CatBoost | Three-component mean | Four-component mean |
|---|---:|---:|---:|---:|
| 1A2 | 0.711147 | 0.725236 | 0.686811 | 0.690163 |
| 2C9 | 0.505252 | 0.512399 | 0.479764 | 0.480066 |
| 2D6 | 0.851538 | 0.862402 | 0.838174 | 0.839139 |
| 3A4 | 0.447035 | 0.452499 | 0.415407 | 0.418796 |

CatBoost is slightly worse on ST-RAE for all four endpoints, and the fixed fourth-component addition also worsens every endpoint. The four-way mean slightly improves ordinary MAE (0.496912 to 0.496456), but that does not offset its worse challenge metric. Preserve CatBoost as a tested model-family control; do not replace the retained blend based on this result.
