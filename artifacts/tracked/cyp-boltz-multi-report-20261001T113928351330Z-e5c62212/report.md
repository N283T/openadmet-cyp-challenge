# Multiple trunk models in the incumbent ensemble

Analysis run: `regression/20261001T113601716233Z-bf06a273`. No base-model refitting, candidate or upload.

All four submitted regression cohorts, labels, confidence bounds and parent-series folds are reused exactly. Each target has 72 trunk recipes: 12 views, with or without molecular features, and three readouts. Pairwise correlations describe pooled OOF.

## Trunk diversity

| Target | Prediction correlation median (10th–90th percentile) | Residual correlation median |
|---|---:|---:|
| cyp1a2 | 0.887 (0.852–0.928) | 0.947 |
| cyp2c9 | 0.898 (0.860–0.945) | 0.926 |
| cyp2d6 | 0.836 (0.762–0.898) | 0.953 |
| cyp3a4 | 0.925 (0.895–0.967) | 0.901 |

## Held-fold results

ST_RAE is lower-is-better. Counts and blending strength use inner held-fold error within the other four folds. Greedy block fitting uses distinct trunk recipes, with block sizes 1, 2, 3, 5 and 10. A block is one prediction column, blended at weight 0, .05, .1, .2, .3, .5 or 1. Multiple additions give each selected recipe one equal-strength slot and preserve incumbent relative weights. Replacement greedily replaces up to three distinct incumbent atoms at their saved weights; its cap is inner-fold selected. D6 replaces at most one branch to retain its hierarchy. Full-bank selection uses Caruana over incumbent atoms plus all 72 trunks, with steps 5, 10, 20 or 50 chosen by inner folds; an incumbent-only control is included.

| Target | Current | Previous one-member placement | Previous one-trunk bank extension | Trunk block added | Multiple added | Multiple replaced | Full bank | Full-bank control |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| cyp1a2 | 0.670833 | 0.669032 | 0.669839 | 0.668365 | 0.667602 | 0.665644 | 0.669440 | 0.674055 |
| cyp2c9 | 0.454711 | 0.459567 | 0.458215 | 0.454248 | 0.457067 | 0.454711 | 0.457663 | 0.456038 |
| cyp2d6 | 0.819829 | 0.819451 | 0.815745 | 0.816442 | 0.816601 | 0.819829 | 0.815874 | 0.820524 |
| cyp3a4 | 0.394567 | 0.395694 | 0.396542 | 0.395140 | 0.394567 | 0.394567 | 0.396687 | 0.394567 |

## Paired uncertainty

Parent-series group bootstrap, 2,000 draws; conditional on saved OOF and fitted selections. Intervals do not cover repeating selection or choosing the best reported method. All improving methods have intervals spanning zero.

| Target | Method | Reference | Delta ST_RAE | Conditional 95% interval |
|---|---|---|---:|---:|
| cyp1a2 | trunk_block_added | current | -0.002468 | [-0.009620, 0.004657] |
| cyp1a2 | multiple_added | current | -0.003231 | [-0.009797, 0.003552] |
| cyp1a2 | multiple_replaced | current | -0.005189 | [-0.013492, 0.002894] |
| cyp1a2 | full_bank | current | -0.001393 | [-0.010908, 0.008437] |
| cyp1a2 | full_bank | original_bank_control | -0.004615 | [-0.013743, 0.004367] |
| cyp2c9 | trunk_block_added | current | -0.000463 | [-0.004203, 0.003083] |
| cyp2c9 | multiple_added | current | 0.002356 | [-0.002252, 0.007219] |
| cyp2c9 | multiple_replaced | current | 0.000000 | [0.000000, 0.000000] |
| cyp2c9 | full_bank | current | 0.002951 | [-0.003356, 0.009062] |
| cyp2c9 | full_bank | original_bank_control | 0.001625 | [-0.003692, 0.007008] |
| cyp3a4 | trunk_block_added | current | 0.000573 | [0.000068, 0.001108] |
| cyp3a4 | multiple_added | current | -0.000000 | [-0.000000, 0.000000] |
| cyp3a4 | multiple_replaced | current | 0.000000 | [0.000000, 0.000000] |
| cyp3a4 | full_bank | current | 0.002120 | [-0.000479, 0.004681] |
| cyp3a4 | full_bank | original_bank_control | 0.002120 | [-0.000479, 0.004681] |
| cyp2d6 | trunk_block_added | current | -0.003387 | [-0.008832, 0.002019] |
| cyp2d6 | multiple_added | current | -0.003228 | [-0.012118, 0.005733] |
| cyp2d6 | multiple_replaced | current | 0.000000 | [0.000000, 0.000000] |
| cyp2d6 | full_bank | current | -0.003955 | [-0.011477, 0.003380] |
| cyp2d6 | full_bank | original_bank_control | -0.004651 | [-0.012170, 0.002804] |

## Selected compositions

Recipe identifiers below are exact reproducibility keys. In the 1A2 replacements, the main pair is ligand-only TabPFN and pocket/HEM bidirectional TabPFN; some folds also add pocket/HEM bidirectional LGBM. No molecular features are used in these replacement recipes.

- cyp1a2, fold 0, block_added, parameter 0.2: `["ligand_only--trunk_molecular--tabpfn_top500", "core_bidirectional_hem--trunk--lgbm_l1"]`
- cyp1a2, fold 0, multiple_added, parameter 2.0: `["ligand_only--trunk_molecular--tabpfn_top500", "core_bidirectional_hem--trunk--lgbm_l1"]`
- cyp1a2, fold 0, multiple_replaced, parameter 2.0: `[["chemeleon_pca256_log2fc_tabicl", "ligand_only--trunk--tabpfn_top500"], ["tabpfn35_fast_anchor59_raw_scores_logfc", "core_bidirectional_hem--trunk--tabpfn_top500"]]`
- cyp1a2, fold 0, full_bank, parameter 0.26: `{"tabpfn35_fast_pub030_raw": 0.02, "lgbm_all_pool_top512": 0.04, "lgbm_jeremy_all_log2fc": 0.3, "lgbm_dense_pool_top256": 0.1, "tabpfn35fast_all_pool_top512": 0.28, "core_bidirectional_hem--trunk--lgbm_l1": 0.02, "core_bidirectional_hem--trunk--tabpfn_top500": 0.08, "ligand_only--trunk--tabpfn_top500": 0.12, "ligand_only--trunk_molecular--tabpfn_top500": 0.04}`
- cyp1a2, fold 0, original_bank_control, parameter 0.0: `{"jeremy_scores_log2fc_tabpfn26": 0.02, "jeremy_finetuned_log2fc_lgbm": 0.06, "chemeleon_pca256_log2fc_tabicl": 0.04, "tabpfn35_fast_pub030_raw": 0.12, "lgbm_all_pool_top512": 0.06, "lgbm_jeremy_all_log2fc": 0.3, "lgbm_dense_pool_top256": 0.14, "tabpfn35fast_all_pool_top512": 0.26}`
- cyp1a2, fold 1, block_added, parameter 0.2: `["core_bidirectional_hem--trunk_molecular--lgbm_l1", "ligand_only--trunk--tabpfn_top500", "core_bidirectional_hem--trunk--tabpfn_top500"]`
- cyp1a2, fold 1, multiple_added, parameter 3.0: `["core_bidirectional_hem--trunk_molecular--lgbm_l1", "ligand_only--trunk--tabpfn_top500", "core_bidirectional_hem--trunk--tabpfn_top500"]`
- cyp1a2, fold 1, multiple_replaced, parameter 2.0: `[["jeremy_finetuned_log2fc_lgbm", "ligand_only--trunk--tabpfn_top500"], ["jeremy_scores_log2fc_tabpfn26", "core_bidirectional_hem--trunk--tabpfn_top500"]]`
- cyp1a2, fold 1, full_bank, parameter 0.2: `{"tabpfn35_fast_pub030_raw": 0.2, "lgbm_all_pool_top512": 0.1, "lgbm_jeremy_all_log2fc": 0.1, "lgbm_dense_pool_top256": 0.1, "tabpfn35fast_all_pool_top512": 0.3, "ligand_only--trunk--tabpfn_top500": 0.2}`
- cyp1a2, fold 1, original_bank_control, parameter 0.0: `{"jeremy_frozen_tabpfn26": 0.05, "chemeleon_pca256_log2fc_tabicl": 0.05, "tabpfn35_fast_pub030_raw": 0.25, "lgbm_all_pool_top512": 0.1, "lgbm_jeremy_all_log2fc": 0.15, "lgbm_dense_pool_top256": 0.2, "tabpfn35fast_all_pool_top512": 0.2}`
- cyp1a2, fold 2, block_added, parameter 0.3: `["core_bidirectional_hem--trunk_molecular--lgbm_l1", "ligand_only--trunk_molecular--tabpfn_top500", "core_bidirectional_hem--trunk--tabpfn_top500", "ligand_only--trunk_molecular--lgbm_l1", "core_bidirectional_hem--trunk--lgbm_l1", "ligand_only--trunk--tabpfn_top500", "core_separate_fe--trunk_molecular--tabpfn_top500", "pxr_protein--trunk_molecular--lgbm_l1", "core_distribution_hem--trunk_molecular--lgbm_l1", "core_separate_fe--trunk--tabpfn_top500"]`
- cyp1a2, fold 2, multiple_added, parameter 5.0: `["core_bidirectional_hem--trunk_molecular--lgbm_l1", "ligand_only--trunk_molecular--tabpfn_top500", "core_bidirectional_hem--trunk--tabpfn_top500", "ligand_only--trunk_molecular--lgbm_l1", "core_bidirectional_hem--trunk--lgbm_l1"]`
- cyp1a2, fold 2, multiple_replaced, parameter 3.0: `[["tabpfn35_fast_anchor59_raw_scores_logfc", "core_bidirectional_hem--trunk--tabpfn_top500"], ["jeremy_frozen_tabpfn26", "core_bidirectional_hem--trunk--lgbm_l1"], ["jeremy_scores_log2fc_tabpfn26", "ligand_only--trunk--tabpfn_top500"]]`
- cyp1a2, fold 2, full_bank, parameter 0.30000000000000004: `{"chemeleon_pca256_log2fc_tabicl": 0.08, "tabpfn35_fast_pub030_raw": 0.06, "lgbm_jeremy_all_log2fc": 0.24, "lgbm_dense_pool_top256": 0.18, "tabpfn35fast_all_pool_top512": 0.14, "core_bidirectional_hem--trunk--lgbm_l1": 0.08, "core_bidirectional_hem--trunk--tabpfn_top500": 0.14, "ligand_only--trunk--tabpfn_top500": 0.08}`
- cyp1a2, fold 2, original_bank_control, parameter 0.0: `{"jeremy_frozen_tabpfn26": 0.04, "jeremy_scores_log2fc_tabpfn26": 0.02, "jeremy_finetuned_log2fc_lgbm": 0.04, "chemeleon_pca256_log2fc_tabicl": 0.08, "tabpfn35_fast_pub030_raw": 0.16, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.02, "lgbm_all_pool_top512": 0.06, "lgbm_jeremy_all_log2fc": 0.26, "lgbm_dense_pool_top256": 0.16, "tabpfn35fast_all_pool_top512": 0.16}`
- cyp1a2, fold 3, block_added, parameter 0.2: `["core_bidirectional_hem--trunk_molecular--lgbm_l1", "ligand_only--trunk--tabpfn_top500", "core_bidirectional_hem--trunk--tabpfn_top500", "ligand_only--trunk_molecular--tabpfn_top500", "core_bidirectional_hem--trunk--lgbm_l1"]`
- cyp1a2, fold 3, multiple_added, parameter 3.0: `["core_bidirectional_hem--trunk_molecular--lgbm_l1", "ligand_only--trunk--tabpfn_top500", "core_bidirectional_hem--trunk--tabpfn_top500"]`
- cyp1a2, fold 3, multiple_replaced, parameter 2.0: `[["jeremy_frozen_tabpfn26", "ligand_only--trunk--tabpfn_top500"], ["tabpfn35_fast_anchor59_raw_scores_logfc", "core_bidirectional_hem--trunk--tabpfn_top500"]]`
- cyp1a2, fold 3, full_bank, parameter 0.30000000000000004: `{"jeremy_finetuned_log2fc_lgbm": 0.1, "chemeleon_pca256_log2fc_tabicl": 0.1, "lgbm_jeremy_all_log2fc": 0.2, "lgbm_dense_pool_top256": 0.1, "tabpfn35fast_all_pool_top512": 0.2, "core_bidirectional_hem--trunk--lgbm_l1": 0.1, "core_bidirectional_hem--trunk--tabpfn_top500": 0.1, "ligand_only--trunk--tabpfn_top500": 0.1}`
- cyp1a2, fold 3, original_bank_control, parameter 0.0: `{"jeremy_scores_log2fc_tabpfn26": 0.08, "jeremy_finetuned_log2fc_lgbm": 0.16, "chemeleon_pca256_log2fc_tabicl": 0.06, "tabpfn35_fast_pub030_raw": 0.08, "lgbm_jeremy_all_log2fc": 0.22, "lgbm_dense_pool_top256": 0.12, "tabpfn35fast_all_pool_top512": 0.28}`
- cyp1a2, fold 4, block_added, parameter 0.3: `["core_bidirectional_hem--trunk_molecular--lgbm_l1", "core_bidirectional_hem--trunk--tabpfn_top500"]`
- cyp1a2, fold 4, multiple_added, parameter 2.0: `["core_bidirectional_hem--trunk_molecular--lgbm_l1", "core_bidirectional_hem--trunk--tabpfn_top500"]`
- cyp1a2, fold 4, multiple_replaced, parameter 3.0: `[["jeremy_frozen_tabpfn26", "core_bidirectional_hem--trunk--tabpfn_top500"], ["tabpfn35_fast_anchor59_raw_scores_logfc", "core_bidirectional_hem--trunk--lgbm_l1"], ["jeremy_scores_log2fc_tabpfn26", "ligand_only--trunk--tabpfn_top500"]]`
- cyp1a2, fold 4, full_bank, parameter 0.30000000000000004: `{"chemeleon_pca256_log2fc_tabicl": 0.05, "tabpfn35_fast_pub030_raw": 0.1, "lgbm_all_pool_top512": 0.1, "lgbm_jeremy_all_log2fc": 0.1, "lgbm_dense_pool_top256": 0.2, "tabpfn35fast_all_pool_top512": 0.15, "core_bidirectional_hem--trunk--lgbm_l1": 0.1, "core_bidirectional_hem--trunk--tabpfn_top500": 0.15, "ligand_only--trunk--tabpfn_top500": 0.05}`
- cyp1a2, fold 4, original_bank_control, parameter 0.0: `{"jeremy_scores_log2fc_tabpfn26": 0.02, "jeremy_finetuned_log2fc_lgbm": 0.04, "chemeleon_pca256_log2fc_tabicl": 0.08, "tabpfn35_fast_pub030_raw": 0.22, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.04, "lgbm_all_pool_top512": 0.08, "lgbm_jeremy_all_log2fc": 0.18, "lgbm_dense_pool_top256": 0.18, "tabpfn35fast_all_pool_top512": 0.16}`
- cyp2c9, fold 0, block_added, parameter 0.1: `["extended_separate_hem--trunk_molecular--tabpfn_top500", "ligand_only--trunk_molecular--lgbm_l1", "core_distribution_hem--trunk_molecular--tabpfn_top500", "ligand_only--trunk_molecular--tabpfn_top500", "pxr_protein--trunk--lgbm_l1", "pxr_protein--trunk_molecular--tabpfn_top500", "pxr_extended--trunk_molecular--lgbm_l1", "pxr_extended--trunk_molecular--tabpfn_top500", "core_separate_fe--trunk--tabpfn_top500", "extended_mixed_hem--trunk_molecular--lgbm_l1"]`
- cyp2c9, fold 0, multiple_added, parameter 1.0: `["extended_separate_hem--trunk_molecular--tabpfn_top500"]`
- cyp2c9, fold 0, multiple_replaced, parameter 0.0: `[]`
- cyp2c9, fold 0, full_bank, parameter 0.21999999999999997: `{"jeremy_scores_log2fc_tabpfn26": 0.16, "jeremy_finetuned_log2fc_lgbm": 0.26, "chemeleon_pca256_log2fc_tabicl": 0.12, "owned_broad5_2c9_specialist_lgbm": 0.04, "tabpfn35_fast_chemeleon_raw_logfc": 0.02, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.04, "tabpfn35_full_pub030_raw": 0.14, "extended_separate_hem--trunk_molecular--tabpfn_top500": 0.04, "ligand_only--trunk--tabpfn_top500": 0.04, "pxr_extended--trunk_molecular--lgbm_l1": 0.08, "pxr_protein--trunk--lgbm_l1": 0.04, "pxr_protein--trunk--tabpfn_top500": 0.02}`
- cyp2c9, fold 0, original_bank_control, parameter 0.0: `{"jeremy_scores_log2fc_tabpfn26": 0.15, "jeremy_finetuned_log2fc_lgbm": 0.3, "chemeleon_pca256_log2fc_tabicl": 0.15, "owned_broad5_2c9_specialist_lgbm": 0.15, "tabpfn35_fast_chemeleon_raw_logfc": 0.05, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.05, "tabpfn35_full_pub030_raw": 0.15}`
- cyp2c9, fold 1, block_added, parameter 0.05: `["extended_separate_hem--trunk_molecular--tabpfn_top500", "ligand_only--trunk_molecular--lgbm_l1", "core_distribution_hem--trunk_molecular--tabpfn_top500", "pxr_extended--trunk_molecular--lgbm_l1", "ligand_only--trunk_molecular--tabpfn_top500"]`
- cyp2c9, fold 1, multiple_added, parameter 0.0: `[]`
- cyp2c9, fold 1, multiple_replaced, parameter 0.0: `[]`
- cyp2c9, fold 1, full_bank, parameter 0.08: `{"jeremy_scores_log2fc_tabpfn26": 0.12, "jeremy_finetuned_log2fc_lgbm": 0.28, "chemeleon_pca256_log2fc_tabicl": 0.16, "owned_broad5_2c9_specialist_lgbm": 0.12, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.06, "tabpfn35_full_anchor59_raw_scores_logfc": 0.02, "tabpfn35_full_pub030_raw": 0.16, "ligand_only--trunk--lgbm_l1": 0.02, "pxr_core--trunk--tabpfn_top500": 0.06}`
- cyp2c9, fold 1, original_bank_control, parameter 0.0: `{"jeremy_scores_log2fc_tabpfn26": 0.15, "jeremy_finetuned_log2fc_lgbm": 0.3, "chemeleon_pca256_log2fc_tabicl": 0.15, "owned_broad5_2c9_specialist_lgbm": 0.2, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.05, "tabpfn35_full_anchor59_raw_scores_logfc": 0.05, "tabpfn35_full_pub030_raw": 0.1}`
- cyp2c9, fold 2, block_added, parameter 0.1: `["extended_separate_hem--trunk_molecular--tabpfn_top500"]`
- cyp2c9, fold 2, multiple_added, parameter 1.0: `["extended_separate_hem--trunk_molecular--tabpfn_top500"]`
- cyp2c9, fold 2, multiple_replaced, parameter 0.0: `[]`
- cyp2c9, fold 2, full_bank, parameter 0.1: `{"jeremy_scores_log2fc_tabpfn26": 0.1, "jeremy_finetuned_log2fc_lgbm": 0.3, "chemeleon_pca256_log2fc_tabicl": 0.2, "owned_broad5_2c9_specialist_lgbm": 0.1, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.1, "tabpfn35_full_pub030_raw": 0.1, "core_separate_hem_s--trunk--tabpfn_top500": 0.1}`
- cyp2c9, fold 2, original_bank_control, parameter 0.0: `{"jeremy_scores_log2fc_tabpfn26": 0.12, "jeremy_finetuned_log2fc_lgbm": 0.28, "chemeleon_pca256_log2fc_tabicl": 0.22, "owned_broad5_2c9_specialist_lgbm": 0.12, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.1, "tabpfn35_full_pub030_raw": 0.16}`
- cyp2c9, fold 3, block_added, parameter 0.1: `["extended_separate_hem--trunk_molecular--tabpfn_top500"]`
- cyp2c9, fold 3, multiple_added, parameter 1.0: `["extended_separate_hem--trunk_molecular--tabpfn_top500"]`
- cyp2c9, fold 3, multiple_replaced, parameter 0.0: `[]`
- cyp2c9, fold 3, full_bank, parameter 0.18: `{"jeremy_scores_log2fc_tabpfn26": 0.14, "jeremy_finetuned_log2fc_lgbm": 0.24, "chemeleon_pca256_log2fc_tabicl": 0.22, "owned_broad5_2c9_specialist_lgbm": 0.08, "tabpfn35_full_chemeleon_raw_logfc": 0.02, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.04, "tabpfn35_full_pub030_raw": 0.08, "ligand_only--trunk--tabpfn_top500": 0.04, "pxr_core--trunk--tabpfn_top500": 0.08, "pxr_protein--trunk--lgbm_l1": 0.06}`
- cyp2c9, fold 3, original_bank_control, parameter 0.0: `{"jeremy_scores_log2fc_tabpfn26": 0.15, "jeremy_finetuned_log2fc_lgbm": 0.25, "chemeleon_pca256_log2fc_tabicl": 0.25, "owned_broad5_2c9_specialist_lgbm": 0.15, "tabpfn35_full_chemeleon_raw_logfc": 0.05, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.05, "tabpfn35_full_pub030_raw": 0.1}`
- cyp2c9, fold 4, block_added, parameter 0.1: `["extended_separate_hem--trunk_molecular--tabpfn_top500"]`
- cyp2c9, fold 4, multiple_added, parameter 0.0: `[]`
- cyp2c9, fold 4, multiple_replaced, parameter 0.0: `[]`
- cyp2c9, fold 4, full_bank, parameter 0.15000000000000002: `{"jeremy_scores_log2fc_tabpfn26": 0.2, "jeremy_finetuned_log2fc_lgbm": 0.3, "chemeleon_pca256_log2fc_tabicl": 0.15, "owned_broad5_2c9_specialist_lgbm": 0.05, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.05, "tabpfn35_full_pub030_raw": 0.1, "pxr_core--trunk--tabpfn_top500": 0.05, "pxr_protein--trunk--lgbm_l1": 0.1}`
- cyp2c9, fold 4, original_bank_control, parameter 0.0: `{"jeremy_scores_log2fc_tabpfn26": 0.2, "jeremy_finetuned_log2fc_lgbm": 0.3, "chemeleon_pca256_log2fc_tabicl": 0.2, "owned_broad5_2c9_specialist_lgbm": 0.15, "tabpfn35_full_anchor59_raw_scores_logfc": 0.05, "tabpfn35_full_pub030_raw": 0.1}`
- cyp3a4, fold 0, block_added, parameter 0.05: `["pxr_extended--trunk_molecular--tabpfn_top500", "ligand_only--trunk_molecular--tabpfn_top500", "extended_separate_hem--trunk_molecular--lgbm_l1", "core_separate_fe--trunk_molecular--tabpfn_top500", "ligand_only--trunk--lgbm_l1"]`
- cyp3a4, fold 0, multiple_added, parameter 0.0: `[]`
- cyp3a4, fold 0, multiple_replaced, parameter 0.0: `[]`
- cyp3a4, fold 0, full_bank, parameter 0.05: `{"class_quantile_065_lgbm": 0.1, "morgan_quantile_065_lgbm": 0.1, "jeremy_scores_log2fc_tabpfn26": 0.1, "jeremy_finetuned_log2fc_lgbm": 0.1, "chemeleon_pca256_log2fc_tabicl": 0.05, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.05, "tabpfn35_full_anchor59_raw_scores_logfc": 0.15, "tdi_continuous_arm_lgbm": 0.25, "tdi_continuous_delta_lgbm": 0.05, "ligand_only--trunk--lgbm_l1": 0.05}`
- cyp3a4, fold 0, original_bank_control, parameter 0.0: `{"class_quantile_065_lgbm": 0.1, "morgan_quantile_065_lgbm": 0.1, "jeremy_scores_log2fc_tabpfn26": 0.1, "jeremy_finetuned_log2fc_lgbm": 0.1, "chemeleon_pca256_log2fc_tabicl": 0.1, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.1, "tabpfn35_full_anchor59_raw_scores_logfc": 0.1, "tdi_continuous_arm_lgbm": 0.2, "tdi_continuous_delta_lgbm": 0.1}`
- cyp3a4, fold 1, block_added, parameter 0.0: `["pxr_extended--trunk_molecular--tabpfn_top500"]`
- cyp3a4, fold 1, multiple_added, parameter 0.0: `[]`
- cyp3a4, fold 1, multiple_replaced, parameter 0.0: `[]`
- cyp3a4, fold 1, full_bank, parameter 0.04: `{"class_quantile_065_lgbm": 0.08, "morgan_quantile_065_lgbm": 0.14, "jeremy_finetuned_log2fc_lgbm": 0.1, "chemeleon_pca256_log2fc_tabicl": 0.12, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.16, "tabpfn35_full_anchor59_raw_scores_logfc": 0.08, "tdi_continuous_arm_lgbm": 0.24, "tdi_continuous_delta_lgbm": 0.04, "ligand_only--trunk--lgbm_l1": 0.04}`
- cyp3a4, fold 1, original_bank_control, parameter 0.0: `{"class_quantile_065_lgbm": 0.08, "morgan_quantile_065_lgbm": 0.16, "jeremy_finetuned_log2fc_lgbm": 0.1, "chemeleon_pca256_log2fc_tabicl": 0.12, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.14, "tabpfn35_full_anchor59_raw_scores_logfc": 0.12, "tdi_continuous_arm_lgbm": 0.24, "tdi_continuous_delta_lgbm": 0.04}`
- cyp3a4, fold 2, block_added, parameter 0.0: `["pxr_extended--trunk_molecular--tabpfn_top500"]`
- cyp3a4, fold 2, multiple_added, parameter 0.0: `[]`
- cyp3a4, fold 2, multiple_replaced, parameter 0.0: `[]`
- cyp3a4, fold 2, full_bank, parameter 0.0: `{"class_quantile_065_lgbm": 0.1, "morgan_quantile_065_lgbm": 0.1, "jeremy_finetuned_log2fc_lgbm": 0.2, "chemeleon_pca256_log2fc_tabicl": 0.1, "tabpfn35_full_anchor59_raw_scores_logfc": 0.3, "tdi_continuous_arm_lgbm": 0.2}`
- cyp3a4, fold 2, original_bank_control, parameter 0.0: `{"class_quantile_065_lgbm": 0.1, "morgan_quantile_065_lgbm": 0.1, "jeremy_finetuned_log2fc_lgbm": 0.2, "chemeleon_pca256_log2fc_tabicl": 0.1, "tabpfn35_full_anchor59_raw_scores_logfc": 0.3, "tdi_continuous_arm_lgbm": 0.2}`
- cyp3a4, fold 3, block_added, parameter 0.0: `["pxr_extended--trunk_molecular--tabpfn_top500"]`
- cyp3a4, fold 3, multiple_added, parameter 0.0: `[]`
- cyp3a4, fold 3, multiple_replaced, parameter 0.0: `[]`
- cyp3a4, fold 3, full_bank, parameter 0.1: `{"class_quantile_065_lgbm": 0.1, "morgan_quantile_065_lgbm": 0.1, "jeremy_finetuned_log2fc_lgbm": 0.1, "chemeleon_pca256_log2fc_tabicl": 0.1, "tabpfn35_full_anchor59_raw_scores_logfc": 0.3, "tdi_continuous_arm_lgbm": 0.2, "ligand_only--trunk--lgbm_l1": 0.1}`
- cyp3a4, fold 3, original_bank_control, parameter 0.0: `{"class_quantile_065_lgbm": 0.1, "morgan_quantile_065_lgbm": 0.1, "jeremy_finetuned_log2fc_lgbm": 0.1, "chemeleon_pca256_log2fc_tabicl": 0.1, "tabpfn35_full_anchor59_raw_scores_logfc": 0.4, "tdi_continuous_arm_lgbm": 0.2}`
- cyp3a4, fold 4, block_added, parameter 0.0: `["pxr_extended--trunk_molecular--tabpfn_top500"]`
- cyp3a4, fold 4, multiple_added, parameter 0.0: `[]`
- cyp3a4, fold 4, multiple_replaced, parameter 0.0: `[]`
- cyp3a4, fold 4, full_bank, parameter 0.0: `{"morgan_quantile_065_lgbm": 0.2, "jeremy_finetuned_log2fc_lgbm": 0.2, "tabpfn35_full_anchor59_raw_scores_logfc": 0.4, "tdi_continuous_arm_lgbm": 0.2}`
- cyp3a4, fold 4, original_bank_control, parameter 0.0: `{"class_quantile_065_lgbm": 0.1, "morgan_quantile_065_lgbm": 0.1, "jeremy_finetuned_log2fc_lgbm": 0.1, "chemeleon_pca256_log2fc_tabicl": 0.1, "tabpfn35_fast_anchor59_raw_scores_logfc": 0.1, "tabpfn35_full_anchor59_raw_scores_logfc": 0.3, "tdi_continuous_arm_lgbm": 0.2}`
- cyp2d6, fold 0, block_added, parameter 0.2: `["pxr_core--trunk_molecular--lgbm_l1", "core_separate_hem--trunk--tabpfn_top500", "core_distribution_hem--trunk_molecular--lgbm_l1", "ligand_only--trunk--tabpfn_top500", "core_separate_fe--trunk_molecular--lgbm_l1", "core_distribution_hem--trunk_molecular--tabpfn_top500", "ligand_only--trunk--lgbm_l1", "core_bidirectional_hem--trunk--lgbm_l1", "core_separate_hem--trunk_molecular--tabpfn_top500", "core_distribution_hem--trunk--lgbm_l1"]`
- cyp2d6, fold 0, multiple_added, parameter 1.0: `["pxr_core--trunk_molecular--lgbm_l1"]`
- cyp2d6, fold 0, multiple_replaced, parameter 0.0: `[]`
- cyp2d6, fold 0, full_bank, parameter 0.30000000000000004: `{"submitted_two_member_core": 0.3, "moment_matched_robust_branch": 0.4, "core_distribution_hem--trunk_molecular--lgbm_l1": 0.1, "core_separate_hem--trunk--tabpfn_top500": 0.2}`
- cyp2d6, fold 0, original_bank_control, parameter 0.0: `{"submitted_two_member_core": 0.6, "moment_matched_robust_branch": 0.4}`
- cyp2d6, fold 1, block_added, parameter 0.2: `["pxr_core--trunk_molecular--lgbm_l1"]`
- cyp2d6, fold 1, multiple_added, parameter 1.0: `["pxr_core--trunk_molecular--lgbm_l1"]`
- cyp2d6, fold 1, multiple_replaced, parameter 0.0: `[]`
- cyp2d6, fold 1, full_bank, parameter 0.2: `{"submitted_two_member_core": 0.4, "moment_matched_robust_branch": 0.4, "core_distribution_hem--trunk_molecular--lgbm_l1": 0.2}`
- cyp2d6, fold 1, original_bank_control, parameter 0.0: `{"submitted_two_member_core": 0.5, "moment_matched_robust_branch": 0.5}`
- cyp2d6, fold 2, block_added, parameter 0.2: `["core_distribution_hem--trunk_molecular--lgbm_l1"]`
- cyp2d6, fold 2, multiple_added, parameter 1.0: `["core_distribution_hem--trunk_molecular--lgbm_l1"]`
- cyp2d6, fold 2, multiple_replaced, parameter 0.0: `[]`
- cyp2d6, fold 2, full_bank, parameter 0.3: `{"submitted_two_member_core": 0.2, "moment_matched_robust_branch": 0.5, "core_distribution_hem--trunk_molecular--lgbm_l1": 0.15, "core_separate_hem--trunk--tabpfn_top500": 0.15}`
- cyp2d6, fold 2, original_bank_control, parameter 0.0: `{"submitted_two_member_core": 0.5, "moment_matched_robust_branch": 0.5}`
- cyp2d6, fold 3, block_added, parameter 0.2: `["core_separate_fe--trunk_molecular--lgbm_l1", "core_distribution_hem--trunk_molecular--lgbm_l1", "core_separate_hem--trunk--tabpfn_top500", "ligand_only--trunk--lgbm_l1", "core_distribution_hem--trunk_molecular--tabpfn_top500", "core_bidirectional_hem--trunk--lgbm_l1", "pxr_core--trunk_molecular--lgbm_l1", "core_distribution_hem--trunk--lgbm_l1", "ligand_only--trunk_molecular--tabpfn_top500", "extended_separate_hem--trunk--lgbm_l1"]`
- cyp2d6, fold 3, multiple_added, parameter 1.0: `["core_separate_fe--trunk_molecular--lgbm_l1"]`
- cyp2d6, fold 3, multiple_replaced, parameter 0.0: `[]`
- cyp2d6, fold 3, full_bank, parameter 0.30000000000000004: `{"submitted_two_member_core": 0.3, "moment_matched_robust_branch": 0.4, "core_bidirectional_hem--trunk--lgbm_l1": 0.05, "core_distribution_hem--trunk--lgbm_l1": 0.1, "core_distribution_hem--trunk_molecular--lgbm_l1": 0.05, "core_separate_hem--trunk--tabpfn_top500": 0.1}`
- cyp2d6, fold 3, original_bank_control, parameter 0.0: `{"submitted_two_member_core": 0.6, "moment_matched_robust_branch": 0.4}`
- cyp2d6, fold 4, block_added, parameter 0.3: `["pxr_core--trunk_molecular--lgbm_l1", "core_separate_hem--trunk--tabpfn_top500", "core_distribution_hem--trunk_molecular--lgbm_l1", "pxr_extended--trunk--lgbm_l1", "ligand_only--trunk_molecular--tabpfn_top500", "core_distribution_hem--trunk--lgbm_l1", "core_bidirectional_hem--trunk--lgbm_l1", "core_distribution_hem--trunk_molecular--tabpfn_top500", "ligand_only--trunk--lgbm_l1", "extended_separate_hem--trunk--lgbm_l1"]`
- cyp2d6, fold 4, multiple_added, parameter 1.0: `["pxr_core--trunk_molecular--lgbm_l1"]`
- cyp2d6, fold 4, multiple_replaced, parameter 0.0: `[]`
- cyp2d6, fold 4, full_bank, parameter 0.35000000000000003: `{"submitted_two_member_core": 0.2, "moment_matched_robust_branch": 0.45, "core_distribution_hem--trunk--lgbm_l1": 0.1, "core_distribution_hem--trunk_molecular--lgbm_l1": 0.05, "core_separate_hem--trunk--tabpfn_top500": 0.15, "pxr_extended--trunk--lgbm_l1": 0.05}`
- cyp2d6, fold 4, original_bank_control, parameter 0.0: `{"submitted_two_member_core": 0.48, "moment_matched_robust_branch": 0.52}`

## Interpretation and limits

1A2 multiple replacement gives 0.665644 versus current 0.670833 and previous single replacement 0.669032. Trunk aggregation gives smaller gains. D6 full-bank selection gives 0.815874, comparable to the previous single-trunk extension 0.815745; extra trunks do not establish further improvement. C9 changes are negligible or adverse; A4 has no benefit. Trunk-only aggregation remains worse than the incumbent. Correlated errors limit the benefit of adding similar views.

This is conditional saved-OOF selection, not fully nested base-model refitting. Base predictions in other folds can have been trained using held-fold labels. No test-set prediction, deployment choice or submission was produced.
