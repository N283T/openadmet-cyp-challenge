# Global Caruana across saved regression model families

The user requested actual selection across all available compatible model families, replacing the prior manual top-level average.

Selection run: `20260912T091825933326Z-1548d56b`. Independent audit: `20260912T091833613985Z-10936af9`.

The frozen inventory contains 159 completed OOF artifacts. Exact endpoint/fold/label alignment and duplicate-vector removal leave 107 candidate identities (100–105 eligible per endpoint).
Candidates include LightGBM, TabPFN, TabICL, Jazzy additions, FP/Mordred, external representations, neural direct predictions, and previously nested Caruana components. Per-fold early FP/Mordred and separate embedding/log2fc readouts are assembled explicitly. The pool also retains prior fixed-ensemble predictions exposed by independent final-package verification runs (m0079 and m0080), alongside individual model candidates. Thus the pool contains both individual models and aggregate predictions; it is not exclusively atomic models. Calibration and other comparison scripts are excluded by the recorded inventory filter, but this is not a claim that every earlier manual blend was excluded. Exact duplicate vectors are collapsed. All source exclusions and duplicate aliases are saved.

For each endpoint, greedy selection with replacement runs 50 steps and retains the best training-loss prefix. Interval MAE is minimized; the ST-RAE denominator is fixed within the endpoint, so this optimizes the official ST-RAE objective. No candidate score prefilter, new hyperparameter search, or clipping is used.

Two outputs are distinct: full-OOF fitted weights (development fit), and leave-one-parent-fold-out weight fitting. The latter excludes evaluation-fold labels from weight fitting only. Base models were not retrained in newly nested splits, and some selector-training predictions come from models trained using that evaluation fold. Therefore it is a weight-crossfit sensitivity check, not an independent nested-CV estimate.

## Scores

| Endpoint | Previous candidate | Weight cross-fit | Full-OOF weight fit |
|---|---:|---:|---:|
| CYP1A2 | 0.670809 | 0.674078 | 0.663594 |
| CYP2C9 | 0.452634 | 0.442189 | 0.435268 |
| CYP2D6 | 0.822681 | 0.825782 | 0.810934 |
| CYP3A4 | 0.401599 | 0.390980 | 0.385022 |
| Macro | 0.586931 | 0.583257 | 0.573705 |

Global weight selection improves the macro point estimate chiefly through 2C9 and 3A4. The weight-crossfit results slightly worsen both focus endpoints, 1A2 and 2D6. This is not evidence that their underlying model problem has been solved. Neighbor macro is also essentially unchanged/slightly worse (0.675304 to 0.676179). MAE and Spearman are preserved in the audit and do not improve uniformly.

Retain the global weights, aligned pool, selection histories, correlations, and keyed predictions as a new ensemble experiment. The lower full-fit score must not be presented as an untouched CV gain. No production pointer, full-data test refit, or submission was changed.

## Full-OOF fitted composition

Weights below are the full-OOF selection recipe, not a new independently validated deployment. Exact cross-fitted weights are in `selected-models.csv` and `weights.json`.

### CYP1A2

| Candidate | Weight | Source |
|---|---:|---|
| `m0096` / `20260912T083421950497Z-036b9595` | 26.53% | focus_ft_primary_fixedparams_v1 — Replace only primary features in selected FT600 LGBM component; fixed old parameters, exact exclusions and25/50blends |
| `m0097` / `20260912T085800716037Z-db5687c1` | 14.29% | jazzy_mean — Controlled DG Jazzy additions to frozen600 and new nested log2fc |
| `m0074` / `20260912T072250002510Z-eb7d9b85` | 10.20% | jeremy_ffn_scores_logfc_lgbm_v1 — Matched fixed-CV evaluation of external representation jeremy_ffn_scores_logfc_lgbm_v1 |
| `m0061` / `20260912T062500451470Z-ad87df0a` | 8.16% | external3_class_emax_v1_ft_tabpfn — Fold-specific adapted external600 readout: external3_class_emax_v1_ft_tabpfn |
| `m0094` / `20260912T083210446408Z-8d59a1f4` | 8.16% | focus_cheme_new_primary_tabicl_v1 — Fixed parent CV TabICL: embedding-only training PCA, then nested log2fc; fit and save all models: focus_cheme_new_primary_tabicl_v1 |
| `m0099` / `20260912T090040918188Z-695fb6ee` | 8.16% | jazzy_mean — Controlled DG Jazzy additions to frozen600 and new nested log2fc |
| `m0042` / `20260912T060513219458Z-97152a22` | 6.12% | mae_st_selection_smoke — Short MAE and official inner ST-RAE FT checkpoint-selection smoke across five fixed folds |
| `m0032` / `20260912T054406239560Z-96fce32b` | 4.08% | nested_caruana_external_lgbm — Honest nested Caruana selection among24 external feature LGBMs using inner group holdout only |
| `m0079` / `20260912T073002574503Z-8be77a87` | 4.08% | caruana_medium_scores_ft_cheme_equal_v1 — Independent final fixed-five OOF reconstruction, official scoring and provenance package |
| `m0092` / `20260912T083032513346Z-568bc526` | 4.08% | focus_primary_new_tabpfn_v1 — Matched fixed-CV evaluation of external representation focus_primary_new_tabpfn_v1 |
| `m0030` / `20260912T054054898981Z-8bde71ca` | 2.04% | external_measured_masked_logfc_lgbm — Matched fixed-CV evaluation of external representation external_measured_masked_logfc_lgbm |
| `m0035` / `20260912T055130250735Z-8bce3e59` | 2.04% | external3_cond_v1_mixed_lgbm — Matched fixed-CV evaluation of external representation external3_cond_v1_mixed_lgbm |
| `m0053` / `20260912T061514330688Z-3ff745e1` | 2.04% | external_az_mixed_v1_mixed_adapted_ft_tabpfn — Fold-specific adapted external600 readout: external_az_mixed_v1_mixed_adapted_ft_tabpfn |

### CYP2C9

| Candidate | Weight | Source |
|---|---:|---|
| `m0094` / `20260912T083210446408Z-8d59a1f4` | 40.43% | focus_cheme_new_primary_tabicl_v1 — Fixed parent CV TabICL: embedding-only training PCA, then nested log2fc; fit and save all models: focus_cheme_new_primary_tabicl_v1 |
| `m0096` / `20260912T083421950497Z-036b9595` | 17.02% | focus_ft_primary_fixedparams_v1 — Replace only primary features in selected FT600 LGBM component; fixed old parameters, exact exclusions and25/50blends |
| `m0072` / `20260912T072056445264Z-e6b6cbcb` | 12.77% | jeremy_ft_logfc_rustuna_v1 — Inner-only Rustuna LGBM tuning of outer-safe Jeremy FT600 and log2fc4 |
| `m0058` / `20260912T062241730532Z-f11cd0a6` | 6.38% | jeremy_medium_mt_pca256_logfc_tabpfn_v1 — Readout of completed fold-specific FT embeddings: jeremy_medium_mt_pca256_logfc_tabpfn_v1 |
| `m0021` / `20260912T053152207404Z-f4f0dfec` | 4.26% | external_tree_logfc_lgbm — Matched fixed-CV evaluation of external representation external_tree_logfc_lgbm |
| `m0071` / `20260912T072034139020Z-92ad0df3` | 4.26% | jeremy_interval_smoke_v1 — Fivefold2epoch smoke officialinterval pIC50 loss plusMSElogfc; strictinnerconfidence joins andfinitegradients |
| `m0077` / `20260912T072355846776Z-2e678154` | 4.26% | jeremy_medium_interval160_v1_ft_direct — External600 -> pIC50/log2fc and configured auxiliary heads; outer compound exclusion: jeremy_medium_interval160_v1 |
| `m0002` / `20260911T111952347762Z-e9f6c457` | 2.13% | lightgbm_full — Compare TabICL compact/full and LightGBM on corrected parent-retained CV |
| `m0005` / `20260911T150557650787Z-4a460d9e` | 2.13% | admetica_emb_tabpfn — Evaluate frozen per-CYP Admetica embeddings alone with TabPFN 2.6 on matched folds |
| `m0017` / `20260912T052601654102Z-6844f96d` | 2.13% | external3_v1_continuous_lgbm — Matched fixed-CV evaluation of external representation external3_v1_continuous_lgbm |
| `m0080` / `20260912T073049029140Z-0c26dd06` | 2.13% | caruana_medium_scores_tuned_ft_cheme_v1 — Independent final fixed-five OOF reconstruction, official scoring and provenance package |
| `m0104` / `20260912T062329198020Z-60f59445` | 2.13% | jeremy_classic_expanded_v1 — Nested Caruana inner-only selection of external LGBMs with Jeremy and train-filtered FP+Mordred candidates |

### CYP2D6

| Candidate | Weight | Source |
|---|---:|---|
| `m0093` / `20260912T083209287272Z-e8487aeb` | 22.00% | focus_scores_new_primary_tabpfn_v1 — Matched fixed-CV evaluation of external representation focus_scores_new_primary_tabpfn_v1 |
| `m0096` / `20260912T083421950497Z-036b9595` | 18.00% | focus_ft_primary_fixedparams_v1 — Replace only primary features in selected FT600 LGBM component; fixed old parameters, exact exclusions and25/50blends |
| `m0106` / `20260912T082618415537Z-79e86a45` | 14.00% | focus_single_CYP2D6_lgbm_logfc — focus_single_CYP2D6_lgbm_logfc_v1 |
| `m0011` / `20260912T051258883571Z-a2ab5b73` | 12.00% | pubchem600_scores_lgbm — Matched fixed-CV evaluation of external representation pubchem600_scores_lgbm |
| `m0034` / `20260912T054534754156Z-3c369a92` | 6.00% | external3_cond_v1_classification_tabpfn — Matched fixed-CV evaluation of external representation external3_cond_v1_classification_tabpfn |
| `m0074` / `20260912T072250002510Z-eb7d9b85` | 6.00% | jeremy_ffn_scores_logfc_lgbm_v1 — Matched fixed-CV evaluation of external representation jeremy_ffn_scores_logfc_lgbm_v1 |
| `m0076` / `20260912T072352902207Z-1b836eeb` | 6.00% | jeremy_ffn_hidden_logfc_tabpfn_v1 — Matched fixed-CV evaluation of external representation jeremy_ffn_hidden_logfc_tabpfn_v1 |
| `m0085` / `20260912T082257132357Z-bf656b4e` | 6.00% | combined_lgbm — Focused1A2/2D6: grouped inner-selected RBF-SVR and combined external features LGBM |
| `m0098` / `20260912T085800716037Z-db5687c1` | 4.00% | jazzy_mean_std — Controlled DG Jazzy additions to frozen600 and new nested log2fc |
| `m0036` / `20260912T055149151415Z-6c457dd0` | 2.00% | external3_cond_v1_mixed_tabpfn — Matched fixed-CV evaluation of external representation external3_cond_v1_mixed_tabpfn |
| `m0039` / `20260912T055825039246Z-36aa5478` | 2.00% | external3_class_mt_long_v1_ft_tabpfn — Fold-specific adapted external600 readout: external3_class_mt_long_v1_ft_tabpfn |
| `m0079` / `20260912T073002574503Z-8be77a87` | 2.00% | caruana_medium_scores_ft_cheme_equal_v1 — Independent final fixed-five OOF reconstruction, official scoring and provenance package |

### CYP3A4

| Candidate | Weight | Source |
|---|---:|---|
| `m0077` / `20260912T072355846776Z-2e678154` | 32.65% | jeremy_medium_interval160_v1_ft_direct — External600 -> pIC50/log2fc and configured auxiliary heads; outer compound exclusion: jeremy_medium_interval160_v1 |
| `m0096` / `20260912T083421950497Z-036b9595` | 16.33% | focus_ft_primary_fixedparams_v1 — Replace only primary features in selected FT600 LGBM component; fixed old parameters, exact exclusions and25/50blends |
| `m0093` / `20260912T083209287272Z-e8487aeb` | 8.16% | focus_scores_new_primary_tabpfn_v1 — Matched fixed-CV evaluation of external representation focus_scores_new_primary_tabpfn_v1 |
| `m0071` / `20260912T072034139020Z-92ad0df3` | 6.12% | jeremy_interval_smoke_v1 — Fivefold2epoch smoke officialinterval pIC50 loss plusMSElogfc; strictinnerconfidence joins andfinitegradients |
| `m0094` / `20260912T083210446408Z-8d59a1f4` | 6.12% | focus_cheme_new_primary_tabicl_v1 — Fixed parent CV TabICL: embedding-only training PCA, then nested log2fc; fit and save all models: focus_cheme_new_primary_tabicl_v1 |
| `m0099` / `20260912T090040918188Z-695fb6ee` | 6.12% | jazzy_mean — Controlled DG Jazzy additions to frozen600 and new nested log2fc |
| `m0032` / `20260912T054406239560Z-96fce32b` | 4.08% | nested_caruana_external_lgbm — Honest nested Caruana selection among24 external feature LGBMs using inner group holdout only |
| `m0070` / `20260912T071619253798Z-73a6fb7e` | 4.08% | jeremy_chemeleon_logfc_tabicl_v1 — Fixed parent CV TabICL: embedding-only training PCA, then nested log2fc; fit and save all models: jeremy_chemeleon_logfc_tabicl_v1 |
| `m0080` / `20260912T073049029140Z-0c26dd06` | 4.08% | caruana_medium_scores_tuned_ft_cheme_v1 — Independent final fixed-five OOF reconstruction, official scoring and provenance package |
| `m0002` / `20260911T111952347762Z-e9f6c457` | 2.04% | lightgbm_full — Compare TabICL compact/full and LightGBM on corrected parent-retained CV |
| `m0030` / `20260912T054054898981Z-8bde71ca` | 2.04% | external_measured_masked_logfc_lgbm — Matched fixed-CV evaluation of external representation external_measured_masked_logfc_lgbm |
| `m0058` / `20260912T062241730532Z-f11cd0a6` | 2.04% | jeremy_medium_mt_pca256_logfc_tabpfn_v1 — Readout of completed fold-specific FT embeddings: jeremy_medium_mt_pca256_logfc_tabpfn_v1 |
| `m0072` / `20260912T072056445264Z-e6b6cbcb` | 2.04% | jeremy_ft_logfc_rustuna_v1 — Inner-only Rustuna LGBM tuning of outer-safe Jeremy FT600 and log2fc4 |
| `m0075` / `20260912T072303592078Z-07f2749a` | 2.04% | jeremy_ffn_scores_logfc_tabpfn_v1 — Matched fixed-CV evaluation of external representation jeremy_ffn_scores_logfc_tabpfn_v1 |
| `m0100` / `20260912T090040918188Z-695fb6ee` | 2.04% | jazzy_mean_std — Controlled DG Jazzy additions to frozen600 and new nested log2fc |

## Verification

The independent audit reconstructs all predictions from nonnegative, sum-one weights; verifies 6,384 exact keys and original confidence bounds; checks all 24 weight sets and selection-fold masks; and recomputes official pooled, fold, and neighbor metrics. Maximum reconstruction error is below 2e-15. The previous 0.5869306288697727 score is reproduced.

Saved artifacts: `catalog.json`, `excluded.json`, `aliases.json`, `aligned-pool.parquet`, `global-oof.parquet`, `weights.json`, `selection-history.json`, `selected-models.csv`, and independent `metrics.csv`/`macro.csv`.

Readable report and selected-model CSV: `20260912T092320697419Z-ab0aeab3`.
This corrects the pool description in initial summary `bdfc7c1d`, which incorrectly
said prior manual blends were excluded. Selection, weights and scores are unchanged;
both reports remain recorded. `pixi run check` passes; tests: 37 passed, 7 skipped.
