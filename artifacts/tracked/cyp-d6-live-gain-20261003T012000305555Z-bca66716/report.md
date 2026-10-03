# Sizing the 2026-10-03 D6 Live gain

Scored subset estimate: 374 of 750 rows, 26 scored positives (prevalence 0.0695).

## Part 1 - operating-point estimators

| state | submitted | calls | mcc | positives_from_recall_sd | calls_from_precision_sd | scored_from_accuracy_sd | prevalence_identity | call_fraction_identity | scored_from_identity | positives_from_identity |
|---|---|---|---|---|---|---|---|---|---|---|
| sep13 | 2026-09-13 11:38 UTC | 49 | 0.2710 | 26.3109 | 29.4392 | 375.9029 | 0.0706 | 0.0795 | 370.1555 | 26.1391 |
| sep14 | 2026-09-14 02:08 UTC | 49 | 0.2786 | 26.3109 | 28.3563 | 369.8419 | 0.0706 | 0.0766 | 369.9458 | 26.1338 |
| sep15 | 2026-09-15 13:00 UTC | 49 | 0.3086 | 25.8792 | 29.3650 | 369.4858 | 0.0706 | 0.0794 | 370.0093 | 26.1272 |
| sep18 | 2026-09-18 09:45 UTC | 49 | 0.1097 | 24.9925 | 28.7009 | 349.6305 | 0.0706 | 0.0817 | 351.2525 | 24.8153 |
| sep23 | 2026-09-23 18:03 UTC | 49 | 0.2348 | 24.5591 | 26.7373 | 366.6731 | 0.0706 | 0.0769 | 347.5368 | 24.5243 |
| oct01 | 2026-10-01 14:13 UTC | 49 | 0.1050 | 23.2511 | 29.8257 | 393.5680 | 0.0703 | 0.0852 | 350.1033 | 24.6243 |
| oct02 | 2026-10-02 12:08 UTC | 175 | 0.1177 | 26.0387 | 90.0333 | 386.0975 | 0.0706 | 0.2387 | 377.1845 | 26.6316 |
| oct03 | 2026-10-03 00:29 UTC | 49 | 0.3642 | 26.3252 | 27.4080 | 402.5987 | 0.0706 | 0.0740 | 370.5997 | 26.1655 |

## Part 1 - joint integer fit

| state | total_calls | scored_calls | true_positive | false_positive | false_negative | fit_loss | mcc_reported | mcc_simulated | mcc_sd_reported | mcc_sd_simulated | precision_reported | precision_simulated | recall_reported | recall_simulated |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sep13 | 49 | 30 | 9 | 21 | 17 | 0.0114 | 0.2710 | 0.2668 | 0.0846 | 0.0856 | 0.3080 | 0.3003 | 0.3469 | 0.3472 |
| sep14 | 49 | 29 | 9 | 20 | 17 | 0.0157 | 0.2786 | 0.2722 | 0.0855 | 0.0864 | 0.3197 | 0.3096 | 0.3469 | 0.3454 |
| sep15 | 49 | 30 | 10 | 20 | 16 | 0.0098 | 0.3086 | 0.3043 | 0.0864 | 0.0866 | 0.3412 | 0.3329 | 0.3835 | 0.3843 |
| sep18 | 49 | 31 | 5 | 26 | 21 | 0.0045 | 0.1097 | 0.1078 | 0.0737 | 0.0721 | 0.1655 | 0.1613 | 0.1914 | 0.1924 |
| sep23 | 49 | 29 | 8 | 21 | 18 | 0.0039 | 0.2348 | 0.2343 | 0.0872 | 0.0846 | 0.2800 | 0.2763 | 0.3052 | 0.3080 |
| oct01 | 49 | 32 | 5 | 27 | 21 | 0.0018 | 0.1050 | 0.1040 | 0.0731 | 0.0715 | 0.1588 | 0.1565 | 0.1924 | 0.1927 |
| oct02 | 175 | 89 | 11 | 78 | 15 | 0.0041 | 0.1177 | 0.1184 | 0.0599 | 0.0604 | 0.1247 | 0.1236 | 0.4215 | 0.4234 |
| oct03 | 49 | 28 | 11 | 17 | 15 | 0.0129 | 0.3642 | 0.3603 | 0.0875 | 0.0896 | 0.4024 | 0.3931 | 0.4215 | 0.4235 |

Alternative (scored, positives) pairs by total loss:

| rank | scored | positives | loss | tp_sep13 | tp_sep14 | tp_sep15 | tp_sep18 | tp_sep23 | tp_oct01 | tp_oct02 | tp_oct03 | calls_sep13 | calls_sep14 | calls_sep15 | calls_sep18 | calls_sep23 | calls_oct01 | calls_oct02 | calls_oct03 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 374 | 26 | 0.0639 | 9 | 9 | 10 | 5 | 8 | 5 | 11 | 11 | 30 | 29 | 30 | 31 | 29 | 32 | 89 | 28 |
| 2 | 376 | 26 | 0.0644 | 9 | 9 | 10 | 5 | 8 | 5 | 11 | 11 | 30 | 29 | 30 | 31 | 29 | 32 | 90 | 28 |
| 3 | 375 | 26 | 0.0647 | 9 | 9 | 10 | 5 | 8 | 5 | 11 | 11 | 30 | 29 | 30 | 31 | 29 | 32 | 90 | 28 |
| 4 | 373 | 26 | 0.0693 | 9 | 9 | 10 | 5 | 8 | 5 | 11 | 11 | 30 | 29 | 30 | 31 | 29 | 32 | 89 | 28 |
| 5 | 367 | 26 | 0.0695 | 9 | 9 | 10 | 5 | 8 | 5 | 11 | 11 | 29 | 28 | 29 | 30 | 29 | 31 | 88 | 27 |
| 6 | 377 | 26 | 0.0705 | 9 | 9 | 10 | 5 | 8 | 5 | 11 | 11 | 30 | 29 | 30 | 31 | 29 | 32 | 90 | 28 |
| 7 | 368 | 26 | 0.0715 | 9 | 9 | 10 | 5 | 8 | 5 | 11 | 11 | 29 | 28 | 29 | 30 | 29 | 31 | 88 | 27 |
| 8 | 366 | 26 | 0.0719 | 9 | 9 | 10 | 5 | 8 | 5 | 11 | 11 | 29 | 28 | 29 | 30 | 28 | 31 | 88 | 27 |
| 9 | 365 | 26 | 0.0742 | 9 | 9 | 10 | 5 | 8 | 5 | 11 | 11 | 29 | 28 | 29 | 30 | 28 | 31 | 87 | 27 |
| 10 | 372 | 26 | 0.0754 | 9 | 9 | 10 | 5 | 8 | 5 | 11 | 11 | 30 | 28 | 30 | 31 | 29 | 32 | 89 | 28 |

## Part 1 - sensitivity band

| state | scored | scored_positives | scored_calls | true_positive | false_positive |
|---|---|---|---|---|---|
| sep15 | 330 | 23.3020 | 26.1898 | 8.9361 | 17.2537 |
| sep15 | 350 | 24.7143 | 27.7770 | 9.4777 | 18.2993 |
| sep15 | 366 | 25.8441 | 29.0468 | 9.9110 | 19.1359 |
| sep15 | 374 | 26.4090 | 29.6817 | 10.1276 | 19.5541 |
| sep15 | 390 | 27.5388 | 30.9515 | 10.5609 | 20.3907 |
| sep15 | 410 | 28.9510 | 32.5388 | 11.1024 | 21.4364 |
| sep15 | 430 | 30.3633 | 34.1261 | 11.6440 | 22.4820 |
| oct03 | 330 | 23.2991 | 24.4054 | 9.8205 | 14.5849 |
| oct03 | 350 | 24.7111 | 25.8845 | 10.4157 | 15.4689 |
| oct03 | 366 | 25.8408 | 27.0678 | 10.8918 | 16.1760 |
| oct03 | 374 | 26.4056 | 27.6595 | 11.1299 | 16.5296 |
| oct03 | 390 | 27.5352 | 28.8428 | 11.6060 | 17.2368 |
| oct03 | 410 | 28.9473 | 30.3219 | 12.2012 | 18.1207 |
| oct03 | 430 | 30.3594 | 31.8010 | 12.7964 | 19.0046 |

## Part 2 - unskilled swap null

| transition | swapped_calls | observed_delta | observed_in_sd | null_mean | null_sd | null_p05 | null_p50 | null_p95 | null_fraction_at_least | null_fraction_as_extreme |
|---|---|---|---|---|---|---|---|---|---|---|
| sep15->oct03 | 17 | 0.0556 | 0.6354 | -0.1068 | 0.0574 | -0.2012 | -0.1045 | -0.0141 | 0.0018 | 0.8225 |
| sep14->sep15 | 9 | 0.0301 | 0.3478 | -0.0507 | 0.0443 | -0.1230 | -0.0458 | 0.0177 | 0.0322 | 0.7149 |
| sep15->sep18 | 27 | -0.1989 | -2.6987 | -0.1734 | 0.0637 | -0.2709 | -0.1744 | -0.0645 | 0.6489 | 0.3511 |

Bootstrap-mean versus point-estimate delta check: `{"draws": 2000, "point_mean": -0.10566234698921842, "point_sd": 0.05905497782754148, "bootstrap_mean": -0.10258637911243329, "bootstrap_sd": 0.05880091617543242, "mean_difference": 0.0030759678767851295}`

## Part 2 - exchangeable split of the changed calls

| transition | swapped_calls | observed_delta | scenarios | pooled_mean | pooled_sd | pooled_fraction_at_least | central_fraction_at_least | min_fraction_at_least | median_fraction_at_least | max_fraction_at_least |
|---|---|---|---|---|---|---|---|---|---|---|
| sep15->oct03 | 17 | 0.0556 | 110 | 0.0270 | 0.0435 | 0.1931 | 0.2381 | 0.0000 | 0.2238 | 0.2595 |
| sep14->sep15 | 9 | 0.0301 | 45 | 0.0159 | 0.0317 | 0.4552 | 0.4273 | 0.3763 | 0.4526 | 0.5019 |
| sep15->sep18 | 27 | -0.1989 | 116 | -0.0999 | 0.0512 | 0.9806 | 0.9973 | 0.9579 | 0.9807 | 1.0000 |

Central scenario per transition: the neutral composition in which the retained,
removed and added groups share one scored fraction and one precision among their
scored calls. Every consistent scenario is in `exchangeable_scenarios.csv`.

| transition | changed_positives | changed_negatives | changed_unscored | removed_true_positive | added_true_positive | removed_false_positive | added_false_positive | fraction_at_least |
|---|---|---|---|---|---|---|---|---|
| sep15->oct03 | 9 | 11 | 14 | 4 | 5 | 7 | 4 | 0.2381 |
| sep14->sep15 | 5 | 6 | 7 | 2 | 3 | 3 | 3 | 0.4273 |
| sep15->sep18 | 5 | 28 | 21 | 5 | 0 | 11 | 17 | 0.9973 |

## Part 3 - attribution between the two score columns

| group | n | legacy_rank | pool_rank | mean_rank | p_active | neighbour_positive_rate | in_pool_top49 |
|---|---|---|---|---|---|---|---|
| retained | 32 | 0.9681 | 0.9478 | 0.9579 | 0.3994 | 0.4661 | 22 |
| removed | 17 | 0.9659 | 0.6098 | 0.7879 | 0.1750 | 0.3483 | 0 |
| added | 17 | 0.9015 | 0.9530 | 0.9273 | 0.3094 | 0.3930 | 13 |

Overlaps: `{"pool_top49_overlap_sep15": 22, "pool_top49_overlap_oct03": 35, "pool_top49_overlap_retained": 22, "additions_in_pool_top49": 13, "removals_in_pool_top49": 0, "closer_vector": "oct03", "legacy_top49_overlap_oct03": 32}`

## Part 4 - the eighth concordance transition

| column | n_removed | n_added | mean_removed | mean_added | difference | permutation_p | sign_agreement | selection_component |
|---|---|---|---|---|---|---|---|---|
| gate_x_tdi.p_active | 17 | 17 | 0.1750 | 0.3094 | 0.1344 | 0.0017 | yes | no |
| gate_x_tdi_chemeleon.p_active | 17 | 17 | 0.1752 | 0.3281 | 0.1530 | 0.0055 | yes | no |
| pool_lgbm.score | 17 | 17 | 0.0653 | 0.4096 | 0.3443 | 0.0001 | yes | yes |
| pool_chemeleon.score | 17 | 17 | 0.5325 | 0.4448 | -0.0877 | 0.3089 | no | no |
| rank_mean_gates.score | 17 | 17 | 0.8356 | 0.8258 | -0.0099 | 0.8309 | no | no |
| gate_x_tdi.p_tdi | 17 | 17 | 0.4749 | 0.3833 | -0.0916 | 0.2476 | no | no |
| official_lgbm.score | 17 | 17 | 0.4749 | 0.3833 | -0.0916 | 0.2448 | no | no |

| column | transitions | sign_agreements | agreements_before | new_transition_agrees | significant | selection_component |
|---|---|---|---|---|---|---|
| gate_x_tdi.p_active | 8 | 7 | 6 | yes | 4 | no |
| gate_x_tdi_chemeleon.p_active | 8 | 7 | 6 | yes | 4 | no |
| pool_lgbm.score | 8 | 7 | 6 | yes | 4 | yes |
| pool_chemeleon.score | 8 | 6 | 6 | no | 2 | no |
| rank_mean_gates.score | 8 | 6 | 6 | no | 4 | no |
| gate_x_tdi.p_tdi | 8 | 5 | 5 | no | 4 | no |
| official_lgbm.score | 8 | 5 | 5 | no | 4 | no |

The deployed legacy score and `pool_lgbm` are the two components of the submitted
rank mean, so their agreement on this transition is partly tautological.
