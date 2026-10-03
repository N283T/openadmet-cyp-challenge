# Diagnosing the D6 activity gate on blind chemistry

## p_active distribution, blind against development strata

| population | n | mean | q10 | q50 | q90 | frac_ge_0.2 | frac_ge_0.3 | frac_ge_0.4 | frac_ge_0.5 |
|---|---|---|---|---|---|---|---|---|---|
| blind_750 | 750 | 0.3094 | 0.0916 | 0.2140 | 0.7205 | 0.5333 | 0.3720 | 0.2627 | 0.1827 |
| scaffold: development_all | 6146 | 0.2778 | 0.0702 | 0.1976 | 0.6431 | 0.4958 | 0.3349 | 0.2317 | 0.1653 |
| scaffold: screened_hit (log2fc < -0.5) | 1553 | 0.4401 | 0.1236 | 0.3895 | 0.8511 | 0.7721 | 0.6188 | 0.4836 | 0.3748 |
| scaffold: screened_non_hit (log2fc >= -0.5) | 2823 | 0.1826 | 0.0605 | 0.1314 | 0.3756 | 0.3036 | 0.1580 | 0.0850 | 0.0496 |
| scaffold: unscreened (no log2fc) | 1770 | 0.2873 | 0.0827 | 0.2265 | 0.6069 | 0.5599 | 0.3678 | 0.2446 | 0.1661 |
| scaffold: role official | 1497 | 0.4431 | 0.1246 | 0.3913 | 0.8544 | 0.7762 | 0.6246 | 0.4883 | 0.3801 |
| scaffold: role pseudo_negative | 2693 | 0.1768 | 0.0598 | 0.1292 | 0.3564 | 0.2896 | 0.1459 | 0.0765 | 0.0446 |
| scaffold: role borderline | 130 | 0.3012 | 0.0986 | 0.2451 | 0.5868 | 0.5923 | 0.4077 | 0.2615 | 0.1538 |
| scaffold: role unlabelled_hit | 58 | 0.3537 | 0.1075 | 0.2714 | 0.7033 | 0.6552 | 0.4483 | 0.3448 | 0.2241 |
| scaffold: role unscreened | 1768 | 0.2874 | 0.0827 | 0.2270 | 0.6070 | 0.5600 | 0.3682 | 0.2449 | 0.1663 |
| umap: development_all | 6146 | 0.2812 | 0.0820 | 0.2103 | 0.5960 | 0.5255 | 0.3529 | 0.2363 | 0.1578 |
| umap: screened_hit (log2fc < -0.5) | 1553 | 0.4172 | 0.1308 | 0.3729 | 0.7875 | 0.7785 | 0.6085 | 0.4623 | 0.3484 |
| umap: screened_non_hit (log2fc >= -0.5) | 2823 | 0.1958 | 0.0704 | 0.1483 | 0.3937 | 0.3429 | 0.1821 | 0.0967 | 0.0489 |
| umap: unscreened (no log2fc) | 1770 | 0.2982 | 0.0952 | 0.2427 | 0.5828 | 0.5949 | 0.4011 | 0.2605 | 0.1644 |
| umap: role official | 1497 | 0.4194 | 0.1327 | 0.3749 | 0.7884 | 0.7829 | 0.6126 | 0.4683 | 0.3534 |
| umap: role pseudo_negative | 2693 | 0.1902 | 0.0700 | 0.1450 | 0.3778 | 0.3301 | 0.1716 | 0.0884 | 0.0434 |
| umap: role borderline | 130 | 0.3124 | 0.1273 | 0.2611 | 0.5968 | 0.6077 | 0.4000 | 0.2692 | 0.1615 |
| umap: role unlabelled_hit | 58 | 0.3541 | 0.0991 | 0.2795 | 0.7402 | 0.6552 | 0.4828 | 0.2931 | 0.2069 |
| umap: role unscreened | 1768 | 0.2983 | 0.0952 | 0.2427 | 0.5832 | 0.5950 | 0.4016 | 0.2607 | 0.1646 |

## Distribution shift of the blind set against each reference

| reference | n_reference | mean_reference | mean_difference | ks_statistic | ks_p | auroc_blind_above_reference |
|---|---|---|---|---|---|---|
| scaffold: development_all | 6146 | 0.2778 | 0.0316 | 0.0949 | 0.0000 | 0.5469 |
| scaffold: screened_hit (log2fc < -0.5) | 1553 | 0.4401 | -0.1307 | 0.2517 | 0.0000 | 0.3431 |
| scaffold: screened_non_hit (log2fc >= -0.5) | 2823 | 0.1826 | 0.1269 | 0.2548 | 0.0000 | 0.6822 |
| scaffold: unscreened (no log2fc) | 1770 | 0.2873 | 0.0221 | 0.0590 | 0.0487 | 0.5098 |
| scaffold: role official | 1497 | 0.4431 | -0.1337 | 0.2573 | 0.0000 | 0.3404 |
| scaffold: role pseudo_negative | 2693 | 0.1768 | 0.1326 | 0.2692 | 0.0000 | 0.6920 |
| scaffold: role borderline | 130 | 0.3012 | 0.0082 | 0.0871 | 0.3480 | 0.4781 |
| scaffold: role unlabelled_hit | 58 | 0.3537 | -0.0442 | 0.1840 | 0.0454 | 0.4224 |
| scaffold: role unscreened | 1768 | 0.2874 | 0.0220 | 0.0590 | 0.0487 | 0.5097 |
| umap: development_all | 6146 | 0.2812 | 0.0282 | 0.0563 | 0.0275 | 0.5238 |
| umap: screened_hit (log2fc < -0.5) | 1553 | 0.4172 | -0.1077 | 0.2513 | 0.0000 | 0.3501 |
| umap: screened_non_hit (log2fc >= -0.5) | 2823 | 0.1958 | 0.1136 | 0.2037 | 0.0000 | 0.6459 |
| umap: unscreened (no log2fc) | 1770 | 0.2982 | 0.0112 | 0.0759 | 0.0043 | 0.4813 |
| umap: role official | 1497 | 0.4194 | -0.1099 | 0.2549 | 0.0000 | 0.3472 |
| umap: role pseudo_negative | 2693 | 0.1902 | 0.1192 | 0.2178 | 0.0000 | 0.6554 |
| umap: role borderline | 130 | 0.3124 | -0.0030 | 0.1440 | 0.0180 | 0.4497 |
| umap: role unlabelled_hit | 58 | 0.3541 | -0.0447 | 0.1593 | 0.1160 | 0.4302 |
| umap: role unscreened | 1768 | 0.2983 | 0.0112 | 0.0755 | 0.0046 | 0.4813 |

## Gate reliability on development held folds

| split_id | bin_low | bin_high | n | predicted | observed | reliability | observed_ci_low | observed_ci_high | predicted_inside_ci |
|---|---|---|---|---|---|---|---|---|---|
| scaffold | 0.0000 | 0.1000 | 1081 | 0.0714 | 0.0444 | -0.0270 | 0.0337 | 0.0584 | no |
| scaffold | 0.1000 | 0.2000 | 1239 | 0.1435 | 0.1340 | -0.0095 | 0.1161 | 0.1541 | yes |
| scaffold | 0.2000 | 0.3000 | 649 | 0.2435 | 0.2666 | 0.0231 | 0.2340 | 0.3019 | yes |
| scaffold | 0.3000 | 0.4000 | 416 | 0.3477 | 0.3822 | 0.0345 | 0.3368 | 0.4298 | yes |
| scaffold | 0.4000 | 0.5000 | 269 | 0.4496 | 0.5279 | 0.0783 | 0.4683 | 0.5867 | no |
| scaffold | 0.5000 | 0.6000 | 197 | 0.5484 | 0.5381 | -0.0103 | 0.4684 | 0.6063 | yes |
| scaffold | 0.6000 | 0.7000 | 148 | 0.6526 | 0.6351 | -0.0174 | 0.5551 | 0.7084 | yes |
| scaffold | 0.7000 | 0.8000 | 138 | 0.7501 | 0.7681 | 0.0180 | 0.6910 | 0.8307 | yes |
| scaffold | 0.8000 | 0.9000 | 130 | 0.8484 | 0.8154 | -0.0330 | 0.7400 | 0.8727 | yes |
| scaffold | 0.9000 | 1.0000 | 109 | 0.9342 | 0.8991 | -0.0351 | 0.8283 | 0.9427 | yes |
| umap | 0.0000 | 0.1000 | 866 | 0.0770 | 0.0520 | -0.0250 | 0.0391 | 0.0688 | no |
| umap | 0.1000 | 0.2000 | 1333 | 0.1447 | 0.1283 | -0.0165 | 0.1114 | 0.1473 | yes |
| umap | 0.2000 | 0.3000 | 718 | 0.2448 | 0.2660 | 0.0212 | 0.2350 | 0.2995 | yes |
| umap | 0.3000 | 0.4000 | 468 | 0.3465 | 0.3590 | 0.0125 | 0.3168 | 0.4034 | yes |
| umap | 0.4000 | 0.5000 | 312 | 0.4485 | 0.4712 | 0.0226 | 0.4165 | 0.5266 | yes |
| umap | 0.5000 | 0.6000 | 236 | 0.5490 | 0.5805 | 0.0315 | 0.5168 | 0.6417 | yes |
| umap | 0.6000 | 0.7000 | 166 | 0.6477 | 0.6807 | 0.0330 | 0.6064 | 0.7469 | yes |
| umap | 0.7000 | 0.8000 | 127 | 0.7474 | 0.7480 | 0.0006 | 0.6660 | 0.8155 | yes |
| umap | 0.8000 | 0.9000 | 100 | 0.8511 | 0.8600 | 0.0089 | 0.7786 | 0.9147 | yes |
| umap | 0.9000 | 1.0000 | 50 | 0.9235 | 0.9000 | -0.0235 | 0.7864 | 0.9565 | yes |

| split_id | n_screened | base_rate | mean_p_active | expected_calibration_error | max_bin_deviation | brier | brier_of_base_rate | auroc | bins_with_predicted_inside_ci | bins |
|---|---|---|---|---|---|---|---|---|---|---|
| scaffold | 4376 | 0.2738 | 0.2740 | 0.0243 | 0.0783 | 0.1418 | 0.1988 | 0.8252 | 8 | 10 |
| umap | 4376 | 0.2738 | 0.2744 | 0.0198 | 0.0330 | 0.1487 | 0.1988 | 0.8045 | 9 | 10 |

## Unscreened development rows

| split_id | n_unscreened | mean_unscreened | mean_screened | difference | auroc_unscreened_above_screened | mannwhitney_p | ks_blind_vs_unscreened | ks_blind_vs_screened |
|---|---|---|---|---|---|---|---|---|
| scaffold | 1770 | 0.2873 | 0.2740 | 0.0133 | 0.5510 | 0.0000 | 0.0590 | 0.1196 |
| umap | 1770 | 0.2982 | 0.2744 | 0.0238 | 0.5594 | 0.0000 | 0.0759 | 0.0697 |

## Officially labelled positive rate by p_active bin

| split_id | bin_low | bin_high | n | mean_p_active | positives | positive_rate | ci_low | ci_high |
|---|---|---|---|---|---|---|---|---|
| scaffold | 0.0000 | 0.2000 | 335 | 0.1297 | 77 | 0.2299 | 0.1880 | 0.2778 |
| scaffold | 0.2000 | 0.3000 | 227 | 0.2465 | 60 | 0.2643 | 0.2112 | 0.3253 |
| scaffold | 0.3000 | 0.5000 | 366 | 0.3933 | 84 | 0.2295 | 0.1894 | 0.2753 |
| scaffold | 0.5000 | 1.0000 | 569 | 0.7382 | 103 | 0.1810 | 0.1516 | 0.2148 |
| umap | 0.0000 | 0.2000 | 325 | 0.1346 | 75 | 0.2308 | 0.1883 | 0.2796 |
| umap | 0.2000 | 0.3000 | 255 | 0.2466 | 70 | 0.2745 | 0.2234 | 0.3323 |
| umap | 0.3000 | 0.5000 | 388 | 0.3936 | 82 | 0.2113 | 0.1736 | 0.2547 |
| umap | 0.5000 | 1.0000 | 529 | 0.6965 | 97 | 0.1834 | 0.1527 | 0.2186 |

| split_id | n_official | positives | base_rate | auroc_p_active_for_tdi | spearman | spearman_p | n_below_0_2 | rate_below_0_2 | n_above_0_5 | rate_above_0_5 |
|---|---|---|---|---|---|---|---|---|---|---|
| scaffold | 1497 | 324 | 0.2164 | 0.4608 | -0.0559 | 0.0304 | 335 | 0.2299 | 569 | 0.1810 |
| umap | 1497 | 324 | 0.2164 | 0.4514 | -0.0694 | 0.0072 | 325 | 0.2308 | 529 | 0.1834 |

## Primary-screen region decomposition of p_active

| split_id | bin_low | bin_high | n | share_strong_hit | share_weak_hit | share_borderline | share_relabelled | implied_positive_rate_clean | implied_positive_rate_pooled_upper_bound | implied_positive_rate_d6_at_blind_base_rate |
|---|---|---|---|---|---|---|---|---|---|---|
| scaffold | 0.0000 | 0.2000 | 2320 | 0.0922 | 0.0603 | 0.0228 | 0.8246 | 0.0334 | 0.0373 | 0.0917 |
| scaffold | 0.2000 | 0.3000 | 649 | 0.2666 | 0.1002 | 0.0370 | 0.5963 | 0.0799 | 0.0827 | 0.1220 |
| scaffold | 0.3000 | 0.5000 | 685 | 0.4394 | 0.1139 | 0.0482 | 0.3985 | 0.1202 | 0.1221 | 0.1484 |
| scaffold | 0.5000 | 1.0000 | 722 | 0.7064 | 0.0997 | 0.0277 | 0.1662 | 0.1743 | 0.1751 | 0.1860 |
| umap | 0.0000 | 0.2000 | 2199 | 0.0982 | 0.0582 | 0.0232 | 0.8204 | 0.0343 | 0.0381 | 0.0922 |
| umap | 0.2000 | 0.3000 | 718 | 0.2660 | 0.1017 | 0.0376 | 0.5947 | 0.0801 | 0.0829 | 0.1221 |
| umap | 0.3000 | 0.5000 | 780 | 0.4038 | 0.1141 | 0.0397 | 0.4423 | 0.1125 | 0.1146 | 0.1438 |
| umap | 0.5000 | 1.0000 | 679 | 0.7010 | 0.0957 | 0.0309 | 0.1723 | 0.1723 | 0.1731 | 0.1844 |

## Implied blind positive rate from the decomposition

| split_id | scenario | relabelled_rate | implied_blind_prevalence | ratio_to_published | implied_rate_blind_below_0_2 | implied_rate_blind_0_2_to_0_5 | implied_rate_blind_above_0_5 |
|---|---|---|---|---|---|---|---|
| scaffold | clean | 0.0000 | 0.0831 | 1.1770 | 0.0334 | 0.1017 | 0.1743 |
| scaffold | pooled_upper_bound | 0.0047 | 0.0859 | 1.2161 | 0.0373 | 0.1040 | 0.1751 |
| scaffold | d6_at_blind_base_rate | 0.0706 | 0.1245 | 1.7638 | 0.0917 | 0.1362 | 0.1860 |
| umap | clean | 0.0000 | 0.0817 | 1.1571 | 0.0343 | 0.0976 | 0.1723 |
| umap | pooled_upper_bound | 0.0047 | 0.0845 | 1.1967 | 0.0381 | 0.1000 | 0.1731 |
| umap | d6_at_blind_base_rate | 0.0706 | 0.1236 | 1.7511 | 0.0922 | 0.1338 | 0.1844 |

## Circularity of a pseudo-negative-labelled metric for a gate filter

| split_id | n_official | n_pseudo_negative | auroc_p_active_official_vs_pseudo | mean_official | mean_pseudo_negative | official_share_above_0.2 | official_share_above_0.3 | official_share_above_0.4 | official_share_above_0.5 |
|---|---|---|---|---|---|---|---|---|---|
| scaffold | 1497 | 2693 | 0.8233 | 0.4431 | 0.1768 | 0.5984 | 0.7041 | 0.7801 | 0.8258 |
| umap | 1497 | 2693 | 0.8047 | 0.4194 | 0.1902 | 0.5687 | 0.6650 | 0.7465 | 0.8189 |

