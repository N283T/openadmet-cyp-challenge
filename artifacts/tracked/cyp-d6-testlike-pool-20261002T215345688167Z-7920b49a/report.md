# D6 test-like pool evaluation

Development molecules evaluated: 6146; smoke: False.

## Pool counts

```json
{
  "counts": {
    "blind": 750,
    "borderline": 130,
    "official": 1497,
    "pseudo_negative": 2693,
    "unlabelled_hit": 58,
    "unscreened": 1768
  },
  "development": 6146,
  "official_positive": 324,
  "official_with_log2fc": 1495,
  "official_log2fc_max": -0.5916683612,
  "official_positive_log2fc_max": -0.5916683612,
  "official_positive_log2fc_ge_gate": 66,
  "testlike_n": 4190,
  "testlike_positive_rate": 0.07732696897374702,
  "screened_development": 4376,
  "gate_active_development": 1198,
  "primary_source_rows": 17504
}
```

## Ranking (pooled OOF)

| split_id | model | auroc_official | ap_official | auroc_testlike | ap_testlike |
|---|---|---|---|---|---|
| umap | gate_x_tdi | 0.516 | 0.228 | 0.601 | 0.116 |
| umap | gate_x_tdi_chemeleon | 0.532 | 0.239 | 0.619 | 0.120 |
| umap | legacy_sep15_refit | 0.586 | 0.301 | 0.542 | 0.101 |
| umap | official_lgbm | 0.534 | 0.251 | 0.509 | 0.085 |
| umap | pool_chemeleon | 0.527 | 0.230 | 0.611 | 0.113 |
| umap | pool_lgbm | 0.528 | 0.239 | 0.639 | 0.133 |
| umap | rank_mean_gates | 0.534 | 0.239 | 0.632 | 0.126 |
| scaffold | gate_x_tdi | 0.566 | 0.262 | 0.644 | 0.144 |
| scaffold | gate_x_tdi_chemeleon | 0.572 | 0.259 | 0.653 | 0.139 |
| scaffold | legacy_sep15_refit | 0.638 | 0.348 | 0.587 | 0.123 |
| scaffold | official_lgbm | 0.584 | 0.280 | 0.536 | 0.095 |
| scaffold | pool_chemeleon | 0.569 | 0.259 | 0.648 | 0.132 |
| scaffold | pool_lgbm | 0.579 | 0.282 | 0.685 | 0.167 |
| scaffold | rank_mean_gates | 0.583 | 0.277 | 0.672 | 0.159 |

## Decision (pooled held-fold calls)

| split_id | model | nested_mcc_testlike | nested_mcc_official | nested_calls | fixed653_mcc_testlike | fixed653_mcc_official | fixed653_calls_borderline | fixed653_calls_unscreened |
|---|---|---|---|---|---|---|---|---|
| umap | gate_x_tdi | 0.078 | 0.012 | 579 | 0.081 | 0.026 | 7 | 137 |
| umap | gate_x_tdi_chemeleon | 0.068 | 0.014 | 554 | 0.076 | 0.031 | 11 | 132 |
| umap | legacy_sep15_refit | 0.068 | 0.128 | 393 | 0.074 | 0.147 | 3 | 143 |
| umap | official_lgbm | 0.042 | 0.093 | 616 | 0.034 | 0.062 | 7 | 145 |
| umap | pool_chemeleon | 0.076 | 0.025 | 491 | 0.050 | -0.004 | 11 | 121 |
| umap | pool_lgbm | 0.102 | 0.028 | 616 | 0.096 | 0.031 | 6 | 129 |
| umap | rank_mean_gates | 0.085 | 0.017 | 517 | 0.089 | 0.033 | 9 | 145 |
| scaffold | gate_x_tdi | 0.112 | 0.066 | 578 | 0.119 | 0.071 | 7 | 138 |
| scaffold | gate_x_tdi_chemeleon | 0.123 | 0.075 | 615 | 0.113 | 0.062 | 9 | 99 |
| scaffold | legacy_sep15_refit | 0.075 | 0.136 | 442 | 0.087 | 0.145 | 7 | 136 |
| scaffold | official_lgbm | 0.053 | 0.099 | 257 | 0.043 | 0.081 | 5 | 137 |
| scaffold | pool_chemeleon | 0.092 | 0.053 | 468 | 0.094 | 0.051 | 9 | 115 |
| scaffold | pool_lgbm | 0.155 | 0.110 | 454 | 0.155 | 0.115 | 9 | 107 |
| scaffold | rank_mean_gates | 0.123 | 0.073 | 590 | 0.121 | 0.074 | 5 | 121 |

## Test-like MCC by call rate (pooled)

| split_id | model | 3.00% | 4.00% | 5.00% | 6.00% | 7.00% | 8.00% | 10.00% | 6.53% |
|---|---|---|---|---|---|---|---|---|---|
| umap | gate_x_tdi | 0.043 | 0.064 | 0.079 | 0.077 | 0.085 | 0.078 | 0.089 | 0.081 |
| umap | gate_x_tdi_chemeleon | 0.041 | 0.068 | 0.076 | 0.080 | 0.080 | 0.086 | 0.091 | 0.076 |
| umap | legacy_sep15_refit | 0.056 | 0.073 | 0.068 | 0.068 | 0.070 | 0.081 | 0.071 | 0.074 |
| umap | official_lgbm | 0.017 | 0.025 | 0.027 | 0.027 | 0.032 | 0.032 | 0.042 | 0.034 |
| umap | pool_chemeleon | 0.030 | 0.036 | 0.047 | 0.052 | 0.065 | 0.076 | 0.071 | 0.050 |
| umap | pool_lgbm | 0.072 | 0.076 | 0.087 | 0.091 | 0.091 | 0.090 | 0.102 | 0.096 |
| umap | rank_mean_gates | 0.064 | 0.080 | 0.086 | 0.091 | 0.098 | 0.096 | 0.097 | 0.089 |
| scaffold | gate_x_tdi | 0.090 | 0.107 | 0.107 | 0.106 | 0.114 | 0.118 | 0.123 | 0.119 |
| scaffold | gate_x_tdi_chemeleon | 0.081 | 0.099 | 0.116 | 0.114 | 0.107 | 0.099 | 0.123 | 0.113 |
| scaffold | legacy_sep15_refit | 0.075 | 0.084 | 0.093 | 0.094 | 0.088 | 0.083 | 0.090 | 0.087 |
| scaffold | official_lgbm | 0.061 | 0.059 | 0.065 | 0.050 | 0.044 | 0.040 | 0.049 | 0.043 |
| scaffold | pool_chemeleon | 0.069 | 0.077 | 0.083 | 0.102 | 0.093 | 0.099 | 0.100 | 0.094 |
| scaffold | pool_lgbm | 0.112 | 0.124 | 0.125 | 0.154 | 0.162 | 0.162 | 0.156 | 0.155 |
| scaffold | rank_mean_gates | 0.117 | 0.098 | 0.113 | 0.109 | 0.124 | 0.128 | 0.132 | 0.121 |

## Paired group bootstrap versus official_lgbm (test-like labels)

| split_id | model | operating_point | delta | ci_low | ci_high | p_delta_gt_0 |
|---|---|---|---|---|---|---|
| umap | gate_x_tdi | nested | 0.036 | 0.010 | 0.066 | 0.994 |
| umap | gate_x_tdi | fixed_6.53pct | 0.047 | 0.014 | 0.083 | 0.998 |
| umap | gate_x_tdi_chemeleon | nested | 0.026 | -0.027 | 0.083 | 0.821 |
| umap | gate_x_tdi_chemeleon | fixed_6.53pct | 0.042 | -0.019 | 0.111 | 0.923 |
| umap | legacy_sep15_refit | nested | 0.026 | -0.018 | 0.074 | 0.873 |
| umap | legacy_sep15_refit | fixed_6.53pct | 0.040 | -0.008 | 0.091 | 0.943 |
| umap | pool_chemeleon | nested | 0.034 | 0.001 | 0.069 | 0.977 |
| umap | pool_chemeleon | fixed_6.53pct | 0.015 | -0.022 | 0.060 | 0.783 |
| umap | pool_lgbm | nested | 0.060 | 0.024 | 0.099 | 0.999 |
| umap | pool_lgbm | fixed_6.53pct | 0.061 | 0.018 | 0.109 | 0.997 |
| umap | rank_mean_gates | nested | 0.043 | 0.011 | 0.079 | 0.995 |
| umap | rank_mean_gates | fixed_6.53pct | 0.054 | 0.010 | 0.108 | 0.992 |
| scaffold | gate_x_tdi | nested | 0.059 | 0.018 | 0.101 | 0.999 |
| scaffold | gate_x_tdi | fixed_6.53pct | 0.076 | 0.039 | 0.121 | 1.000 |
| scaffold | gate_x_tdi_chemeleon | nested | 0.070 | 0.020 | 0.119 | 0.993 |
| scaffold | gate_x_tdi_chemeleon | fixed_6.53pct | 0.070 | 0.020 | 0.119 | 0.997 |
| scaffold | legacy_sep15_refit | nested | 0.022 | -0.023 | 0.065 | 0.828 |
| scaffold | legacy_sep15_refit | fixed_6.53pct | 0.045 | 0.002 | 0.085 | 0.980 |
| scaffold | pool_chemeleon | nested | 0.039 | -0.009 | 0.091 | 0.939 |
| scaffold | pool_chemeleon | fixed_6.53pct | 0.051 | 0.001 | 0.104 | 0.977 |
| scaffold | pool_lgbm | nested | 0.102 | 0.054 | 0.156 | 1.000 |
| scaffold | pool_lgbm | fixed_6.53pct | 0.113 | 0.068 | 0.163 | 1.000 |
| scaffold | rank_mean_gates | nested | 0.070 | 0.026 | 0.117 | 0.999 |
| scaffold | rank_mean_gates | fixed_6.53pct | 0.078 | 0.035 | 0.127 | 0.999 |

## Blind diagnostics

| model | sep15_positives | top49_overlap_sep15 | p_active_gt_0_5 | mean_p_active | umap_rate | umap_calls | umap_overlap_sep15 | umap_oof_threshold | umap_threshold_calls | umap_nested_rates | scaffold_rate | scaffold_calls | scaffold_overlap_sep15 | scaffold_oof_threshold | scaffold_threshold_calls | scaffold_nested_rates |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gate_x_tdi | 49 | 27 | 137 | 0.309 | 0.100 | 75 | 33 | 0.106 | 52 | 0.10,0.10,0.10,0.07,0.10 | 0.100 | 75 | 33 | 0.109 | 49 | 0.10,0.07,0.10,0.10,0.10 |
| gate_x_tdi_chemeleon | 49 | 24 | 144 | 0.301 | 0.100 | 75 | 30 | 0.043 | 119 | 0.10,0.08,0.10,0.10,0.07 | 0.100 | 75 | 30 | 0.048 | 106 | 0.10,0.10,0.10,0.10,0.10 |
| official_lgbm | 49 | 31 | – | – | 0.100 | 75 | 36 | 0.413 | 70 | 0.10,0.10,0.10,0.10,0.10 | 0.050 | 38 | 25 | 0.573 | 36 | 0.03,0.05,0.05,0.03,0.05 |
| pool_chemeleon | 49 | 17 | – | – | 0.080 | 60 | 20 | 0.660 | 79 | 0.08,0.08,0.08,0.08,0.08 | 0.060 | 45 | 16 | 0.733 | 59 | 0.06,0.06,0.10,0.10,0.06 |
| pool_lgbm | 49 | 22 | – | – | 0.100 | 75 | 25 | 0.172 | 89 | 0.10,0.10,0.10,0.10,0.10 | 0.080 | 60 | 23 | 0.204 | 76 | 0.08,0.07,0.07,0.08,0.07 |
| rank_mean_gates | 49 | 27 | – | – | 0.070 | 53 | 30 | 0.886 | 56 | 0.07,0.10,0.08,0.10,0.07 | 0.100 | 75 | 34 | 0.848 | 74 | 0.08,0.10,0.10,0.10,0.10 |

Legacy comparator audit: `{"used": true, "missing": 0, "fold_mismatch": 0, "seed_counts": [2], "fold_inconsistent_across_seeds": 0, "label_mismatch": 0, "rows": 12292}`
