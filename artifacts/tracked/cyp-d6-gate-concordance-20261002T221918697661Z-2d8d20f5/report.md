# D6 activity-gate historical concordance and legacy x gate OOF variants

Smoke: False. Rate grid: 3.00%, 4.00%, 5.00%, 6.53%, 8.00%, 10.00%.

## Part 1 - submitted payloads

| payload | positives | live_mcc | mean_gate_x_tdi.p_active | mean_gate_x_tdi_chemeleon.p_active | n_below_0.3_gate_x_tdi.p_active |
|---|---|---|---|---|---|
| sep13 | 49 | 0.271 | 0.301 | 0.277 | 28 |
| sep14 | 49 | 0.279 | 0.295 | 0.274 | 29 |
| sep15 | 49 | 0.309 | 0.322 | 0.304 | 27 |
| sep18 | 49 | 0.110 | 0.283 | 0.254 | 29 |
| sep23 | 49 | 0.235 | 0.290 | 0.271 | 32 |
| oct01 | 49 | 0.105 | 0.236 | 0.219 | 37 |
| oct02 | 175 | 0.118 | 0.215 | 0.216 | 139 |

## Part 1 - transition concordance

| transition | lb_delta | column | n_removed | n_added | mean_removed | mean_added | difference | permutation_p | sign_agreement |
|---|---|---|---|---|---|---|---|---|---|
| sep13->sep14 | 0.008 | gate_x_tdi.p_active | 6 | 6 | 0.278 | 0.223 | -0.055 | 0.558 | no |
| sep13->sep14 | 0.008 | gate_x_tdi_chemeleon.p_active | 6 | 6 | 0.250 | 0.218 | -0.032 | 0.693 | no |
| sep13->sep14 | 0.008 | pool_lgbm.score | 6 | 6 | 0.083 | 0.153 | 0.070 | 0.166 | yes |
| sep13->sep14 | 0.008 | pool_chemeleon.score | 6 | 6 | 0.517 | 0.283 | -0.233 | 0.174 | no |
| sep13->sep14 | 0.008 | rank_mean_gates.score | 6 | 6 | 0.886 | 0.793 | -0.094 | 0.216 | no |
| sep13->sep14 | 0.008 | gate_x_tdi.p_tdi | 6 | 6 | 0.413 | 0.370 | -0.044 | 0.660 | no |
| sep13->sep14 | 0.008 | official_lgbm.score | 6 | 6 | 0.413 | 0.370 | -0.044 | 0.660 | no |
| sep14->sep15 | 0.030 | gate_x_tdi.p_active | 9 | 9 | 0.164 | 0.310 | 0.146 | 0.041 | yes |
| sep14->sep15 | 0.030 | gate_x_tdi_chemeleon.p_active | 9 | 9 | 0.161 | 0.325 | 0.164 | 0.010 | yes |
| sep14->sep15 | 0.030 | pool_lgbm.score | 9 | 9 | 0.098 | 0.271 | 0.173 | 0.058 | yes |
| sep14->sep15 | 0.030 | pool_chemeleon.score | 9 | 9 | 0.339 | 0.554 | 0.215 | 0.183 | yes |
| sep14->sep15 | 0.030 | rank_mean_gates.score | 9 | 9 | 0.656 | 0.904 | 0.247 | 0.004 | yes |
| sep14->sep15 | 0.030 | gate_x_tdi.p_tdi | 9 | 9 | 0.236 | 0.464 | 0.228 | 0.040 | yes |
| sep14->sep15 | 0.030 | official_lgbm.score | 9 | 9 | 0.236 | 0.464 | 0.228 | 0.039 | yes |
| sep15->sep18 | -0.199 | gate_x_tdi.p_active | 27 | 27 | 0.315 | 0.245 | -0.070 | 0.151 | yes |
| sep15->sep18 | -0.199 | gate_x_tdi_chemeleon.p_active | 27 | 27 | 0.312 | 0.221 | -0.091 | 0.056 | yes |
| sep15->sep18 | -0.199 | pool_lgbm.score | 27 | 27 | 0.328 | 0.125 | -0.203 | 0.003 | yes |
| sep15->sep18 | -0.199 | pool_chemeleon.score | 27 | 27 | 0.498 | 0.375 | -0.123 | 0.153 | yes |
| sep15->sep18 | -0.199 | rank_mean_gates.score | 27 | 27 | 0.835 | 0.729 | -0.105 | 0.025 | yes |
| sep15->sep18 | -0.199 | gate_x_tdi.p_tdi | 27 | 27 | 0.495 | 0.301 | -0.193 | 0.003 | yes |
| sep15->sep18 | -0.199 | official_lgbm.score | 27 | 27 | 0.495 | 0.301 | -0.193 | 0.003 | yes |
| sep15->sep23 | -0.074 | gate_x_tdi.p_active | 15 | 15 | 0.294 | 0.191 | -0.103 | 0.104 | yes |
| sep15->sep23 | -0.074 | gate_x_tdi_chemeleon.p_active | 15 | 15 | 0.311 | 0.206 | -0.105 | 0.124 | yes |
| sep15->sep23 | -0.074 | pool_lgbm.score | 15 | 15 | 0.209 | 0.268 | 0.059 | 0.545 | no |
| sep15->sep23 | -0.074 | pool_chemeleon.score | 15 | 15 | 0.539 | 0.450 | -0.089 | 0.333 | yes |
| sep15->sep23 | -0.074 | rank_mean_gates.score | 15 | 15 | 0.806 | 0.721 | -0.085 | 0.221 | yes |
| sep15->sep23 | -0.074 | gate_x_tdi.p_tdi | 15 | 15 | 0.367 | 0.377 | 0.011 | 0.902 | no |
| sep15->sep23 | -0.074 | official_lgbm.score | 15 | 15 | 0.367 | 0.377 | 0.011 | 0.907 | no |
| sep18->sep23 | 0.125 | gate_x_tdi.p_active | 21 | 21 | 0.254 | 0.271 | 0.017 | 0.733 | yes |
| sep18->sep23 | 0.125 | gate_x_tdi_chemeleon.p_active | 21 | 21 | 0.231 | 0.272 | 0.042 | 0.373 | yes |
| sep18->sep23 | 0.125 | pool_lgbm.score | 21 | 21 | 0.097 | 0.400 | 0.303 | 0.000 | yes |
| sep18->sep23 | 0.125 | pool_chemeleon.score | 21 | 21 | 0.358 | 0.453 | 0.095 | 0.323 | yes |
| sep18->sep23 | 0.125 | rank_mean_gates.score | 21 | 21 | 0.728 | 0.802 | 0.075 | 0.186 | yes |
| sep18->sep23 | 0.125 | gate_x_tdi.p_tdi | 21 | 21 | 0.285 | 0.541 | 0.257 | 0.001 | yes |
| sep18->sep23 | 0.125 | official_lgbm.score | 21 | 21 | 0.285 | 0.541 | 0.257 | 0.001 | yes |
| sep15->oct01 | -0.204 | gate_x_tdi.p_active | 21 | 21 | 0.365 | 0.164 | -0.200 | 0.000 | yes |
| sep15->oct01 | -0.204 | gate_x_tdi_chemeleon.p_active | 21 | 21 | 0.365 | 0.168 | -0.197 | 0.000 | yes |
| sep15->oct01 | -0.204 | pool_lgbm.score | 21 | 21 | 0.364 | 0.232 | -0.133 | 0.121 | yes |
| sep15->oct01 | -0.204 | pool_chemeleon.score | 21 | 21 | 0.523 | 0.313 | -0.210 | 0.022 | yes |
| sep15->oct01 | -0.204 | rank_mean_gates.score | 21 | 21 | 0.877 | 0.677 | -0.200 | 0.000 | yes |
| sep15->oct01 | -0.204 | gate_x_tdi.p_tdi | 21 | 21 | 0.493 | 0.367 | -0.126 | 0.136 | yes |
| sep15->oct01 | -0.204 | official_lgbm.score | 21 | 21 | 0.493 | 0.367 | -0.126 | 0.135 | yes |
| sep15->oct02 | -0.191 | gate_x_tdi.p_active | 6 | 132 | 0.435 | 0.185 | -0.249 | 0.000 | yes |
| sep15->oct02 | -0.191 | gate_x_tdi_chemeleon.p_active | 6 | 132 | 0.435 | 0.193 | -0.241 | 0.000 | yes |
| sep15->oct02 | -0.191 | pool_lgbm.score | 6 | 132 | 0.340 | 0.106 | -0.234 | 0.006 | yes |
| sep15->oct02 | -0.191 | pool_chemeleon.score | 6 | 132 | 0.641 | 0.306 | -0.335 | 0.004 | yes |
| sep15->oct02 | -0.191 | rank_mean_gates.score | 6 | 132 | 0.945 | 0.606 | -0.339 | 0.000 | yes |
| sep15->oct02 | -0.191 | gate_x_tdi.p_tdi | 6 | 132 | 0.496 | 0.234 | -0.262 | 0.004 | yes |
| sep15->oct02 | -0.191 | official_lgbm.score | 6 | 132 | 0.496 | 0.234 | -0.262 | 0.005 | yes |

## Part 1 - sign agreement by score column

| column | agreeing_transitions | transitions | min_permutation_p | significant_transitions |
|---|---|---|---|---|
| gate_x_tdi.p_active | 6 | 7 | 0.000 | 3 |
| gate_x_tdi_chemeleon.p_active | 6 | 7 | 0.000 | 3 |
| pool_lgbm.score | 6 | 7 | 0.000 | 3 |
| pool_chemeleon.score | 6 | 7 | 0.004 | 2 |
| rank_mean_gates.score | 6 | 7 | 0.000 | 4 |
| official_lgbm.score | 5 | 7 | 0.001 | 4 |
| gate_x_tdi.p_tdi | 5 | 7 | 0.001 | 4 |

## Part 2 - OOF pools

| split_id | molecules | official | testlike | testlike_strict | legacy_gate_spearman | gate_ge_0_2 | gate_ge_0_3 | gate_ge_0_4 | legacy_rank_mean |
|---|---|---|---|---|---|---|---|---|---|
| umap | 6146 | 1497 | 4190 | 4017 | -0.054 | 3230 | 2169 | 1452 | 0.500 |
| scaffold | 6146 | 1497 | 4190 | 4017 | -0.062 | 3047 | 2058 | 1424 | 0.500 |

## Part 2 - nested selection (pooled held folds)

| split_id | label_variant | model | n_calls | mcc_testlike | precision_testlike | tp_testlike | fp_testlike | mcc_official | precision_official | tp_official | fp_official |
|---|---|---|---|---|---|---|---|---|---|---|---|
| umap | testlike | L0_legacy | 393 | 0.068 | 0.152 | 36 | 201 | 0.128 | 0.434 | 36 | 47 |
| umap | testlike | L1_legacy_x_pa_morgan | 616 | 0.106 | 0.168 | 63 | 312 | 0.015 | 0.229 | 63 | 212 |
| umap | testlike | L2_legacy_x_sqrt_pa_morgan | 616 | 0.128 | 0.187 | 69 | 300 | 0.062 | 0.273 | 69 | 184 |
| umap | testlike | L3_legacy_x_sqrt_pa_both | 616 | 0.114 | 0.175 | 65 | 306 | 0.013 | 0.227 | 65 | 221 |
| umap | testlike | L4_rankmean_legacy_pool | 591 | 0.107 | 0.170 | 62 | 302 | 0.106 | 0.332 | 62 | 125 |
| umap | testlike | L5_gate_filter_tau0.2 | 554 | 0.111 | 0.180 | 58 | 264 | 0.098 | 0.326 | 58 | 120 |
| umap | testlike | L5_gate_filter_tau0.3 | 616 | 0.105 | 0.169 | 61 | 300 | 0.057 | 0.272 | 61 | 163 |
| umap | testlike | L5_gate_filter_tau0.4 | 566 | 0.087 | 0.158 | 51 | 272 | 0.010 | 0.227 | 51 | 174 |
| umap | testlike | B_pool_lgbm | 616 | 0.102 | 0.162 | 64 | 331 | 0.028 | 0.242 | 64 | 201 |
| umap | testlike_strict | L0_legacy | 393 | 0.071 | 0.159 | 36 | 190 | 0.128 | 0.434 | 36 | 47 |
| umap | testlike_strict | L1_legacy_x_pa_morgan | 616 | 0.108 | 0.175 | 63 | 298 | 0.015 | 0.229 | 63 | 212 |
| umap | testlike_strict | L2_legacy_x_sqrt_pa_morgan | 616 | 0.130 | 0.194 | 69 | 287 | 0.062 | 0.273 | 69 | 184 |
| umap | testlike_strict | L3_legacy_x_sqrt_pa_both | 616 | 0.116 | 0.182 | 65 | 293 | 0.013 | 0.227 | 65 | 221 |
| umap | testlike_strict | L4_rankmean_legacy_pool | 591 | 0.108 | 0.176 | 62 | 291 | 0.106 | 0.332 | 62 | 125 |
| umap | testlike_strict | L5_gate_filter_tau0.2 | 554 | 0.113 | 0.186 | 58 | 253 | 0.098 | 0.326 | 58 | 120 |
| umap | testlike_strict | L5_gate_filter_tau0.3 | 616 | 0.107 | 0.175 | 61 | 288 | 0.057 | 0.272 | 61 | 163 |
| umap | testlike_strict | L5_gate_filter_tau0.4 | 566 | 0.089 | 0.165 | 51 | 258 | 0.010 | 0.227 | 51 | 174 |
| umap | testlike_strict | B_pool_lgbm | 616 | 0.103 | 0.167 | 64 | 319 | 0.028 | 0.242 | 64 | 201 |
| scaffold | testlike | L0_legacy | 491 | 0.072 | 0.146 | 44 | 257 | 0.131 | 0.411 | 44 | 63 |
| scaffold | testlike | L1_legacy_x_pa_morgan | 590 | 0.159 | 0.215 | 78 | 284 | 0.070 | 0.277 | 78 | 204 |
| scaffold | testlike | L2_legacy_x_sqrt_pa_morgan | 615 | 0.170 | 0.223 | 83 | 289 | 0.112 | 0.317 | 83 | 179 |
| scaffold | testlike | L3_legacy_x_sqrt_pa_both | 615 | 0.176 | 0.226 | 86 | 295 | 0.084 | 0.286 | 86 | 215 |
| scaffold | testlike | L4_rankmean_legacy_pool | 424 | 0.116 | 0.200 | 50 | 200 | 0.124 | 0.382 | 50 | 81 |
| scaffold | testlike | L5_gate_filter_tau0.2 | 572 | 0.161 | 0.223 | 75 | 262 | 0.172 | 0.405 | 75 | 110 |
| scaffold | testlike | L5_gate_filter_tau0.3 | 362 | 0.153 | 0.256 | 53 | 154 | 0.140 | 0.402 | 53 | 79 |
| scaffold | testlike | L5_gate_filter_tau0.4 | 454 | 0.128 | 0.208 | 56 | 213 | 0.057 | 0.276 | 56 | 147 |
| scaffold | testlike | B_pool_lgbm | 472 | 0.154 | 0.219 | 71 | 253 | 0.110 | 0.326 | 71 | 147 |
| scaffold | testlike_strict | L0_legacy | 491 | 0.072 | 0.151 | 44 | 247 | 0.131 | 0.411 | 44 | 63 |
| scaffold | testlike_strict | L1_legacy_x_pa_morgan | 590 | 0.161 | 0.223 | 78 | 272 | 0.070 | 0.277 | 78 | 204 |
| scaffold | testlike_strict | L2_legacy_x_sqrt_pa_morgan | 615 | 0.173 | 0.231 | 83 | 277 | 0.112 | 0.317 | 83 | 179 |
| scaffold | testlike_strict | L3_legacy_x_sqrt_pa_both | 615 | 0.179 | 0.234 | 86 | 281 | 0.084 | 0.286 | 86 | 215 |
| scaffold | testlike_strict | L4_rankmean_legacy_pool | 424 | 0.119 | 0.210 | 50 | 188 | 0.124 | 0.382 | 50 | 81 |
| scaffold | testlike_strict | L5_gate_filter_tau0.2 | 572 | 0.165 | 0.232 | 75 | 248 | 0.172 | 0.405 | 75 | 110 |
| scaffold | testlike_strict | L5_gate_filter_tau0.3 | 362 | 0.156 | 0.266 | 53 | 146 | 0.140 | 0.402 | 53 | 79 |
| scaffold | testlike_strict | L5_gate_filter_tau0.4 | 478 | 0.145 | 0.228 | 61 | 206 | 0.071 | 0.289 | 61 | 150 |
| scaffold | testlike_strict | B_pool_lgbm | 472 | 0.155 | 0.225 | 71 | 244 | 0.110 | 0.326 | 71 | 147 |

## Part 2 - paired group bootstrap versus the legacy control

| split_id | model | operating_point | label_set | delta | ci_low | ci_high | p_delta_gt_0 |
|---|---|---|---|---|---|---|---|
| umap | L1_legacy_x_pa_morgan | nested | testlike | 0.038 | -0.002 | 0.080 | 0.969 |
| umap | L1_legacy_x_pa_morgan | nested | official | -0.113 | -0.161 | -0.058 | 0.000 |
| umap | L1_legacy_x_pa_morgan | historical | testlike | 0.015 | -0.032 | 0.066 | 0.734 |
| umap | L1_legacy_x_pa_morgan | historical | official | -0.142 | -0.206 | -0.070 | 0.001 |
| umap | L2_legacy_x_sqrt_pa_morgan | nested | testlike | 0.059 | 0.016 | 0.105 | 0.997 |
| umap | L2_legacy_x_sqrt_pa_morgan | nested | official | -0.066 | -0.118 | -0.007 | 0.013 |
| umap | L2_legacy_x_sqrt_pa_morgan | historical | testlike | 0.027 | -0.015 | 0.072 | 0.889 |
| umap | L2_legacy_x_sqrt_pa_morgan | historical | official | -0.107 | -0.165 | -0.045 | 0.001 |
| umap | L3_legacy_x_sqrt_pa_both | nested | testlike | 0.046 | -0.010 | 0.102 | 0.954 |
| umap | L3_legacy_x_sqrt_pa_both | nested | official | -0.115 | -0.180 | -0.041 | 0.002 |
| umap | L3_legacy_x_sqrt_pa_both | historical | testlike | 0.017 | -0.030 | 0.067 | 0.771 |
| umap | L3_legacy_x_sqrt_pa_both | historical | official | -0.150 | -0.209 | -0.084 | 0.001 |
| umap | L4_rankmean_legacy_pool | nested | testlike | 0.039 | 0.006 | 0.072 | 0.993 |
| umap | L4_rankmean_legacy_pool | nested | official | -0.022 | -0.076 | 0.040 | 0.234 |
| umap | L4_rankmean_legacy_pool | historical | testlike | 0.029 | -0.014 | 0.073 | 0.916 |
| umap | L4_rankmean_legacy_pool | historical | official | -0.041 | -0.102 | 0.031 | 0.115 |
| umap | L5_gate_filter_tau0.2 | nested | testlike | 0.043 | 0.010 | 0.078 | 0.994 |
| umap | L5_gate_filter_tau0.2 | nested | official | -0.030 | -0.074 | 0.018 | 0.107 |
| umap | L5_gate_filter_tau0.2 | historical | testlike | 0.020 | -0.013 | 0.055 | 0.881 |
| umap | L5_gate_filter_tau0.2 | historical | official | -0.051 | -0.097 | -0.001 | 0.024 |
| umap | L5_gate_filter_tau0.3 | nested | testlike | 0.037 | -0.003 | 0.077 | 0.970 |
| umap | L5_gate_filter_tau0.3 | nested | official | -0.071 | -0.125 | -0.016 | 0.006 |
| umap | L5_gate_filter_tau0.3 | historical | testlike | 0.011 | -0.031 | 0.054 | 0.699 |
| umap | L5_gate_filter_tau0.3 | historical | official | -0.100 | -0.159 | -0.037 | 0.001 |
| umap | L5_gate_filter_tau0.4 | nested | testlike | 0.019 | -0.026 | 0.067 | 0.809 |
| umap | L5_gate_filter_tau0.4 | nested | official | -0.117 | -0.171 | -0.053 | 0.001 |
| umap | L5_gate_filter_tau0.4 | historical | testlike | 0.021 | -0.016 | 0.065 | 0.852 |
| umap | L5_gate_filter_tau0.4 | historical | official | -0.114 | -0.166 | -0.055 | 0.001 |
| umap | B_pool_lgbm | nested | testlike | 0.034 | -0.015 | 0.086 | 0.905 |
| umap | B_pool_lgbm | nested | official | -0.100 | -0.160 | -0.027 | 0.005 |
| umap | B_pool_lgbm | historical | testlike | 0.022 | -0.028 | 0.077 | 0.804 |
| umap | B_pool_lgbm | historical | official | -0.116 | -0.182 | -0.043 | 0.001 |
| scaffold | L1_legacy_x_pa_morgan | nested | testlike | 0.087 | 0.036 | 0.140 | 1.000 |
| scaffold | L1_legacy_x_pa_morgan | nested | official | -0.061 | -0.129 | 0.008 | 0.041 |
| scaffold | L1_legacy_x_pa_morgan | historical | testlike | 0.050 | -0.004 | 0.109 | 0.965 |
| scaffold | L1_legacy_x_pa_morgan | historical | official | -0.088 | -0.158 | -0.020 | 0.006 |
| scaffold | L2_legacy_x_sqrt_pa_morgan | nested | testlike | 0.099 | 0.050 | 0.149 | 1.000 |
| scaffold | L2_legacy_x_sqrt_pa_morgan | nested | official | -0.019 | -0.084 | 0.045 | 0.281 |
| scaffold | L2_legacy_x_sqrt_pa_morgan | historical | testlike | 0.071 | 0.018 | 0.123 | 0.998 |
| scaffold | L2_legacy_x_sqrt_pa_morgan | historical | official | -0.041 | -0.110 | 0.023 | 0.116 |
| scaffold | L3_legacy_x_sqrt_pa_both | nested | testlike | 0.104 | 0.053 | 0.159 | 1.000 |
| scaffold | L3_legacy_x_sqrt_pa_both | nested | official | -0.047 | -0.116 | 0.023 | 0.095 |
| scaffold | L3_legacy_x_sqrt_pa_both | historical | testlike | 0.053 | -0.004 | 0.113 | 0.965 |
| scaffold | L3_legacy_x_sqrt_pa_both | historical | official | -0.087 | -0.158 | -0.012 | 0.011 |
| scaffold | L4_rankmean_legacy_pool | nested | testlike | 0.044 | 0.010 | 0.081 | 0.993 |
| scaffold | L4_rankmean_legacy_pool | nested | official | -0.007 | -0.057 | 0.045 | 0.387 |
| scaffold | L4_rankmean_legacy_pool | historical | testlike | 0.054 | 0.014 | 0.096 | 0.995 |
| scaffold | L4_rankmean_legacy_pool | historical | official | 0.012 | -0.044 | 0.068 | 0.638 |
| scaffold | L5_gate_filter_tau0.2 | nested | testlike | 0.089 | 0.048 | 0.134 | 1.000 |
| scaffold | L5_gate_filter_tau0.2 | nested | official | 0.041 | -0.020 | 0.098 | 0.899 |
| scaffold | L5_gate_filter_tau0.2 | historical | testlike | 0.073 | 0.034 | 0.115 | 1.000 |
| scaffold | L5_gate_filter_tau0.2 | historical | official | 0.033 | -0.019 | 0.086 | 0.881 |
| scaffold | L5_gate_filter_tau0.3 | nested | testlike | 0.081 | 0.035 | 0.132 | 1.000 |
| scaffold | L5_gate_filter_tau0.3 | nested | official | 0.009 | -0.053 | 0.073 | 0.592 |
| scaffold | L5_gate_filter_tau0.3 | historical | testlike | 0.079 | 0.031 | 0.132 | 1.000 |
| scaffold | L5_gate_filter_tau0.3 | historical | official | 0.010 | -0.055 | 0.073 | 0.629 |
| scaffold | L5_gate_filter_tau0.4 | nested | testlike | 0.057 | 0.007 | 0.107 | 0.988 |
| scaffold | L5_gate_filter_tau0.4 | nested | official | -0.074 | -0.139 | -0.007 | 0.017 |
| scaffold | L5_gate_filter_tau0.4 | historical | testlike | 0.061 | 0.010 | 0.114 | 0.989 |
| scaffold | L5_gate_filter_tau0.4 | historical | official | -0.055 | -0.123 | 0.012 | 0.059 |
| scaffold | B_pool_lgbm | nested | testlike | 0.082 | 0.032 | 0.132 | 0.999 |
| scaffold | B_pool_lgbm | nested | official | -0.022 | -0.092 | 0.047 | 0.256 |
| scaffold | B_pool_lgbm | historical | testlike | 0.068 | 0.018 | 0.118 | 0.995 |
| scaffold | B_pool_lgbm | historical | official | -0.030 | -0.104 | 0.038 | 0.204 |

## Part 3 - blind swaps against the Sep 15 vector

Verification: `{"deployed_rows": 750, "deployed_calls": 49, "top49_equals_sep15_payload": true, "boundary_gap": 0.0013315579227696217, "blind_gate_ge_tau": 279, "blind_tau": 0.3}`

| variant | n_calls | overlap_sep15 | added | removed | mean_pa_morgan_called | mean_pa_morgan_added | mean_pa_morgan_removed |
|---|---|---|---|---|---|---|---|
| L0_legacy | 49 | 49 | 0 | 0 | 0.322 | – | – |
| L0_legacy | 60 | 49 | 11 | 0 | 0.323 | 0.331 | – |
| L1_legacy_x_pa_morgan | 49 | 21 | 28 | 28 | 0.559 | 0.594 | 0.177 |
| L1_legacy_x_pa_morgan | 60 | 22 | 38 | 27 | 0.543 | 0.565 | 0.172 |
| L2_legacy_x_sqrt_pa_morgan | 49 | 23 | 26 | 26 | 0.526 | 0.552 | 0.167 |
| L2_legacy_x_sqrt_pa_morgan | 60 | 27 | 33 | 22 | 0.491 | 0.514 | 0.147 |
| L3_legacy_x_sqrt_pa_both | 49 | 19 | 30 | 30 | 0.541 | 0.564 | 0.206 |
| L3_legacy_x_sqrt_pa_both | 60 | 21 | 39 | 28 | 0.522 | 0.540 | 0.197 |
| L4_rankmean_legacy_pool | 49 | 32 | 17 | 17 | 0.368 | 0.309 | 0.175 |
| L4_rankmean_legacy_pool | 60 | 34 | 26 | 15 | 0.358 | 0.327 | 0.184 |
| L5_gate_filter_tau0.3 | 49 | 22 | 27 | 27 | 0.478 | 0.456 | 0.172 |
| L5_gate_filter_tau0.3 | 60 | 22 | 38 | 27 | 0.483 | 0.470 | 0.172 |

Sep 15 gate profile: `{"n_positives": 49, "mean_p_active_morgan": 0.32153316914219016, "mean_p_active_chemeleon": 0.3036504669347515, "below_0.2": 18, "below_0.3": 27, "below_0.5": 41}`

## Part 4 - sensitivity to the stricter pseudo-negative set

| split_id | testlike_order | testlike_strict_order | identical | top_model_changed | spearman_of_ranks |
|---|---|---|---|---|---|
| umap | L2_legacy_x_sqrt_pa_morgan > L3_legacy_x_sqrt_pa_both > L5_gate_filter_tau0.2 > L4_rankmean_legacy_pool > L1_legacy_x_pa_morgan > L5_gate_filter_tau0.3 > B_pool_lgbm > L5_gate_filter_tau0.4 > L0_legacy | L2_legacy_x_sqrt_pa_morgan > L3_legacy_x_sqrt_pa_both > L5_gate_filter_tau0.2 > L1_legacy_x_pa_morgan > L4_rankmean_legacy_pool > L5_gate_filter_tau0.3 > B_pool_lgbm > L5_gate_filter_tau0.4 > L0_legacy | no | no | 0.983 |
| scaffold | L3_legacy_x_sqrt_pa_both > L2_legacy_x_sqrt_pa_morgan > L5_gate_filter_tau0.2 > L1_legacy_x_pa_morgan > B_pool_lgbm > L5_gate_filter_tau0.3 > L5_gate_filter_tau0.4 > L4_rankmean_legacy_pool > L0_legacy | L3_legacy_x_sqrt_pa_both > L2_legacy_x_sqrt_pa_morgan > L5_gate_filter_tau0.2 > L1_legacy_x_pa_morgan > L5_gate_filter_tau0.3 > B_pool_lgbm > L5_gate_filter_tau0.4 > L4_rankmean_legacy_pool > L0_legacy | no | no | 0.983 |

Paired group bootstrap under the stricter pseudo-negative set:

| split_id | model | operating_point | label_set | delta | ci_low | ci_high | p_delta_gt_0 |
|---|---|---|---|---|---|---|---|
| umap | L1_legacy_x_pa_morgan | nested | testlike | 0.038 | -0.003 | 0.081 | 0.967 |
| umap | L1_legacy_x_pa_morgan | nested | official | -0.113 | -0.161 | -0.058 | 0.000 |
| umap | L1_legacy_x_pa_morgan | historical | testlike | 0.016 | -0.032 | 0.068 | 0.742 |
| umap | L1_legacy_x_pa_morgan | historical | official | -0.142 | -0.206 | -0.070 | 0.001 |
| umap | L2_legacy_x_sqrt_pa_morgan | nested | testlike | 0.059 | 0.014 | 0.107 | 0.996 |
| umap | L2_legacy_x_sqrt_pa_morgan | nested | official | -0.066 | -0.118 | -0.007 | 0.013 |
| umap | L2_legacy_x_sqrt_pa_morgan | historical | testlike | 0.029 | -0.015 | 0.074 | 0.895 |
| umap | L2_legacy_x_sqrt_pa_morgan | historical | official | -0.107 | -0.165 | -0.045 | 0.001 |
| umap | L3_legacy_x_sqrt_pa_both | nested | testlike | 0.045 | -0.011 | 0.103 | 0.948 |
| umap | L3_legacy_x_sqrt_pa_both | nested | official | -0.115 | -0.180 | -0.041 | 0.002 |
| umap | L3_legacy_x_sqrt_pa_both | historical | testlike | 0.016 | -0.031 | 0.067 | 0.750 |
| umap | L3_legacy_x_sqrt_pa_both | historical | official | -0.150 | -0.209 | -0.084 | 0.001 |
| umap | L4_rankmean_legacy_pool | nested | testlike | 0.038 | 0.004 | 0.072 | 0.987 |
| umap | L4_rankmean_legacy_pool | nested | official | -0.022 | -0.076 | 0.040 | 0.234 |
| umap | L4_rankmean_legacy_pool | historical | testlike | 0.030 | -0.014 | 0.075 | 0.915 |
| umap | L4_rankmean_legacy_pool | historical | official | -0.041 | -0.102 | 0.031 | 0.115 |
| umap | L5_gate_filter_tau0.2 | nested | testlike | 0.042 | 0.008 | 0.078 | 0.993 |
| umap | L5_gate_filter_tau0.2 | nested | official | -0.030 | -0.074 | 0.018 | 0.107 |
| umap | L5_gate_filter_tau0.2 | historical | testlike | 0.020 | -0.013 | 0.055 | 0.877 |
| umap | L5_gate_filter_tau0.2 | historical | official | -0.051 | -0.097 | -0.001 | 0.024 |
| umap | L5_gate_filter_tau0.3 | nested | testlike | 0.036 | -0.004 | 0.077 | 0.962 |
| umap | L5_gate_filter_tau0.3 | nested | official | -0.071 | -0.125 | -0.016 | 0.006 |
| umap | L5_gate_filter_tau0.3 | historical | testlike | 0.013 | -0.031 | 0.056 | 0.721 |
| umap | L5_gate_filter_tau0.3 | historical | official | -0.100 | -0.159 | -0.037 | 0.001 |
| umap | L5_gate_filter_tau0.4 | nested | testlike | 0.019 | -0.028 | 0.068 | 0.805 |
| umap | L5_gate_filter_tau0.4 | nested | official | -0.117 | -0.171 | -0.053 | 0.001 |
| umap | L5_gate_filter_tau0.4 | historical | testlike | 0.022 | -0.015 | 0.067 | 0.867 |
| umap | L5_gate_filter_tau0.4 | historical | official | -0.114 | -0.166 | -0.055 | 0.001 |
| umap | B_pool_lgbm | nested | testlike | 0.033 | -0.017 | 0.086 | 0.888 |
| umap | B_pool_lgbm | nested | official | -0.100 | -0.160 | -0.027 | 0.005 |
| umap | B_pool_lgbm | historical | testlike | 0.022 | -0.028 | 0.078 | 0.807 |
| umap | B_pool_lgbm | historical | official | -0.116 | -0.182 | -0.043 | 0.001 |
| scaffold | L1_legacy_x_pa_morgan | nested | testlike | 0.089 | 0.035 | 0.142 | 1.000 |
| scaffold | L1_legacy_x_pa_morgan | nested | official | -0.061 | -0.129 | 0.008 | 0.041 |
| scaffold | L1_legacy_x_pa_morgan | historical | testlike | 0.051 | -0.006 | 0.106 | 0.962 |
| scaffold | L1_legacy_x_pa_morgan | historical | official | -0.088 | -0.158 | -0.020 | 0.006 |
| scaffold | L2_legacy_x_sqrt_pa_morgan | nested | testlike | 0.100 | 0.051 | 0.150 | 1.000 |
| scaffold | L2_legacy_x_sqrt_pa_morgan | nested | official | -0.019 | -0.084 | 0.045 | 0.281 |
| scaffold | L2_legacy_x_sqrt_pa_morgan | historical | testlike | 0.073 | 0.023 | 0.126 | 0.998 |
| scaffold | L2_legacy_x_sqrt_pa_morgan | historical | official | -0.041 | -0.110 | 0.023 | 0.116 |
| scaffold | L3_legacy_x_sqrt_pa_both | nested | testlike | 0.107 | 0.051 | 0.159 | 1.000 |
| scaffold | L3_legacy_x_sqrt_pa_both | nested | official | -0.047 | -0.116 | 0.023 | 0.095 |
| scaffold | L3_legacy_x_sqrt_pa_both | historical | testlike | 0.054 | -0.007 | 0.112 | 0.958 |
| scaffold | L3_legacy_x_sqrt_pa_both | historical | official | -0.087 | -0.158 | -0.012 | 0.011 |
| scaffold | L4_rankmean_legacy_pool | nested | testlike | 0.047 | 0.010 | 0.086 | 0.992 |
| scaffold | L4_rankmean_legacy_pool | nested | official | -0.007 | -0.057 | 0.045 | 0.387 |
| scaffold | L4_rankmean_legacy_pool | historical | testlike | 0.058 | 0.016 | 0.099 | 0.999 |
| scaffold | L4_rankmean_legacy_pool | historical | official | 0.012 | -0.044 | 0.068 | 0.638 |
| scaffold | L5_gate_filter_tau0.2 | nested | testlike | 0.092 | 0.046 | 0.137 | 1.000 |
| scaffold | L5_gate_filter_tau0.2 | nested | official | 0.041 | -0.020 | 0.098 | 0.899 |
| scaffold | L5_gate_filter_tau0.2 | historical | testlike | 0.076 | 0.035 | 0.119 | 1.000 |
| scaffold | L5_gate_filter_tau0.2 | historical | official | 0.033 | -0.019 | 0.086 | 0.881 |
| scaffold | L5_gate_filter_tau0.3 | nested | testlike | 0.083 | 0.033 | 0.134 | 1.000 |
| scaffold | L5_gate_filter_tau0.3 | nested | official | 0.009 | -0.053 | 0.073 | 0.592 |
| scaffold | L5_gate_filter_tau0.3 | historical | testlike | 0.082 | 0.031 | 0.134 | 1.000 |
| scaffold | L5_gate_filter_tau0.3 | historical | official | 0.010 | -0.055 | 0.073 | 0.629 |
| scaffold | L5_gate_filter_tau0.4 | nested | testlike | 0.072 | 0.020 | 0.124 | 0.996 |
| scaffold | L5_gate_filter_tau0.4 | nested | official | -0.060 | -0.124 | 0.007 | 0.037 |
| scaffold | L5_gate_filter_tau0.4 | historical | testlike | 0.063 | 0.009 | 0.117 | 0.990 |
| scaffold | L5_gate_filter_tau0.4 | historical | official | -0.055 | -0.123 | 0.012 | 0.059 |
| scaffold | B_pool_lgbm | nested | testlike | 0.083 | 0.031 | 0.134 | 0.998 |
| scaffold | B_pool_lgbm | nested | official | -0.022 | -0.092 | 0.047 | 0.256 |
| scaffold | B_pool_lgbm | historical | testlike | 0.067 | 0.012 | 0.121 | 0.993 |
| scaffold | B_pool_lgbm | historical | official | -0.030 | -0.104 | 0.038 | 0.204 |

Transitions are few, share one hidden label set and one bootstrap LB estimate;
sign agreement here is calibration evidence about the score column, not proof
that any variant would improve the Live MCC.
