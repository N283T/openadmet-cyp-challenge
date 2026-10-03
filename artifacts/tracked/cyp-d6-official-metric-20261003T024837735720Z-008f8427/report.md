# D6 recipes under metrics that do not score pseudo-negative labels

## Mixture fraction identified from the published prevalence

| scored_prevalence | p_strong | p_weak | phi | interpretable |
|---|---|---|---|---|
| 0.050 | 0.216 | 0.000 | 0.231 | True |
| 0.071 | 0.216 | 0.000 | 0.326 | True |
| 0.090 | 0.216 | 0.000 | 0.416 | True |
| 0.050 | 0.216 | 0.005 | 0.214 | True |
| 0.071 | 0.216 | 0.005 | 0.311 | True |
| 0.090 | 0.216 | 0.005 | 0.403 | True |
| 0.050 | 0.216 | 0.006 | 0.210 | True |
| 0.071 | 0.216 | 0.006 | 0.308 | True |
| 0.090 | 0.216 | 0.006 | 0.400 | True |

## Rate optima

| split_id | metric | parameter | model | best_rate | best_mcc | best_calls_of_750 | at_grid_edge |
|---|---|---|---|---|---|---|---|
| umap | O | 0.050 | gate_x_tdi | 0.500 | 0.026 | 375 | True |
| umap | O | 0.050 | gate_x_tdi_chemeleon | 0.500 | 0.034 | 375 | True |
| umap | O | 0.050 | legacy_sep15_refit | 0.080 | 0.098 | 60 | False |
| umap | O | 0.050 | official_lgbm | 0.125 | 0.058 | 94 | False |
| umap | O | 0.050 | pool_chemeleon | 0.350 | 0.040 | 262 | False |
| umap | O | 0.050 | pool_lgbm | 0.150 | 0.039 | 112 | False |
| umap | O | 0.050 | rank_mean_gates | 0.500 | 0.029 | 375 | True |
| umap | O | 0.050 | l4_rank_mean | 0.150 | 0.070 | 112 | False |
| umap | O | 0.071 | gate_x_tdi | 0.500 | 0.030 | 375 | True |
| umap | O | 0.071 | gate_x_tdi_chemeleon | 0.500 | 0.039 | 375 | True |
| umap | O | 0.071 | legacy_sep15_refit | 0.080 | 0.113 | 60 | False |
| umap | O | 0.071 | official_lgbm | 0.125 | 0.067 | 94 | False |
| umap | O | 0.071 | pool_chemeleon | 0.350 | 0.047 | 262 | False |
| umap | O | 0.071 | pool_lgbm | 0.150 | 0.046 | 112 | False |
| umap | O | 0.071 | rank_mean_gates | 0.500 | 0.034 | 375 | True |
| umap | O | 0.071 | l4_rank_mean | 0.150 | 0.082 | 112 | False |
| umap | O | 0.090 | gate_x_tdi | 0.500 | 0.034 | 375 | True |
| umap | O | 0.090 | gate_x_tdi_chemeleon | 0.500 | 0.044 | 375 | True |
| umap | O | 0.090 | legacy_sep15_refit | 0.080 | 0.124 | 60 | False |
| umap | O | 0.090 | official_lgbm | 0.125 | 0.075 | 94 | False |
| umap | O | 0.090 | pool_chemeleon | 0.350 | 0.052 | 262 | False |
| umap | O | 0.090 | pool_lgbm | 0.150 | 0.052 | 112 | False |
| umap | O | 0.090 | rank_mean_gates | 0.500 | 0.038 | 375 | True |
| umap | O | 0.090 | l4_rank_mean | 0.150 | 0.091 | 112 | False |
| umap | M | 0.000 | gate_x_tdi | 0.010 | 0.000 | 8 | True |
| umap | M | 0.000 | gate_x_tdi_chemeleon | 0.010 | 0.000 | 8 | True |
| umap | M | 0.000 | legacy_sep15_refit | 0.010 | 0.000 | 8 | True |
| umap | M | 0.000 | official_lgbm | 0.010 | 0.000 | 8 | True |
| umap | M | 0.000 | pool_chemeleon | 0.010 | 0.000 | 8 | True |
| umap | M | 0.000 | pool_lgbm | 0.010 | 0.000 | 8 | True |
| umap | M | 0.000 | rank_mean_gates | 0.010 | 0.000 | 8 | True |
| umap | M | 0.000 | l4_rank_mean | 0.010 | 0.000 | 8 | True |
| umap | M | 0.100 | gate_x_tdi | 0.125 | 0.079 | 94 | False |
| umap | M | 0.100 | gate_x_tdi_chemeleon | 0.400 | 0.075 | 300 | False |
| umap | M | 0.100 | legacy_sep15_refit | 0.080 | 0.040 | 60 | False |
| umap | M | 0.100 | official_lgbm | 0.125 | 0.022 | 94 | False |
| umap | M | 0.100 | pool_chemeleon | 0.350 | 0.076 | 262 | False |
| umap | M | 0.100 | pool_lgbm | 0.175 | 0.104 | 131 | False |
| umap | M | 0.100 | rank_mean_gates | 0.125 | 0.086 | 94 | False |
| umap | M | 0.100 | l4_rank_mean | 0.125 | 0.078 | 94 | False |
| umap | M | 0.200 | gate_x_tdi | 0.150 | 0.096 | 112 | False |
| umap | M | 0.200 | gate_x_tdi_chemeleon | 0.400 | 0.097 | 300 | False |
| umap | M | 0.200 | legacy_sep15_refit | 0.080 | 0.058 | 60 | False |
| umap | M | 0.200 | official_lgbm | 0.125 | 0.032 | 94 | False |
| umap | M | 0.200 | pool_chemeleon | 0.350 | 0.099 | 262 | False |
| umap | M | 0.200 | pool_lgbm | 0.175 | 0.129 | 131 | False |
| umap | M | 0.200 | rank_mean_gates | 0.300 | 0.104 | 225 | False |
| umap | M | 0.200 | l4_rank_mean | 0.125 | 0.102 | 94 | False |
| umap | M | 0.210 | gate_x_tdi | 0.150 | 0.097 | 112 | False |
| umap | M | 0.210 | gate_x_tdi_chemeleon | 0.400 | 0.098 | 300 | False |
| umap | M | 0.210 | legacy_sep15_refit | 0.080 | 0.060 | 60 | False |
| umap | M | 0.210 | official_lgbm | 0.125 | 0.033 | 94 | False |
| umap | M | 0.210 | pool_chemeleon | 0.350 | 0.101 | 262 | False |
| umap | M | 0.210 | pool_lgbm | 0.175 | 0.130 | 131 | False |
| umap | M | 0.210 | rank_mean_gates | 0.300 | 0.105 | 225 | False |
| umap | M | 0.210 | l4_rank_mean | 0.125 | 0.104 | 94 | False |
| umap | M | 0.214 | gate_x_tdi | 0.150 | 0.097 | 112 | False |
| umap | M | 0.214 | gate_x_tdi_chemeleon | 0.400 | 0.099 | 300 | False |
| umap | M | 0.214 | legacy_sep15_refit | 0.080 | 0.060 | 60 | False |
| umap | M | 0.214 | official_lgbm | 0.125 | 0.034 | 94 | False |
| umap | M | 0.214 | pool_chemeleon | 0.350 | 0.101 | 262 | False |
| umap | M | 0.214 | pool_lgbm | 0.175 | 0.131 | 131 | False |
| umap | M | 0.214 | rank_mean_gates | 0.300 | 0.106 | 225 | False |
| umap | M | 0.214 | l4_rank_mean | 0.125 | 0.105 | 94 | False |
| umap | M | 0.231 | gate_x_tdi | 0.150 | 0.098 | 112 | False |
| umap | M | 0.231 | gate_x_tdi_chemeleon | 0.400 | 0.101 | 300 | False |
| umap | M | 0.231 | legacy_sep15_refit | 0.080 | 0.063 | 60 | False |
| umap | M | 0.231 | official_lgbm | 0.125 | 0.035 | 94 | False |
| umap | M | 0.231 | pool_chemeleon | 0.350 | 0.104 | 262 | False |
| umap | M | 0.231 | pool_lgbm | 0.175 | 0.133 | 131 | False |
| umap | M | 0.231 | rank_mean_gates | 0.300 | 0.108 | 225 | False |
| umap | M | 0.231 | l4_rank_mean | 0.150 | 0.108 | 112 | False |
| umap | M | 0.250 | gate_x_tdi | 0.150 | 0.099 | 112 | False |
| umap | M | 0.250 | gate_x_tdi_chemeleon | 0.400 | 0.104 | 300 | False |
| umap | M | 0.250 | legacy_sep15_refit | 0.080 | 0.066 | 60 | False |
| umap | M | 0.250 | official_lgbm | 0.125 | 0.037 | 94 | False |
| umap | M | 0.250 | pool_chemeleon | 0.350 | 0.106 | 262 | False |
| umap | M | 0.250 | pool_lgbm | 0.175 | 0.135 | 131 | False |
| umap | M | 0.250 | rank_mean_gates | 0.300 | 0.109 | 225 | False |
| umap | M | 0.250 | l4_rank_mean | 0.150 | 0.111 | 112 | False |
| umap | M | 0.300 | gate_x_tdi | 0.150 | 0.101 | 112 | False |
| umap | M | 0.300 | gate_x_tdi_chemeleon | 0.400 | 0.108 | 300 | False |
| umap | M | 0.300 | legacy_sep15_refit | 0.080 | 0.073 | 60 | False |
| umap | M | 0.300 | official_lgbm | 0.125 | 0.042 | 94 | False |
| umap | M | 0.300 | pool_chemeleon | 0.350 | 0.111 | 262 | False |
| umap | M | 0.300 | pool_lgbm | 0.250 | 0.139 | 188 | False |
| umap | M | 0.300 | rank_mean_gates | 0.300 | 0.113 | 225 | False |
| umap | M | 0.300 | l4_rank_mean | 0.150 | 0.117 | 112 | False |
| umap | M | 0.308 | gate_x_tdi | 0.150 | 0.102 | 112 | False |
| umap | M | 0.308 | gate_x_tdi_chemeleon | 0.400 | 0.109 | 300 | False |
| umap | M | 0.308 | legacy_sep15_refit | 0.080 | 0.074 | 60 | False |
| umap | M | 0.308 | official_lgbm | 0.125 | 0.042 | 94 | False |
| umap | M | 0.308 | pool_chemeleon | 0.350 | 0.111 | 262 | False |
| umap | M | 0.308 | pool_lgbm | 0.250 | 0.139 | 188 | False |
| umap | M | 0.308 | rank_mean_gates | 0.300 | 0.113 | 225 | False |
| umap | M | 0.308 | l4_rank_mean | 0.150 | 0.118 | 112 | False |
| umap | M | 0.311 | gate_x_tdi | 0.150 | 0.102 | 112 | False |
| umap | M | 0.311 | gate_x_tdi_chemeleon | 0.400 | 0.109 | 300 | False |
| umap | M | 0.311 | legacy_sep15_refit | 0.080 | 0.075 | 60 | False |
| umap | M | 0.311 | official_lgbm | 0.125 | 0.043 | 94 | False |
| umap | M | 0.311 | pool_chemeleon | 0.350 | 0.112 | 262 | False |
| umap | M | 0.311 | pool_lgbm | 0.250 | 0.140 | 188 | False |
| umap | M | 0.311 | rank_mean_gates | 0.300 | 0.113 | 225 | False |
| umap | M | 0.311 | l4_rank_mean | 0.150 | 0.118 | 112 | False |
| umap | M | 0.326 | gate_x_tdi | 0.150 | 0.102 | 112 | False |
| umap | M | 0.326 | gate_x_tdi_chemeleon | 0.400 | 0.110 | 300 | False |
| umap | M | 0.326 | legacy_sep15_refit | 0.080 | 0.077 | 60 | False |
| umap | M | 0.326 | official_lgbm | 0.125 | 0.044 | 94 | False |
| umap | M | 0.326 | pool_chemeleon | 0.350 | 0.113 | 262 | False |
| umap | M | 0.326 | pool_lgbm | 0.250 | 0.140 | 188 | False |
| umap | M | 0.326 | rank_mean_gates | 0.300 | 0.114 | 225 | False |
| umap | M | 0.326 | l4_rank_mean | 0.150 | 0.120 | 112 | False |
| umap | M | 0.400 | gate_x_tdi | 0.150 | 0.101 | 112 | False |
| umap | M | 0.400 | gate_x_tdi_chemeleon | 0.400 | 0.112 | 300 | False |
| umap | M | 0.400 | legacy_sep15_refit | 0.080 | 0.087 | 60 | False |
| umap | M | 0.400 | official_lgbm | 0.125 | 0.050 | 94 | False |
| umap | M | 0.400 | pool_chemeleon | 0.350 | 0.116 | 262 | False |
| umap | M | 0.400 | pool_lgbm | 0.250 | 0.142 | 188 | False |
| umap | M | 0.400 | rank_mean_gates | 0.300 | 0.115 | 225 | False |
| umap | M | 0.400 | l4_rank_mean | 0.150 | 0.126 | 112 | False |
| umap | M | 0.400 | gate_x_tdi | 0.150 | 0.101 | 112 | False |
| umap | M | 0.400 | gate_x_tdi_chemeleon | 0.400 | 0.113 | 300 | False |
| umap | M | 0.400 | legacy_sep15_refit | 0.080 | 0.087 | 60 | False |
| umap | M | 0.400 | official_lgbm | 0.125 | 0.050 | 94 | False |
| umap | M | 0.400 | pool_chemeleon | 0.350 | 0.116 | 262 | False |
| umap | M | 0.400 | pool_lgbm | 0.250 | 0.142 | 188 | False |
| umap | M | 0.400 | rank_mean_gates | 0.300 | 0.115 | 225 | False |
| umap | M | 0.400 | l4_rank_mean | 0.150 | 0.126 | 112 | False |
| umap | M | 0.403 | gate_x_tdi | 0.150 | 0.101 | 112 | False |
| umap | M | 0.403 | gate_x_tdi_chemeleon | 0.400 | 0.113 | 300 | False |
| umap | M | 0.403 | legacy_sep15_refit | 0.080 | 0.087 | 60 | False |
| umap | M | 0.403 | official_lgbm | 0.125 | 0.050 | 94 | False |
| umap | M | 0.403 | pool_chemeleon | 0.350 | 0.116 | 262 | False |
| umap | M | 0.403 | pool_lgbm | 0.250 | 0.141 | 188 | False |
| umap | M | 0.403 | rank_mean_gates | 0.300 | 0.115 | 225 | False |
| umap | M | 0.403 | l4_rank_mean | 0.150 | 0.126 | 112 | False |
| umap | M | 0.416 | gate_x_tdi | 0.150 | 0.100 | 112 | False |
| umap | M | 0.416 | gate_x_tdi_chemeleon | 0.400 | 0.113 | 300 | False |
| umap | M | 0.416 | legacy_sep15_refit | 0.080 | 0.089 | 60 | False |
| umap | M | 0.416 | official_lgbm | 0.125 | 0.052 | 94 | False |
| umap | M | 0.416 | pool_chemeleon | 0.350 | 0.116 | 262 | False |
| umap | M | 0.416 | pool_lgbm | 0.250 | 0.141 | 188 | False |
| umap | M | 0.416 | rank_mean_gates | 0.300 | 0.115 | 225 | False |
| umap | M | 0.416 | l4_rank_mean | 0.150 | 0.127 | 112 | False |
| umap | M | 0.500 | gate_x_tdi | 0.500 | 0.097 | 375 | True |
| umap | M | 0.500 | gate_x_tdi_chemeleon | 0.400 | 0.112 | 300 | False |
| umap | M | 0.500 | legacy_sep15_refit | 0.080 | 0.100 | 60 | False |
| umap | M | 0.500 | official_lgbm | 0.125 | 0.059 | 94 | False |
| umap | M | 0.500 | pool_chemeleon | 0.350 | 0.117 | 262 | False |
| umap | M | 0.500 | pool_lgbm | 0.250 | 0.138 | 188 | False |
| umap | M | 0.500 | rank_mean_gates | 0.300 | 0.112 | 225 | False |
| umap | M | 0.500 | l4_rank_mean | 0.150 | 0.132 | 112 | False |
| umap | M | 0.600 | gate_x_tdi | 0.500 | 0.093 | 375 | True |
| umap | M | 0.600 | gate_x_tdi_chemeleon | 0.400 | 0.108 | 300 | False |
| umap | M | 0.600 | legacy_sep15_refit | 0.080 | 0.112 | 60 | False |
| umap | M | 0.600 | official_lgbm | 0.125 | 0.067 | 94 | False |
| umap | M | 0.600 | pool_chemeleon | 0.350 | 0.114 | 262 | False |
| umap | M | 0.600 | pool_lgbm | 0.250 | 0.131 | 188 | False |
| umap | M | 0.600 | rank_mean_gates | 0.300 | 0.105 | 225 | False |
| umap | M | 0.600 | l4_rank_mean | 0.150 | 0.134 | 112 | False |
| umap | M | 0.700 | gate_x_tdi | 0.500 | 0.086 | 375 | True |
| umap | M | 0.700 | gate_x_tdi_chemeleon | 0.400 | 0.101 | 300 | False |
| umap | M | 0.700 | legacy_sep15_refit | 0.080 | 0.125 | 60 | False |
| umap | M | 0.700 | official_lgbm | 0.125 | 0.076 | 94 | False |
| umap | M | 0.700 | pool_chemeleon | 0.350 | 0.108 | 262 | False |
| umap | M | 0.700 | pool_lgbm | 0.250 | 0.120 | 188 | False |
| umap | M | 0.700 | rank_mean_gates | 0.500 | 0.096 | 375 | True |
| umap | M | 0.700 | l4_rank_mean | 0.150 | 0.135 | 112 | False |
| umap | M | 0.800 | gate_x_tdi | 0.500 | 0.076 | 375 | True |
| umap | M | 0.800 | gate_x_tdi_chemeleon | 0.400 | 0.091 | 300 | False |
| umap | M | 0.800 | legacy_sep15_refit | 0.080 | 0.137 | 60 | False |
| umap | M | 0.800 | official_lgbm | 0.125 | 0.085 | 94 | False |
| umap | M | 0.800 | pool_chemeleon | 0.350 | 0.099 | 262 | False |
| umap | M | 0.800 | pool_lgbm | 0.250 | 0.107 | 188 | False |
| umap | M | 0.800 | rank_mean_gates | 0.500 | 0.085 | 375 | True |
| umap | M | 0.800 | l4_rank_mean | 0.150 | 0.133 | 112 | False |
| umap | M | 0.900 | gate_x_tdi | 0.500 | 0.064 | 375 | True |
| umap | M | 0.900 | gate_x_tdi_chemeleon | 0.500 | 0.078 | 375 | True |
| umap | M | 0.900 | legacy_sep15_refit | 0.080 | 0.150 | 60 | False |
| umap | M | 0.900 | official_lgbm | 0.125 | 0.094 | 94 | False |
| umap | M | 0.900 | pool_chemeleon | 0.350 | 0.088 | 262 | False |
| umap | M | 0.900 | pool_lgbm | 0.250 | 0.090 | 188 | False |
| umap | M | 0.900 | rank_mean_gates | 0.500 | 0.071 | 375 | True |
| umap | M | 0.900 | l4_rank_mean | 0.150 | 0.131 | 112 | False |
| umap | M | 1.000 | gate_x_tdi | 0.500 | 0.049 | 375 | True |
| umap | M | 1.000 | gate_x_tdi_chemeleon | 0.500 | 0.064 | 375 | True |
| umap | M | 1.000 | legacy_sep15_refit | 0.080 | 0.163 | 60 | False |
| umap | M | 1.000 | official_lgbm | 0.125 | 0.103 | 94 | False |
| umap | M | 1.000 | pool_chemeleon | 0.350 | 0.075 | 262 | False |
| umap | M | 1.000 | pool_lgbm | 0.150 | 0.073 | 112 | False |
| umap | M | 1.000 | rank_mean_gates | 0.500 | 0.054 | 375 | True |
| umap | M | 1.000 | l4_rank_mean | 0.150 | 0.126 | 112 | False |
| scaffold | O | 0.050 | gate_x_tdi | 0.125 | 0.055 | 94 | False |
| scaffold | O | 0.050 | gate_x_tdi_chemeleon | 0.400 | 0.057 | 300 | False |
| scaffold | O | 0.050 | legacy_sep15_refit | 0.125 | 0.120 | 94 | False |
| scaffold | O | 0.050 | official_lgbm | 0.030 | 0.067 | 22 | False |
| scaffold | O | 0.050 | pool_chemeleon | 0.250 | 0.056 | 188 | False |
| scaffold | O | 0.050 | pool_lgbm | 0.070 | 0.067 | 53 | False |
| scaffold | O | 0.050 | rank_mean_gates | 0.200 | 0.069 | 150 | False |
| scaffold | O | 0.050 | l4_rank_mean | 0.275 | 0.099 | 206 | False |
| scaffold | O | 0.071 | gate_x_tdi | 0.125 | 0.064 | 94 | False |
| scaffold | O | 0.071 | gate_x_tdi_chemeleon | 0.400 | 0.067 | 300 | False |
| scaffold | O | 0.071 | legacy_sep15_refit | 0.125 | 0.139 | 94 | False |
| scaffold | O | 0.071 | official_lgbm | 0.030 | 0.078 | 22 | False |
| scaffold | O | 0.071 | pool_chemeleon | 0.250 | 0.066 | 188 | False |
| scaffold | O | 0.071 | pool_lgbm | 0.070 | 0.078 | 53 | False |
| scaffold | O | 0.071 | rank_mean_gates | 0.200 | 0.080 | 150 | False |
| scaffold | O | 0.071 | l4_rank_mean | 0.275 | 0.116 | 206 | False |
| scaffold | O | 0.090 | gate_x_tdi | 0.125 | 0.072 | 94 | False |
| scaffold | O | 0.090 | gate_x_tdi_chemeleon | 0.400 | 0.075 | 300 | False |
| scaffold | O | 0.090 | legacy_sep15_refit | 0.125 | 0.153 | 94 | False |
| scaffold | O | 0.090 | official_lgbm | 0.030 | 0.085 | 22 | False |
| scaffold | O | 0.090 | pool_chemeleon | 0.250 | 0.074 | 188 | False |
| scaffold | O | 0.090 | pool_lgbm | 0.070 | 0.087 | 53 | False |
| scaffold | O | 0.090 | rank_mean_gates | 0.200 | 0.089 | 150 | False |
| scaffold | O | 0.090 | l4_rank_mean | 0.275 | 0.129 | 206 | False |
| scaffold | M | 0.000 | gate_x_tdi | 0.010 | 0.000 | 8 | True |
| scaffold | M | 0.000 | gate_x_tdi_chemeleon | 0.010 | 0.000 | 8 | True |
| scaffold | M | 0.000 | legacy_sep15_refit | 0.010 | 0.000 | 8 | True |
| scaffold | M | 0.000 | official_lgbm | 0.010 | 0.000 | 8 | True |
| scaffold | M | 0.000 | pool_chemeleon | 0.010 | 0.000 | 8 | True |
| scaffold | M | 0.000 | pool_lgbm | 0.010 | 0.000 | 8 | True |
| scaffold | M | 0.000 | rank_mean_gates | 0.010 | 0.000 | 8 | True |
| scaffold | M | 0.000 | l4_rank_mean | 0.010 | 0.000 | 8 | True |
| scaffold | M | 0.100 | gate_x_tdi | 0.200 | 0.099 | 150 | False |
| scaffold | M | 0.100 | gate_x_tdi_chemeleon | 0.150 | 0.101 | 112 | False |
| scaffold | M | 0.100 | legacy_sep15_refit | 0.125 | 0.053 | 94 | False |
| scaffold | M | 0.100 | official_lgbm | 0.020 | 0.036 | 15 | False |
| scaffold | M | 0.100 | pool_chemeleon | 0.250 | 0.095 | 188 | False |
| scaffold | M | 0.100 | pool_lgbm | 0.070 | 0.122 | 53 | False |
| scaffold | M | 0.100 | rank_mean_gates | 0.200 | 0.118 | 150 | False |
| scaffold | M | 0.100 | l4_rank_mean | 0.275 | 0.091 | 206 | False |
| scaffold | M | 0.200 | gate_x_tdi | 0.200 | 0.126 | 150 | False |
| scaffold | M | 0.200 | gate_x_tdi_chemeleon | 0.225 | 0.128 | 169 | False |
| scaffold | M | 0.200 | legacy_sep15_refit | 0.125 | 0.076 | 94 | False |
| scaffold | M | 0.200 | official_lgbm | 0.020 | 0.050 | 15 | False |
| scaffold | M | 0.200 | pool_chemeleon | 0.250 | 0.122 | 188 | False |
| scaffold | M | 0.200 | pool_lgbm | 0.070 | 0.149 | 53 | False |
| scaffold | M | 0.200 | rank_mean_gates | 0.200 | 0.150 | 150 | False |
| scaffold | M | 0.200 | l4_rank_mean | 0.275 | 0.123 | 206 | False |
| scaffold | M | 0.210 | gate_x_tdi | 0.200 | 0.128 | 150 | False |
| scaffold | M | 0.210 | gate_x_tdi_chemeleon | 0.225 | 0.130 | 169 | False |
| scaffold | M | 0.210 | legacy_sep15_refit | 0.125 | 0.078 | 94 | False |
| scaffold | M | 0.210 | official_lgbm | 0.020 | 0.051 | 15 | False |
| scaffold | M | 0.210 | pool_chemeleon | 0.250 | 0.124 | 188 | False |
| scaffold | M | 0.210 | pool_lgbm | 0.070 | 0.150 | 53 | False |
| scaffold | M | 0.210 | rank_mean_gates | 0.200 | 0.152 | 150 | False |
| scaffold | M | 0.210 | l4_rank_mean | 0.275 | 0.125 | 206 | False |
| scaffold | M | 0.214 | gate_x_tdi | 0.200 | 0.129 | 150 | False |
| scaffold | M | 0.214 | gate_x_tdi_chemeleon | 0.225 | 0.131 | 169 | False |
| scaffold | M | 0.214 | legacy_sep15_refit | 0.125 | 0.079 | 94 | False |
| scaffold | M | 0.214 | official_lgbm | 0.020 | 0.051 | 15 | False |
| scaffold | M | 0.214 | pool_chemeleon | 0.250 | 0.125 | 188 | False |
| scaffold | M | 0.214 | pool_lgbm | 0.300 | 0.151 | 225 | False |
| scaffold | M | 0.214 | rank_mean_gates | 0.200 | 0.153 | 150 | False |
| scaffold | M | 0.214 | l4_rank_mean | 0.275 | 0.126 | 206 | False |
| scaffold | M | 0.231 | gate_x_tdi | 0.200 | 0.131 | 150 | False |
| scaffold | M | 0.231 | gate_x_tdi_chemeleon | 0.225 | 0.133 | 169 | False |
| scaffold | M | 0.231 | legacy_sep15_refit | 0.125 | 0.083 | 94 | False |
| scaffold | M | 0.231 | official_lgbm | 0.020 | 0.053 | 15 | False |
| scaffold | M | 0.231 | pool_chemeleon | 0.250 | 0.128 | 188 | False |
| scaffold | M | 0.231 | pool_lgbm | 0.300 | 0.154 | 225 | False |
| scaffold | M | 0.231 | rank_mean_gates | 0.200 | 0.156 | 150 | False |
| scaffold | M | 0.231 | l4_rank_mean | 0.275 | 0.130 | 206 | False |
| scaffold | M | 0.250 | gate_x_tdi | 0.200 | 0.134 | 150 | False |
| scaffold | M | 0.250 | gate_x_tdi_chemeleon | 0.225 | 0.136 | 169 | False |
| scaffold | M | 0.250 | legacy_sep15_refit | 0.125 | 0.086 | 94 | False |
| scaffold | M | 0.250 | official_lgbm | 0.020 | 0.055 | 15 | False |
| scaffold | M | 0.250 | pool_chemeleon | 0.250 | 0.131 | 188 | False |
| scaffold | M | 0.250 | pool_lgbm | 0.300 | 0.157 | 225 | False |
| scaffold | M | 0.250 | rank_mean_gates | 0.200 | 0.159 | 150 | False |
| scaffold | M | 0.250 | l4_rank_mean | 0.275 | 0.134 | 206 | False |
| scaffold | M | 0.300 | gate_x_tdi | 0.200 | 0.139 | 150 | False |
| scaffold | M | 0.300 | gate_x_tdi_chemeleon | 0.225 | 0.142 | 169 | False |
| scaffold | M | 0.300 | legacy_sep15_refit | 0.125 | 0.095 | 94 | False |
| scaffold | M | 0.300 | official_lgbm | 0.020 | 0.060 | 15 | False |
| scaffold | M | 0.300 | pool_chemeleon | 0.250 | 0.137 | 188 | False |
| scaffold | M | 0.300 | pool_lgbm | 0.300 | 0.163 | 225 | False |
| scaffold | M | 0.300 | rank_mean_gates | 0.200 | 0.166 | 150 | False |
| scaffold | M | 0.300 | l4_rank_mean | 0.275 | 0.143 | 206 | False |
| scaffold | M | 0.308 | gate_x_tdi | 0.200 | 0.140 | 150 | False |
| scaffold | M | 0.308 | gate_x_tdi_chemeleon | 0.225 | 0.143 | 169 | False |
| scaffold | M | 0.308 | legacy_sep15_refit | 0.125 | 0.097 | 94 | False |
| scaffold | M | 0.308 | official_lgbm | 0.020 | 0.061 | 15 | False |
| scaffold | M | 0.308 | pool_chemeleon | 0.250 | 0.137 | 188 | False |
| scaffold | M | 0.308 | pool_lgbm | 0.300 | 0.164 | 225 | False |
| scaffold | M | 0.308 | rank_mean_gates | 0.200 | 0.166 | 150 | False |
| scaffold | M | 0.308 | l4_rank_mean | 0.275 | 0.145 | 206 | False |
| scaffold | M | 0.311 | gate_x_tdi | 0.200 | 0.140 | 150 | False |
| scaffold | M | 0.311 | gate_x_tdi_chemeleon | 0.225 | 0.143 | 169 | False |
| scaffold | M | 0.311 | legacy_sep15_refit | 0.125 | 0.097 | 94 | False |
| scaffold | M | 0.311 | official_lgbm | 0.020 | 0.061 | 15 | False |
| scaffold | M | 0.311 | pool_chemeleon | 0.250 | 0.138 | 188 | False |
| scaffold | M | 0.311 | pool_lgbm | 0.300 | 0.164 | 225 | False |
| scaffold | M | 0.311 | rank_mean_gates | 0.200 | 0.167 | 150 | False |
| scaffold | M | 0.311 | l4_rank_mean | 0.275 | 0.145 | 206 | False |
| scaffold | M | 0.326 | gate_x_tdi | 0.200 | 0.141 | 150 | False |
| scaffold | M | 0.326 | gate_x_tdi_chemeleon | 0.225 | 0.144 | 169 | False |
| scaffold | M | 0.326 | legacy_sep15_refit | 0.125 | 0.100 | 94 | False |
| scaffold | M | 0.326 | official_lgbm | 0.020 | 0.063 | 15 | False |
| scaffold | M | 0.326 | pool_chemeleon | 0.250 | 0.139 | 188 | False |
| scaffold | M | 0.326 | pool_lgbm | 0.300 | 0.166 | 225 | False |
| scaffold | M | 0.326 | rank_mean_gates | 0.200 | 0.168 | 150 | False |
| scaffold | M | 0.326 | l4_rank_mean | 0.275 | 0.148 | 206 | False |
| scaffold | M | 0.400 | gate_x_tdi | 0.200 | 0.144 | 150 | False |
| scaffold | M | 0.400 | gate_x_tdi_chemeleon | 0.225 | 0.148 | 169 | False |
| scaffold | M | 0.400 | legacy_sep15_refit | 0.125 | 0.112 | 94 | False |
| scaffold | M | 0.400 | official_lgbm | 0.020 | 0.069 | 15 | False |
| scaffold | M | 0.400 | pool_chemeleon | 0.250 | 0.143 | 188 | False |
| scaffold | M | 0.400 | pool_lgbm | 0.300 | 0.170 | 225 | False |
| scaffold | M | 0.400 | rank_mean_gates | 0.200 | 0.172 | 150 | False |
| scaffold | M | 0.400 | l4_rank_mean | 0.275 | 0.158 | 206 | False |
| scaffold | M | 0.400 | gate_x_tdi | 0.200 | 0.144 | 150 | False |
| scaffold | M | 0.400 | gate_x_tdi_chemeleon | 0.225 | 0.148 | 169 | False |
| scaffold | M | 0.400 | legacy_sep15_refit | 0.125 | 0.112 | 94 | False |
| scaffold | M | 0.400 | official_lgbm | 0.020 | 0.069 | 15 | False |
| scaffold | M | 0.400 | pool_chemeleon | 0.250 | 0.143 | 188 | False |
| scaffold | M | 0.400 | pool_lgbm | 0.300 | 0.170 | 225 | False |
| scaffold | M | 0.400 | rank_mean_gates | 0.200 | 0.172 | 150 | False |
| scaffold | M | 0.400 | l4_rank_mean | 0.275 | 0.158 | 206 | False |
| scaffold | M | 0.403 | gate_x_tdi | 0.200 | 0.144 | 150 | False |
| scaffold | M | 0.403 | gate_x_tdi_chemeleon | 0.225 | 0.148 | 169 | False |
| scaffold | M | 0.403 | legacy_sep15_refit | 0.125 | 0.113 | 94 | False |
| scaffold | M | 0.403 | official_lgbm | 0.020 | 0.069 | 15 | False |
| scaffold | M | 0.403 | pool_chemeleon | 0.250 | 0.143 | 188 | False |
| scaffold | M | 0.403 | pool_lgbm | 0.300 | 0.170 | 225 | False |
| scaffold | M | 0.403 | rank_mean_gates | 0.200 | 0.172 | 150 | False |
| scaffold | M | 0.403 | l4_rank_mean | 0.275 | 0.158 | 206 | False |
| scaffold | M | 0.416 | gate_x_tdi | 0.200 | 0.144 | 150 | False |
| scaffold | M | 0.416 | gate_x_tdi_chemeleon | 0.225 | 0.148 | 169 | False |
| scaffold | M | 0.416 | legacy_sep15_refit | 0.125 | 0.115 | 94 | False |
| scaffold | M | 0.416 | official_lgbm | 0.050 | 0.070 | 38 | False |
| scaffold | M | 0.416 | pool_chemeleon | 0.250 | 0.144 | 188 | False |
| scaffold | M | 0.416 | pool_lgbm | 0.300 | 0.170 | 225 | False |
| scaffold | M | 0.416 | rank_mean_gates | 0.200 | 0.172 | 150 | False |
| scaffold | M | 0.416 | l4_rank_mean | 0.275 | 0.160 | 206 | False |
| scaffold | M | 0.500 | gate_x_tdi | 0.200 | 0.144 | 150 | False |
| scaffold | M | 0.500 | gate_x_tdi_chemeleon | 0.225 | 0.148 | 169 | False |
| scaffold | M | 0.500 | legacy_sep15_refit | 0.125 | 0.128 | 94 | False |
| scaffold | M | 0.500 | official_lgbm | 0.050 | 0.077 | 38 | False |
| scaffold | M | 0.500 | pool_chemeleon | 0.250 | 0.144 | 188 | False |
| scaffold | M | 0.500 | pool_lgbm | 0.300 | 0.170 | 225 | False |
| scaffold | M | 0.500 | rank_mean_gates | 0.200 | 0.173 | 150 | False |
| scaffold | M | 0.500 | l4_rank_mean | 0.275 | 0.168 | 206 | False |
| scaffold | M | 0.600 | gate_x_tdi | 0.200 | 0.140 | 150 | False |
| scaffold | M | 0.600 | gate_x_tdi_chemeleon | 0.225 | 0.145 | 169 | False |
| scaffold | M | 0.600 | legacy_sep15_refit | 0.125 | 0.144 | 94 | False |
| scaffold | M | 0.600 | official_lgbm | 0.050 | 0.084 | 38 | False |
| scaffold | M | 0.600 | pool_chemeleon | 0.250 | 0.142 | 188 | False |
| scaffold | M | 0.600 | pool_lgbm | 0.300 | 0.165 | 225 | False |
| scaffold | M | 0.600 | rank_mean_gates | 0.200 | 0.169 | 150 | False |
| scaffold | M | 0.600 | l4_rank_mean | 0.275 | 0.175 | 206 | False |
| scaffold | M | 0.700 | gate_x_tdi | 0.200 | 0.133 | 150 | False |
| scaffold | M | 0.700 | gate_x_tdi_chemeleon | 0.225 | 0.138 | 169 | False |
| scaffold | M | 0.700 | legacy_sep15_refit | 0.125 | 0.158 | 94 | False |
| scaffold | M | 0.700 | official_lgbm | 0.050 | 0.091 | 38 | False |
| scaffold | M | 0.700 | pool_chemeleon | 0.250 | 0.136 | 188 | False |
| scaffold | M | 0.700 | pool_lgbm | 0.300 | 0.156 | 225 | False |
| scaffold | M | 0.700 | rank_mean_gates | 0.200 | 0.162 | 150 | False |
| scaffold | M | 0.700 | l4_rank_mean | 0.275 | 0.179 | 206 | False |
| scaffold | M | 0.800 | gate_x_tdi | 0.200 | 0.124 | 150 | False |
| scaffold | M | 0.800 | gate_x_tdi_chemeleon | 0.400 | 0.130 | 300 | False |
| scaffold | M | 0.800 | legacy_sep15_refit | 0.125 | 0.173 | 94 | False |
| scaffold | M | 0.800 | official_lgbm | 0.050 | 0.097 | 38 | False |
| scaffold | M | 0.800 | pool_chemeleon | 0.250 | 0.128 | 188 | False |
| scaffold | M | 0.800 | pool_lgbm | 0.300 | 0.145 | 225 | False |
| scaffold | M | 0.800 | rank_mean_gates | 0.200 | 0.153 | 150 | False |
| scaffold | M | 0.800 | l4_rank_mean | 0.275 | 0.182 | 206 | False |
| scaffold | M | 0.900 | gate_x_tdi | 0.200 | 0.112 | 150 | False |
| scaffold | M | 0.900 | gate_x_tdi_chemeleon | 0.400 | 0.120 | 300 | False |
| scaffold | M | 0.900 | legacy_sep15_refit | 0.125 | 0.188 | 94 | False |
| scaffold | M | 0.900 | official_lgbm | 0.050 | 0.103 | 38 | False |
| scaffold | M | 0.900 | pool_chemeleon | 0.250 | 0.118 | 188 | False |
| scaffold | M | 0.900 | pool_lgbm | 0.070 | 0.131 | 53 | False |
| scaffold | M | 0.900 | rank_mean_gates | 0.200 | 0.141 | 150 | False |
| scaffold | M | 0.900 | l4_rank_mean | 0.275 | 0.182 | 206 | False |
| scaffold | M | 1.000 | gate_x_tdi | 0.125 | 0.100 | 94 | False |
| scaffold | M | 1.000 | gate_x_tdi_chemeleon | 0.400 | 0.108 | 300 | False |
| scaffold | M | 1.000 | legacy_sep15_refit | 0.125 | 0.202 | 94 | False |
| scaffold | M | 1.000 | official_lgbm | 0.030 | 0.110 | 22 | False |
| scaffold | M | 1.000 | pool_chemeleon | 0.250 | 0.105 | 188 | False |
| scaffold | M | 1.000 | pool_lgbm | 0.070 | 0.120 | 53 | False |
| scaffold | M | 1.000 | rank_mean_gates | 0.200 | 0.127 | 150 | False |
| scaffold | M | 1.000 | l4_rank_mean | 0.275 | 0.181 | 206 | False |

## Paired group bootstrap versus the official-only control

| split_id | metric | parameter | model | rate | delta | ci_low | ci_high | p_delta_gt_0 |
|---|---|---|---|---|---|---|---|---|
| umap | O | 0.050 | gate_x_tdi | 0.065 | -0.020 | -0.051 | 0.009 | 0.096 |
| umap | O | 0.050 | gate_x_tdi | 0.100 | -0.039 | -0.070 | -0.012 | 0.001 |
| umap | O | 0.050 | gate_x_tdi | 0.200 | -0.035 | -0.054 | -0.013 | 0.002 |
| umap | O | 0.050 | gate_x_tdi_chemeleon | 0.065 | -0.018 | -0.067 | 0.039 | 0.261 |
| umap | O | 0.050 | gate_x_tdi_chemeleon | 0.100 | -0.029 | -0.076 | 0.014 | 0.101 |
| umap | O | 0.050 | gate_x_tdi_chemeleon | 0.200 | -0.026 | -0.063 | 0.019 | 0.110 |
| umap | O | 0.050 | legacy_sep15_refit | 0.065 | 0.054 | 0.004 | 0.107 | 0.986 |
| umap | O | 0.050 | legacy_sep15_refit | 0.100 | 0.025 | -0.022 | 0.072 | 0.852 |
| umap | O | 0.050 | legacy_sep15_refit | 0.200 | 0.022 | -0.004 | 0.051 | 0.947 |
| umap | O | 0.050 | pool_chemeleon | 0.065 | -0.037 | -0.069 | -0.001 | 0.021 |
| umap | O | 0.050 | pool_chemeleon | 0.100 | -0.047 | -0.079 | -0.014 | 0.003 |
| umap | O | 0.050 | pool_chemeleon | 0.200 | -0.026 | -0.057 | 0.015 | 0.107 |
| umap | O | 0.050 | pool_lgbm | 0.065 | -0.018 | -0.054 | 0.022 | 0.172 |
| umap | O | 0.050 | pool_lgbm | 0.100 | -0.037 | -0.067 | -0.006 | 0.013 |
| umap | O | 0.050 | pool_lgbm | 0.200 | -0.018 | -0.051 | 0.020 | 0.166 |
| umap | O | 0.050 | rank_mean_gates | 0.065 | -0.017 | -0.056 | 0.027 | 0.212 |
| umap | O | 0.050 | rank_mean_gates | 0.100 | -0.039 | -0.072 | -0.011 | 0.007 |
| umap | O | 0.050 | rank_mean_gates | 0.200 | -0.031 | -0.063 | 0.005 | 0.050 |
| umap | O | 0.050 | l4_rank_mean | 0.065 | 0.026 | -0.007 | 0.062 | 0.941 |
| umap | O | 0.050 | l4_rank_mean | 0.100 | 0.013 | -0.014 | 0.047 | 0.802 |
| umap | O | 0.050 | l4_rank_mean | 0.200 | 0.010 | -0.022 | 0.042 | 0.724 |
| umap | O | 0.071 | gate_x_tdi | 0.065 | -0.024 | -0.059 | 0.010 | 0.096 |
| umap | O | 0.071 | gate_x_tdi | 0.100 | -0.045 | -0.080 | -0.014 | 0.001 |
| umap | O | 0.071 | gate_x_tdi | 0.200 | -0.041 | -0.063 | -0.016 | 0.002 |
| umap | O | 0.071 | gate_x_tdi_chemeleon | 0.065 | -0.021 | -0.078 | 0.045 | 0.261 |
| umap | O | 0.071 | gate_x_tdi_chemeleon | 0.100 | -0.034 | -0.087 | 0.017 | 0.101 |
| umap | O | 0.071 | gate_x_tdi_chemeleon | 0.200 | -0.031 | -0.074 | 0.022 | 0.110 |
| umap | O | 0.071 | legacy_sep15_refit | 0.065 | 0.062 | 0.005 | 0.122 | 0.985 |
| umap | O | 0.071 | legacy_sep15_refit | 0.100 | 0.029 | -0.025 | 0.083 | 0.851 |
| umap | O | 0.071 | legacy_sep15_refit | 0.200 | 0.026 | -0.005 | 0.059 | 0.947 |
| umap | O | 0.071 | pool_chemeleon | 0.065 | -0.043 | -0.080 | -0.001 | 0.021 |
| umap | O | 0.071 | pool_chemeleon | 0.100 | -0.055 | -0.091 | -0.016 | 0.003 |
| umap | O | 0.071 | pool_chemeleon | 0.200 | -0.031 | -0.067 | 0.017 | 0.107 |
| umap | O | 0.071 | pool_lgbm | 0.065 | -0.021 | -0.062 | 0.026 | 0.172 |
| umap | O | 0.071 | pool_lgbm | 0.100 | -0.043 | -0.078 | -0.007 | 0.013 |
| umap | O | 0.071 | pool_lgbm | 0.200 | -0.021 | -0.059 | 0.024 | 0.167 |
| umap | O | 0.071 | rank_mean_gates | 0.065 | -0.020 | -0.064 | 0.032 | 0.212 |
| umap | O | 0.071 | rank_mean_gates | 0.100 | -0.046 | -0.084 | -0.012 | 0.007 |
| umap | O | 0.071 | rank_mean_gates | 0.200 | -0.036 | -0.074 | 0.006 | 0.050 |
| umap | O | 0.071 | l4_rank_mean | 0.065 | 0.030 | -0.008 | 0.071 | 0.942 |
| umap | O | 0.071 | l4_rank_mean | 0.100 | 0.015 | -0.016 | 0.055 | 0.806 |
| umap | O | 0.071 | l4_rank_mean | 0.200 | 0.011 | -0.025 | 0.049 | 0.724 |
| umap | O | 0.090 | gate_x_tdi | 0.065 | -0.026 | -0.065 | 0.011 | 0.096 |
| umap | O | 0.090 | gate_x_tdi | 0.100 | -0.050 | -0.088 | -0.015 | 0.001 |
| umap | O | 0.090 | gate_x_tdi | 0.200 | -0.046 | -0.070 | -0.018 | 0.002 |
| umap | O | 0.090 | gate_x_tdi_chemeleon | 0.065 | -0.023 | -0.085 | 0.051 | 0.261 |
| umap | O | 0.090 | gate_x_tdi_chemeleon | 0.100 | -0.037 | -0.096 | 0.019 | 0.101 |
| umap | O | 0.090 | gate_x_tdi_chemeleon | 0.200 | -0.034 | -0.082 | 0.025 | 0.110 |
| umap | O | 0.090 | legacy_sep15_refit | 0.065 | 0.068 | 0.005 | 0.132 | 0.985 |
| umap | O | 0.090 | legacy_sep15_refit | 0.100 | 0.032 | -0.027 | 0.091 | 0.851 |
| umap | O | 0.090 | legacy_sep15_refit | 0.200 | 0.029 | -0.005 | 0.065 | 0.947 |
| umap | O | 0.090 | pool_chemeleon | 0.065 | -0.048 | -0.089 | -0.002 | 0.021 |
| umap | O | 0.090 | pool_chemeleon | 0.100 | -0.061 | -0.101 | -0.018 | 0.003 |
| umap | O | 0.090 | pool_chemeleon | 0.200 | -0.034 | -0.075 | 0.019 | 0.108 |
| umap | O | 0.090 | pool_lgbm | 0.065 | -0.023 | -0.069 | 0.029 | 0.172 |
| umap | O | 0.090 | pool_lgbm | 0.100 | -0.048 | -0.085 | -0.007 | 0.013 |
| umap | O | 0.090 | pool_lgbm | 0.200 | -0.024 | -0.066 | 0.026 | 0.167 |
| umap | O | 0.090 | rank_mean_gates | 0.065 | -0.022 | -0.070 | 0.035 | 0.213 |
| umap | O | 0.090 | rank_mean_gates | 0.100 | -0.051 | -0.093 | -0.014 | 0.007 |
| umap | O | 0.090 | rank_mean_gates | 0.200 | -0.040 | -0.082 | 0.007 | 0.050 |
| umap | O | 0.090 | l4_rank_mean | 0.065 | 0.033 | -0.008 | 0.079 | 0.942 |
| umap | O | 0.090 | l4_rank_mean | 0.100 | 0.017 | -0.017 | 0.060 | 0.808 |
| umap | O | 0.090 | l4_rank_mean | 0.200 | 0.013 | -0.028 | 0.055 | 0.724 |
| umap | M | 0.000 | gate_x_tdi | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | gate_x_tdi | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | gate_x_tdi | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | gate_x_tdi_chemeleon | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | gate_x_tdi_chemeleon | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | gate_x_tdi_chemeleon | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | legacy_sep15_refit | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | legacy_sep15_refit | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | legacy_sep15_refit | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | pool_chemeleon | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | pool_chemeleon | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | pool_chemeleon | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | pool_lgbm | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | pool_lgbm | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | pool_lgbm | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | rank_mean_gates | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | rank_mean_gates | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | rank_mean_gates | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | l4_rank_mean | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | l4_rank_mean | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.000 | l4_rank_mean | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | M | 0.100 | gate_x_tdi | 0.065 | 0.049 | 0.027 | 0.074 | 1.000 |
| umap | M | 0.100 | gate_x_tdi | 0.100 | 0.052 | 0.037 | 0.070 | 1.000 |
| umap | M | 0.100 | gate_x_tdi | 0.200 | 0.050 | 0.035 | 0.067 | 1.000 |
| umap | M | 0.100 | gate_x_tdi_chemeleon | 0.065 | 0.042 | 0.003 | 0.086 | 0.983 |
| umap | M | 0.100 | gate_x_tdi_chemeleon | 0.100 | 0.048 | 0.016 | 0.081 | 0.998 |
| umap | M | 0.100 | gate_x_tdi_chemeleon | 0.200 | 0.050 | 0.022 | 0.081 | 1.000 |
| umap | M | 0.100 | legacy_sep15_refit | 0.065 | 0.019 | -0.007 | 0.047 | 0.904 |
| umap | M | 0.100 | legacy_sep15_refit | 0.100 | 0.016 | -0.008 | 0.041 | 0.899 |
| umap | M | 0.100 | legacy_sep15_refit | 0.200 | 0.015 | -0.002 | 0.033 | 0.960 |
| umap | M | 0.100 | pool_chemeleon | 0.065 | 0.025 | 0.001 | 0.051 | 0.978 |
| umap | M | 0.100 | pool_chemeleon | 0.100 | 0.038 | 0.017 | 0.061 | 1.000 |
| umap | M | 0.100 | pool_chemeleon | 0.200 | 0.048 | 0.025 | 0.076 | 1.000 |
| umap | M | 0.100 | pool_lgbm | 0.065 | 0.064 | 0.036 | 0.094 | 1.000 |
| umap | M | 0.100 | pool_lgbm | 0.100 | 0.065 | 0.041 | 0.089 | 1.000 |
| umap | M | 0.100 | pool_lgbm | 0.200 | 0.078 | 0.047 | 0.111 | 1.000 |
| umap | M | 0.100 | rank_mean_gates | 0.065 | 0.055 | 0.028 | 0.087 | 1.000 |
| umap | M | 0.100 | rank_mean_gates | 0.100 | 0.061 | 0.040 | 0.082 | 1.000 |
| umap | M | 0.100 | rank_mean_gates | 0.200 | 0.061 | 0.036 | 0.085 | 1.000 |
| umap | M | 0.100 | l4_rank_mean | 0.065 | 0.048 | 0.024 | 0.077 | 1.000 |
| umap | M | 0.100 | l4_rank_mean | 0.100 | 0.056 | 0.035 | 0.079 | 1.000 |
| umap | M | 0.100 | l4_rank_mean | 0.200 | 0.052 | 0.029 | 0.076 | 1.000 |
| umap | M | 0.200 | gate_x_tdi | 0.065 | 0.054 | 0.027 | 0.084 | 1.000 |
| umap | M | 0.200 | gate_x_tdi | 0.100 | 0.057 | 0.038 | 0.079 | 1.000 |
| umap | M | 0.200 | gate_x_tdi | 0.200 | 0.058 | 0.038 | 0.078 | 1.000 |
| umap | M | 0.200 | gate_x_tdi_chemeleon | 0.065 | 0.047 | -0.004 | 0.103 | 0.966 |
| umap | M | 0.200 | gate_x_tdi_chemeleon | 0.100 | 0.055 | 0.012 | 0.098 | 0.995 |
| umap | M | 0.200 | gate_x_tdi_chemeleon | 0.200 | 0.059 | 0.021 | 0.101 | 0.998 |
| umap | M | 0.200 | legacy_sep15_refit | 0.065 | 0.028 | -0.008 | 0.066 | 0.920 |
| umap | M | 0.200 | legacy_sep15_refit | 0.100 | 0.022 | -0.012 | 0.058 | 0.897 |
| umap | M | 0.200 | legacy_sep15_refit | 0.200 | 0.021 | -0.003 | 0.046 | 0.960 |
| umap | M | 0.200 | pool_chemeleon | 0.065 | 0.025 | -0.007 | 0.060 | 0.944 |
| umap | M | 0.200 | pool_chemeleon | 0.100 | 0.040 | 0.013 | 0.071 | 0.999 |
| umap | M | 0.200 | pool_chemeleon | 0.200 | 0.056 | 0.024 | 0.094 | 1.000 |
| umap | M | 0.200 | pool_lgbm | 0.065 | 0.070 | 0.034 | 0.109 | 1.000 |
| umap | M | 0.200 | pool_lgbm | 0.100 | 0.071 | 0.040 | 0.102 | 1.000 |
| umap | M | 0.200 | pool_lgbm | 0.200 | 0.092 | 0.050 | 0.134 | 1.000 |
| umap | M | 0.200 | rank_mean_gates | 0.065 | 0.061 | 0.025 | 0.104 | 0.998 |
| umap | M | 0.200 | rank_mean_gates | 0.100 | 0.066 | 0.041 | 0.091 | 1.000 |
| umap | M | 0.200 | rank_mean_gates | 0.200 | 0.070 | 0.038 | 0.104 | 1.000 |
| umap | M | 0.200 | l4_rank_mean | 0.065 | 0.062 | 0.029 | 0.098 | 1.000 |
| umap | M | 0.200 | l4_rank_mean | 0.100 | 0.070 | 0.042 | 0.099 | 1.000 |
| umap | M | 0.200 | l4_rank_mean | 0.200 | 0.066 | 0.034 | 0.098 | 1.000 |
| umap | M | 0.210 | gate_x_tdi | 0.065 | 0.054 | 0.026 | 0.085 | 1.000 |
| umap | M | 0.210 | gate_x_tdi | 0.100 | 0.057 | 0.037 | 0.079 | 1.000 |
| umap | M | 0.210 | gate_x_tdi | 0.200 | 0.058 | 0.038 | 0.079 | 1.000 |
| umap | M | 0.210 | gate_x_tdi_chemeleon | 0.065 | 0.047 | -0.005 | 0.105 | 0.966 |
| umap | M | 0.210 | gate_x_tdi_chemeleon | 0.100 | 0.055 | 0.012 | 0.099 | 0.994 |
| umap | M | 0.210 | gate_x_tdi_chemeleon | 0.200 | 0.059 | 0.020 | 0.102 | 0.998 |
| umap | M | 0.210 | legacy_sep15_refit | 0.065 | 0.029 | -0.009 | 0.068 | 0.922 |
| umap | M | 0.210 | legacy_sep15_refit | 0.100 | 0.023 | -0.012 | 0.059 | 0.897 |
| umap | M | 0.210 | legacy_sep15_refit | 0.200 | 0.021 | -0.003 | 0.047 | 0.960 |
| umap | M | 0.210 | pool_chemeleon | 0.065 | 0.025 | -0.008 | 0.060 | 0.939 |
| umap | M | 0.210 | pool_chemeleon | 0.100 | 0.040 | 0.013 | 0.071 | 0.999 |
| umap | M | 0.210 | pool_chemeleon | 0.200 | 0.056 | 0.024 | 0.095 | 1.000 |
| umap | M | 0.210 | pool_lgbm | 0.065 | 0.070 | 0.033 | 0.109 | 1.000 |
| umap | M | 0.210 | pool_lgbm | 0.100 | 0.071 | 0.039 | 0.102 | 1.000 |
| umap | M | 0.210 | pool_lgbm | 0.200 | 0.092 | 0.050 | 0.135 | 1.000 |
| umap | M | 0.210 | rank_mean_gates | 0.065 | 0.061 | 0.024 | 0.105 | 0.998 |
| umap | M | 0.210 | rank_mean_gates | 0.100 | 0.066 | 0.040 | 0.092 | 1.000 |
| umap | M | 0.210 | rank_mean_gates | 0.200 | 0.070 | 0.037 | 0.105 | 1.000 |
| umap | M | 0.210 | l4_rank_mean | 0.065 | 0.062 | 0.029 | 0.100 | 1.000 |
| umap | M | 0.210 | l4_rank_mean | 0.100 | 0.070 | 0.042 | 0.100 | 1.000 |
| umap | M | 0.210 | l4_rank_mean | 0.200 | 0.066 | 0.035 | 0.099 | 1.000 |
| umap | M | 0.214 | gate_x_tdi | 0.065 | 0.054 | 0.026 | 0.085 | 1.000 |
| umap | M | 0.214 | gate_x_tdi | 0.100 | 0.057 | 0.037 | 0.079 | 1.000 |
| umap | M | 0.214 | gate_x_tdi | 0.200 | 0.058 | 0.038 | 0.079 | 1.000 |
| umap | M | 0.214 | gate_x_tdi_chemeleon | 0.065 | 0.047 | -0.005 | 0.105 | 0.965 |
| umap | M | 0.214 | gate_x_tdi_chemeleon | 0.100 | 0.055 | 0.012 | 0.099 | 0.994 |
| umap | M | 0.214 | gate_x_tdi_chemeleon | 0.200 | 0.060 | 0.020 | 0.102 | 0.998 |
| umap | M | 0.214 | legacy_sep15_refit | 0.065 | 0.029 | -0.009 | 0.068 | 0.922 |
| umap | M | 0.214 | legacy_sep15_refit | 0.100 | 0.023 | -0.012 | 0.060 | 0.896 |
| umap | M | 0.214 | legacy_sep15_refit | 0.200 | 0.021 | -0.003 | 0.048 | 0.960 |
| umap | M | 0.214 | pool_chemeleon | 0.065 | 0.025 | -0.008 | 0.060 | 0.936 |
| umap | M | 0.214 | pool_chemeleon | 0.100 | 0.039 | 0.012 | 0.071 | 0.999 |
| umap | M | 0.214 | pool_chemeleon | 0.200 | 0.056 | 0.024 | 0.095 | 1.000 |
| umap | M | 0.214 | pool_lgbm | 0.065 | 0.070 | 0.033 | 0.109 | 1.000 |
| umap | M | 0.214 | pool_lgbm | 0.100 | 0.071 | 0.039 | 0.102 | 1.000 |
| umap | M | 0.214 | pool_lgbm | 0.200 | 0.092 | 0.050 | 0.136 | 1.000 |
| umap | M | 0.214 | rank_mean_gates | 0.065 | 0.061 | 0.024 | 0.106 | 0.998 |
| umap | M | 0.214 | rank_mean_gates | 0.100 | 0.066 | 0.040 | 0.092 | 1.000 |
| umap | M | 0.214 | rank_mean_gates | 0.200 | 0.071 | 0.037 | 0.105 | 1.000 |
| umap | M | 0.214 | l4_rank_mean | 0.065 | 0.063 | 0.029 | 0.100 | 1.000 |
| umap | M | 0.214 | l4_rank_mean | 0.100 | 0.071 | 0.042 | 0.101 | 1.000 |
| umap | M | 0.214 | l4_rank_mean | 0.200 | 0.067 | 0.035 | 0.100 | 1.000 |
| umap | M | 0.231 | gate_x_tdi | 0.065 | 0.053 | 0.025 | 0.085 | 1.000 |
| umap | M | 0.231 | gate_x_tdi | 0.100 | 0.057 | 0.036 | 0.079 | 1.000 |
| umap | M | 0.231 | gate_x_tdi | 0.200 | 0.058 | 0.037 | 0.080 | 1.000 |
| umap | M | 0.231 | gate_x_tdi_chemeleon | 0.065 | 0.047 | -0.007 | 0.107 | 0.963 |
| umap | M | 0.231 | gate_x_tdi_chemeleon | 0.100 | 0.055 | 0.010 | 0.100 | 0.993 |
| umap | M | 0.231 | gate_x_tdi_chemeleon | 0.200 | 0.060 | 0.019 | 0.104 | 0.998 |
| umap | M | 0.231 | legacy_sep15_refit | 0.065 | 0.030 | -0.009 | 0.071 | 0.923 |
| umap | M | 0.231 | legacy_sep15_refit | 0.100 | 0.024 | -0.013 | 0.062 | 0.896 |
| umap | M | 0.231 | legacy_sep15_refit | 0.200 | 0.022 | -0.003 | 0.049 | 0.960 |
| umap | M | 0.231 | pool_chemeleon | 0.065 | 0.024 | -0.010 | 0.061 | 0.924 |
| umap | M | 0.231 | pool_chemeleon | 0.100 | 0.039 | 0.011 | 0.072 | 0.998 |
| umap | M | 0.231 | pool_chemeleon | 0.200 | 0.056 | 0.023 | 0.097 | 1.000 |
| umap | M | 0.231 | pool_lgbm | 0.065 | 0.069 | 0.031 | 0.110 | 1.000 |
| umap | M | 0.231 | pool_lgbm | 0.100 | 0.070 | 0.037 | 0.103 | 1.000 |
| umap | M | 0.231 | pool_lgbm | 0.200 | 0.093 | 0.050 | 0.137 | 1.000 |
| umap | M | 0.231 | rank_mean_gates | 0.065 | 0.061 | 0.023 | 0.107 | 0.998 |
| umap | M | 0.231 | rank_mean_gates | 0.100 | 0.065 | 0.039 | 0.092 | 1.000 |
| umap | M | 0.231 | rank_mean_gates | 0.200 | 0.071 | 0.036 | 0.107 | 0.999 |
| umap | M | 0.231 | l4_rank_mean | 0.065 | 0.064 | 0.030 | 0.103 | 1.000 |
| umap | M | 0.231 | l4_rank_mean | 0.100 | 0.072 | 0.043 | 0.103 | 1.000 |
| umap | M | 0.231 | l4_rank_mean | 0.200 | 0.068 | 0.035 | 0.102 | 1.000 |
| umap | M | 0.250 | gate_x_tdi | 0.065 | 0.053 | 0.024 | 0.085 | 1.000 |
| umap | M | 0.250 | gate_x_tdi | 0.100 | 0.056 | 0.034 | 0.079 | 1.000 |
| umap | M | 0.250 | gate_x_tdi | 0.200 | 0.057 | 0.037 | 0.080 | 1.000 |
| umap | M | 0.250 | gate_x_tdi_chemeleon | 0.065 | 0.046 | -0.008 | 0.108 | 0.962 |
| umap | M | 0.250 | gate_x_tdi_chemeleon | 0.100 | 0.054 | 0.009 | 0.101 | 0.989 |
| umap | M | 0.250 | gate_x_tdi_chemeleon | 0.200 | 0.060 | 0.018 | 0.105 | 0.997 |
| umap | M | 0.250 | legacy_sep15_refit | 0.065 | 0.032 | -0.009 | 0.074 | 0.926 |
| umap | M | 0.250 | legacy_sep15_refit | 0.100 | 0.024 | -0.014 | 0.065 | 0.897 |
| umap | M | 0.250 | legacy_sep15_refit | 0.200 | 0.023 | -0.003 | 0.051 | 0.959 |
| umap | M | 0.250 | pool_chemeleon | 0.065 | 0.023 | -0.011 | 0.061 | 0.909 |
| umap | M | 0.250 | pool_chemeleon | 0.100 | 0.038 | 0.009 | 0.072 | 0.995 |
| umap | M | 0.250 | pool_chemeleon | 0.200 | 0.056 | 0.021 | 0.098 | 0.999 |
| umap | M | 0.250 | pool_lgbm | 0.065 | 0.069 | 0.029 | 0.110 | 1.000 |
| umap | M | 0.250 | pool_lgbm | 0.100 | 0.069 | 0.035 | 0.103 | 1.000 |
| umap | M | 0.250 | pool_lgbm | 0.200 | 0.093 | 0.049 | 0.139 | 1.000 |
| umap | M | 0.250 | rank_mean_gates | 0.065 | 0.060 | 0.022 | 0.107 | 0.998 |
| umap | M | 0.250 | rank_mean_gates | 0.100 | 0.064 | 0.037 | 0.092 | 1.000 |
| umap | M | 0.250 | rank_mean_gates | 0.200 | 0.070 | 0.034 | 0.108 | 0.999 |
| umap | M | 0.250 | l4_rank_mean | 0.065 | 0.065 | 0.031 | 0.105 | 1.000 |
| umap | M | 0.250 | l4_rank_mean | 0.100 | 0.073 | 0.043 | 0.105 | 1.000 |
| umap | M | 0.250 | l4_rank_mean | 0.200 | 0.069 | 0.034 | 0.104 | 1.000 |
| umap | M | 0.300 | gate_x_tdi | 0.065 | 0.051 | 0.020 | 0.085 | 1.000 |
| umap | M | 0.300 | gate_x_tdi | 0.100 | 0.053 | 0.030 | 0.078 | 1.000 |
| umap | M | 0.300 | gate_x_tdi | 0.200 | 0.055 | 0.034 | 0.080 | 1.000 |
| umap | M | 0.300 | gate_x_tdi_chemeleon | 0.065 | 0.045 | -0.013 | 0.110 | 0.953 |
| umap | M | 0.300 | gate_x_tdi_chemeleon | 0.100 | 0.052 | 0.004 | 0.103 | 0.983 |
| umap | M | 0.300 | gate_x_tdi_chemeleon | 0.200 | 0.059 | 0.014 | 0.108 | 0.994 |
| umap | M | 0.300 | legacy_sep15_refit | 0.065 | 0.036 | -0.009 | 0.082 | 0.937 |
| umap | M | 0.300 | legacy_sep15_refit | 0.100 | 0.027 | -0.015 | 0.070 | 0.897 |
| umap | M | 0.300 | legacy_sep15_refit | 0.200 | 0.025 | -0.003 | 0.055 | 0.960 |
| umap | M | 0.300 | pool_chemeleon | 0.065 | 0.020 | -0.016 | 0.062 | 0.870 |
| umap | M | 0.300 | pool_chemeleon | 0.100 | 0.034 | 0.004 | 0.071 | 0.983 |
| umap | M | 0.300 | pool_chemeleon | 0.200 | 0.055 | 0.016 | 0.101 | 0.996 |
| umap | M | 0.300 | pool_lgbm | 0.065 | 0.066 | 0.024 | 0.111 | 0.999 |
| umap | M | 0.300 | pool_lgbm | 0.100 | 0.066 | 0.030 | 0.102 | 1.000 |
| umap | M | 0.300 | pool_lgbm | 0.200 | 0.092 | 0.046 | 0.141 | 1.000 |
| umap | M | 0.300 | rank_mean_gates | 0.065 | 0.058 | 0.016 | 0.108 | 0.995 |
| umap | M | 0.300 | rank_mean_gates | 0.100 | 0.061 | 0.032 | 0.090 | 1.000 |
| umap | M | 0.300 | rank_mean_gates | 0.200 | 0.069 | 0.029 | 0.109 | 0.999 |
| umap | M | 0.300 | l4_rank_mean | 0.065 | 0.067 | 0.031 | 0.110 | 0.999 |
| umap | M | 0.300 | l4_rank_mean | 0.100 | 0.075 | 0.043 | 0.109 | 1.000 |
| umap | M | 0.300 | l4_rank_mean | 0.200 | 0.071 | 0.034 | 0.109 | 1.000 |
| umap | M | 0.308 | gate_x_tdi | 0.065 | 0.050 | 0.019 | 0.085 | 1.000 |
| umap | M | 0.308 | gate_x_tdi | 0.100 | 0.052 | 0.029 | 0.078 | 1.000 |
| umap | M | 0.308 | gate_x_tdi | 0.200 | 0.055 | 0.033 | 0.080 | 1.000 |
| umap | M | 0.308 | gate_x_tdi_chemeleon | 0.065 | 0.044 | -0.014 | 0.110 | 0.951 |
| umap | M | 0.308 | gate_x_tdi_chemeleon | 0.100 | 0.052 | 0.003 | 0.104 | 0.981 |
| umap | M | 0.308 | gate_x_tdi_chemeleon | 0.200 | 0.058 | 0.014 | 0.108 | 0.993 |
| umap | M | 0.308 | legacy_sep15_refit | 0.065 | 0.036 | -0.009 | 0.083 | 0.939 |
| umap | M | 0.308 | legacy_sep15_refit | 0.100 | 0.027 | -0.015 | 0.070 | 0.898 |
| umap | M | 0.308 | legacy_sep15_refit | 0.200 | 0.025 | -0.003 | 0.056 | 0.961 |
| umap | M | 0.308 | pool_chemeleon | 0.065 | 0.019 | -0.017 | 0.061 | 0.858 |
| umap | M | 0.308 | pool_chemeleon | 0.100 | 0.034 | 0.002 | 0.071 | 0.981 |
| umap | M | 0.308 | pool_chemeleon | 0.200 | 0.055 | 0.016 | 0.102 | 0.996 |
| umap | M | 0.308 | pool_lgbm | 0.065 | 0.065 | 0.024 | 0.111 | 0.999 |
| umap | M | 0.308 | pool_lgbm | 0.100 | 0.065 | 0.029 | 0.102 | 1.000 |
| umap | M | 0.308 | pool_lgbm | 0.200 | 0.092 | 0.045 | 0.141 | 1.000 |
| umap | M | 0.308 | rank_mean_gates | 0.065 | 0.058 | 0.016 | 0.108 | 0.995 |
| umap | M | 0.308 | rank_mean_gates | 0.100 | 0.060 | 0.031 | 0.090 | 1.000 |
| umap | M | 0.308 | rank_mean_gates | 0.200 | 0.068 | 0.028 | 0.109 | 0.999 |
| umap | M | 0.308 | l4_rank_mean | 0.065 | 0.068 | 0.031 | 0.110 | 0.999 |
| umap | M | 0.308 | l4_rank_mean | 0.100 | 0.075 | 0.043 | 0.109 | 1.000 |
| umap | M | 0.308 | l4_rank_mean | 0.200 | 0.071 | 0.034 | 0.109 | 1.000 |
| umap | M | 0.311 | gate_x_tdi | 0.065 | 0.050 | 0.019 | 0.085 | 1.000 |
| umap | M | 0.311 | gate_x_tdi | 0.100 | 0.052 | 0.028 | 0.078 | 1.000 |
| umap | M | 0.311 | gate_x_tdi | 0.200 | 0.055 | 0.033 | 0.080 | 1.000 |
| umap | M | 0.311 | gate_x_tdi_chemeleon | 0.065 | 0.044 | -0.015 | 0.110 | 0.948 |
| umap | M | 0.311 | gate_x_tdi_chemeleon | 0.100 | 0.052 | 0.003 | 0.104 | 0.980 |
| umap | M | 0.311 | gate_x_tdi_chemeleon | 0.200 | 0.058 | 0.013 | 0.108 | 0.993 |
| umap | M | 0.311 | legacy_sep15_refit | 0.065 | 0.036 | -0.009 | 0.084 | 0.940 |
| umap | M | 0.311 | legacy_sep15_refit | 0.100 | 0.027 | -0.015 | 0.071 | 0.898 |
| umap | M | 0.311 | legacy_sep15_refit | 0.200 | 0.025 | -0.003 | 0.056 | 0.961 |
| umap | M | 0.311 | pool_chemeleon | 0.065 | 0.019 | -0.017 | 0.061 | 0.855 |
| umap | M | 0.311 | pool_chemeleon | 0.100 | 0.033 | 0.002 | 0.071 | 0.981 |
| umap | M | 0.311 | pool_chemeleon | 0.200 | 0.055 | 0.015 | 0.102 | 0.996 |
| umap | M | 0.311 | pool_lgbm | 0.065 | 0.065 | 0.023 | 0.111 | 0.999 |
| umap | M | 0.311 | pool_lgbm | 0.100 | 0.065 | 0.029 | 0.102 | 1.000 |
| umap | M | 0.311 | pool_lgbm | 0.200 | 0.092 | 0.045 | 0.141 | 1.000 |
| umap | M | 0.311 | rank_mean_gates | 0.065 | 0.057 | 0.015 | 0.108 | 0.995 |
| umap | M | 0.311 | rank_mean_gates | 0.100 | 0.060 | 0.030 | 0.090 | 1.000 |
| umap | M | 0.311 | rank_mean_gates | 0.200 | 0.068 | 0.028 | 0.109 | 0.999 |
| umap | M | 0.311 | l4_rank_mean | 0.065 | 0.068 | 0.031 | 0.111 | 0.999 |
| umap | M | 0.311 | l4_rank_mean | 0.100 | 0.075 | 0.043 | 0.109 | 1.000 |
| umap | M | 0.311 | l4_rank_mean | 0.200 | 0.071 | 0.033 | 0.110 | 1.000 |
| umap | M | 0.326 | gate_x_tdi | 0.065 | 0.049 | 0.017 | 0.084 | 0.999 |
| umap | M | 0.326 | gate_x_tdi | 0.100 | 0.050 | 0.027 | 0.077 | 1.000 |
| umap | M | 0.326 | gate_x_tdi | 0.200 | 0.054 | 0.032 | 0.079 | 1.000 |
| umap | M | 0.326 | gate_x_tdi_chemeleon | 0.065 | 0.044 | -0.016 | 0.111 | 0.937 |
| umap | M | 0.326 | gate_x_tdi_chemeleon | 0.100 | 0.051 | 0.002 | 0.104 | 0.977 |
| umap | M | 0.326 | gate_x_tdi_chemeleon | 0.200 | 0.058 | 0.012 | 0.108 | 0.991 |
| umap | M | 0.326 | legacy_sep15_refit | 0.065 | 0.037 | -0.008 | 0.086 | 0.941 |
| umap | M | 0.326 | legacy_sep15_refit | 0.100 | 0.028 | -0.015 | 0.073 | 0.898 |
| umap | M | 0.326 | legacy_sep15_refit | 0.200 | 0.026 | -0.003 | 0.057 | 0.961 |
| umap | M | 0.326 | pool_chemeleon | 0.065 | 0.018 | -0.019 | 0.061 | 0.833 |
| umap | M | 0.326 | pool_chemeleon | 0.100 | 0.032 | -0.000 | 0.070 | 0.974 |
| umap | M | 0.326 | pool_chemeleon | 0.200 | 0.054 | 0.014 | 0.102 | 0.996 |
| umap | M | 0.326 | pool_lgbm | 0.065 | 0.064 | 0.021 | 0.110 | 0.999 |
| umap | M | 0.326 | pool_lgbm | 0.100 | 0.064 | 0.027 | 0.101 | 1.000 |
| umap | M | 0.326 | pool_lgbm | 0.200 | 0.091 | 0.043 | 0.142 | 1.000 |
| umap | M | 0.326 | rank_mean_gates | 0.065 | 0.056 | 0.014 | 0.108 | 0.994 |
| umap | M | 0.326 | rank_mean_gates | 0.100 | 0.058 | 0.028 | 0.088 | 1.000 |
| umap | M | 0.326 | rank_mean_gates | 0.200 | 0.067 | 0.026 | 0.109 | 0.999 |
| umap | M | 0.326 | l4_rank_mean | 0.065 | 0.068 | 0.030 | 0.112 | 0.999 |
| umap | M | 0.326 | l4_rank_mean | 0.100 | 0.075 | 0.042 | 0.110 | 1.000 |
| umap | M | 0.326 | l4_rank_mean | 0.200 | 0.071 | 0.033 | 0.111 | 1.000 |
| umap | M | 0.400 | gate_x_tdi | 0.065 | 0.043 | 0.009 | 0.081 | 0.995 |
| umap | M | 0.400 | gate_x_tdi | 0.100 | 0.043 | 0.017 | 0.072 | 1.000 |
| umap | M | 0.400 | gate_x_tdi | 0.200 | 0.048 | 0.023 | 0.075 | 1.000 |
| umap | M | 0.400 | gate_x_tdi_chemeleon | 0.065 | 0.039 | -0.024 | 0.112 | 0.900 |
| umap | M | 0.400 | gate_x_tdi_chemeleon | 0.100 | 0.045 | -0.007 | 0.102 | 0.944 |
| umap | M | 0.400 | gate_x_tdi_chemeleon | 0.200 | 0.053 | 0.003 | 0.109 | 0.978 |
| umap | M | 0.400 | legacy_sep15_refit | 0.065 | 0.043 | -0.008 | 0.095 | 0.945 |
| umap | M | 0.400 | legacy_sep15_refit | 0.100 | 0.030 | -0.017 | 0.080 | 0.896 |
| umap | M | 0.400 | legacy_sep15_refit | 0.200 | 0.028 | -0.003 | 0.061 | 0.961 |
| umap | M | 0.400 | pool_chemeleon | 0.065 | 0.012 | -0.028 | 0.059 | 0.721 |
| umap | M | 0.400 | pool_chemeleon | 0.100 | 0.024 | -0.010 | 0.066 | 0.918 |
| umap | M | 0.400 | pool_chemeleon | 0.200 | 0.050 | 0.005 | 0.104 | 0.987 |
| umap | M | 0.400 | pool_lgbm | 0.065 | 0.057 | 0.012 | 0.107 | 0.994 |
| umap | M | 0.400 | pool_lgbm | 0.100 | 0.055 | 0.017 | 0.095 | 0.998 |
| umap | M | 0.400 | pool_lgbm | 0.200 | 0.086 | 0.034 | 0.140 | 0.999 |
| umap | M | 0.400 | rank_mean_gates | 0.065 | 0.051 | 0.005 | 0.106 | 0.983 |
| umap | M | 0.400 | rank_mean_gates | 0.100 | 0.050 | 0.017 | 0.083 | 1.000 |
| umap | M | 0.400 | rank_mean_gates | 0.200 | 0.061 | 0.016 | 0.107 | 0.995 |
| umap | M | 0.400 | l4_rank_mean | 0.065 | 0.070 | 0.028 | 0.115 | 0.999 |
| umap | M | 0.400 | l4_rank_mean | 0.100 | 0.075 | 0.039 | 0.114 | 1.000 |
| umap | M | 0.400 | l4_rank_mean | 0.200 | 0.071 | 0.030 | 0.113 | 1.000 |
| umap | M | 0.400 | gate_x_tdi | 0.065 | 0.043 | 0.009 | 0.081 | 0.995 |
| umap | M | 0.400 | gate_x_tdi | 0.100 | 0.043 | 0.017 | 0.072 | 1.000 |
| umap | M | 0.400 | gate_x_tdi | 0.200 | 0.048 | 0.023 | 0.075 | 1.000 |
| umap | M | 0.400 | gate_x_tdi_chemeleon | 0.065 | 0.039 | -0.024 | 0.112 | 0.900 |
| umap | M | 0.400 | gate_x_tdi_chemeleon | 0.100 | 0.045 | -0.007 | 0.102 | 0.943 |
| umap | M | 0.400 | gate_x_tdi_chemeleon | 0.200 | 0.053 | 0.003 | 0.109 | 0.978 |
| umap | M | 0.400 | legacy_sep15_refit | 0.065 | 0.043 | -0.008 | 0.095 | 0.945 |
| umap | M | 0.400 | legacy_sep15_refit | 0.100 | 0.030 | -0.017 | 0.080 | 0.896 |
| umap | M | 0.400 | legacy_sep15_refit | 0.200 | 0.028 | -0.003 | 0.061 | 0.961 |
| umap | M | 0.400 | pool_chemeleon | 0.065 | 0.012 | -0.028 | 0.059 | 0.721 |
| umap | M | 0.400 | pool_chemeleon | 0.100 | 0.024 | -0.010 | 0.066 | 0.918 |
| umap | M | 0.400 | pool_chemeleon | 0.200 | 0.050 | 0.005 | 0.104 | 0.987 |
| umap | M | 0.400 | pool_lgbm | 0.065 | 0.057 | 0.012 | 0.107 | 0.994 |
| umap | M | 0.400 | pool_lgbm | 0.100 | 0.055 | 0.017 | 0.095 | 0.998 |
| umap | M | 0.400 | pool_lgbm | 0.200 | 0.085 | 0.034 | 0.140 | 0.999 |
| umap | M | 0.400 | rank_mean_gates | 0.065 | 0.051 | 0.005 | 0.106 | 0.983 |
| umap | M | 0.400 | rank_mean_gates | 0.100 | 0.050 | 0.017 | 0.083 | 1.000 |
| umap | M | 0.400 | rank_mean_gates | 0.200 | 0.061 | 0.016 | 0.107 | 0.995 |
| umap | M | 0.400 | l4_rank_mean | 0.065 | 0.070 | 0.028 | 0.115 | 0.999 |
| umap | M | 0.400 | l4_rank_mean | 0.100 | 0.075 | 0.039 | 0.114 | 1.000 |
| umap | M | 0.400 | l4_rank_mean | 0.200 | 0.071 | 0.030 | 0.113 | 1.000 |
| umap | M | 0.403 | gate_x_tdi | 0.065 | 0.043 | 0.009 | 0.081 | 0.994 |
| umap | M | 0.403 | gate_x_tdi | 0.100 | 0.042 | 0.016 | 0.072 | 1.000 |
| umap | M | 0.403 | gate_x_tdi | 0.200 | 0.047 | 0.023 | 0.075 | 1.000 |
| umap | M | 0.403 | gate_x_tdi_chemeleon | 0.065 | 0.039 | -0.025 | 0.112 | 0.899 |
| umap | M | 0.403 | gate_x_tdi_chemeleon | 0.100 | 0.045 | -0.007 | 0.102 | 0.942 |
| umap | M | 0.403 | gate_x_tdi_chemeleon | 0.200 | 0.052 | 0.002 | 0.109 | 0.978 |
| umap | M | 0.403 | legacy_sep15_refit | 0.065 | 0.043 | -0.008 | 0.095 | 0.945 |
| umap | M | 0.403 | legacy_sep15_refit | 0.100 | 0.030 | -0.017 | 0.080 | 0.897 |
| umap | M | 0.403 | legacy_sep15_refit | 0.200 | 0.028 | -0.003 | 0.062 | 0.961 |
| umap | M | 0.403 | pool_chemeleon | 0.065 | 0.011 | -0.028 | 0.058 | 0.714 |
| umap | M | 0.403 | pool_chemeleon | 0.100 | 0.024 | -0.011 | 0.065 | 0.912 |
| umap | M | 0.403 | pool_chemeleon | 0.200 | 0.049 | 0.005 | 0.104 | 0.987 |
| umap | M | 0.403 | pool_lgbm | 0.065 | 0.057 | 0.012 | 0.107 | 0.994 |
| umap | M | 0.403 | pool_lgbm | 0.100 | 0.055 | 0.017 | 0.095 | 0.998 |
| umap | M | 0.403 | pool_lgbm | 0.200 | 0.085 | 0.034 | 0.140 | 0.999 |
| umap | M | 0.403 | rank_mean_gates | 0.065 | 0.050 | 0.005 | 0.106 | 0.982 |
| umap | M | 0.403 | rank_mean_gates | 0.100 | 0.049 | 0.016 | 0.083 | 1.000 |
| umap | M | 0.403 | rank_mean_gates | 0.200 | 0.061 | 0.015 | 0.106 | 0.994 |
| umap | M | 0.403 | l4_rank_mean | 0.065 | 0.070 | 0.028 | 0.115 | 0.999 |
| umap | M | 0.403 | l4_rank_mean | 0.100 | 0.075 | 0.039 | 0.115 | 1.000 |
| umap | M | 0.403 | l4_rank_mean | 0.200 | 0.071 | 0.029 | 0.113 | 1.000 |
| umap | M | 0.416 | gate_x_tdi | 0.065 | 0.042 | 0.008 | 0.080 | 0.993 |
| umap | M | 0.416 | gate_x_tdi | 0.100 | 0.041 | 0.015 | 0.071 | 0.998 |
| umap | M | 0.416 | gate_x_tdi | 0.200 | 0.046 | 0.021 | 0.074 | 1.000 |
| umap | M | 0.416 | gate_x_tdi_chemeleon | 0.065 | 0.038 | -0.026 | 0.112 | 0.886 |
| umap | M | 0.416 | gate_x_tdi_chemeleon | 0.100 | 0.044 | -0.009 | 0.102 | 0.936 |
| umap | M | 0.416 | gate_x_tdi_chemeleon | 0.200 | 0.051 | 0.001 | 0.109 | 0.977 |
| umap | M | 0.416 | legacy_sep15_refit | 0.065 | 0.044 | -0.007 | 0.096 | 0.946 |
| umap | M | 0.416 | legacy_sep15_refit | 0.100 | 0.031 | -0.017 | 0.081 | 0.895 |
| umap | M | 0.416 | legacy_sep15_refit | 0.200 | 0.029 | -0.003 | 0.062 | 0.961 |
| umap | M | 0.416 | pool_chemeleon | 0.065 | 0.010 | -0.030 | 0.058 | 0.689 |
| umap | M | 0.416 | pool_chemeleon | 0.100 | 0.022 | -0.013 | 0.064 | 0.890 |
| umap | M | 0.416 | pool_chemeleon | 0.200 | 0.048 | 0.003 | 0.104 | 0.982 |
| umap | M | 0.416 | pool_lgbm | 0.065 | 0.056 | 0.010 | 0.106 | 0.993 |
| umap | M | 0.416 | pool_lgbm | 0.100 | 0.053 | 0.014 | 0.094 | 0.998 |
| umap | M | 0.416 | pool_lgbm | 0.200 | 0.084 | 0.032 | 0.139 | 0.999 |
| umap | M | 0.416 | rank_mean_gates | 0.065 | 0.049 | 0.003 | 0.106 | 0.980 |
| umap | M | 0.416 | rank_mean_gates | 0.100 | 0.048 | 0.014 | 0.082 | 1.000 |
| umap | M | 0.416 | rank_mean_gates | 0.200 | 0.059 | 0.013 | 0.106 | 0.993 |
| umap | M | 0.416 | l4_rank_mean | 0.065 | 0.070 | 0.028 | 0.115 | 0.999 |
| umap | M | 0.416 | l4_rank_mean | 0.100 | 0.075 | 0.038 | 0.115 | 1.000 |
| umap | M | 0.416 | l4_rank_mean | 0.200 | 0.071 | 0.029 | 0.114 | 1.000 |
| umap | M | 0.500 | gate_x_tdi | 0.065 | 0.034 | -0.003 | 0.073 | 0.959 |
| umap | M | 0.500 | gate_x_tdi | 0.100 | 0.029 | 0.002 | 0.062 | 0.980 |
| umap | M | 0.500 | gate_x_tdi | 0.200 | 0.036 | 0.009 | 0.066 | 0.997 |
| umap | M | 0.500 | gate_x_tdi_chemeleon | 0.065 | 0.031 | -0.036 | 0.110 | 0.797 |
| umap | M | 0.500 | gate_x_tdi_chemeleon | 0.100 | 0.035 | -0.023 | 0.097 | 0.874 |
| umap | M | 0.500 | gate_x_tdi_chemeleon | 0.200 | 0.043 | -0.012 | 0.103 | 0.934 |
| umap | M | 0.500 | legacy_sep15_refit | 0.065 | 0.050 | -0.006 | 0.107 | 0.953 |
| umap | M | 0.500 | legacy_sep15_refit | 0.100 | 0.033 | -0.021 | 0.088 | 0.891 |
| umap | M | 0.500 | legacy_sep15_refit | 0.200 | 0.031 | -0.003 | 0.067 | 0.961 |
| umap | M | 0.500 | pool_chemeleon | 0.065 | 0.002 | -0.041 | 0.053 | 0.555 |
| umap | M | 0.500 | pool_chemeleon | 0.100 | 0.011 | -0.026 | 0.055 | 0.727 |
| umap | M | 0.500 | pool_chemeleon | 0.200 | 0.040 | -0.008 | 0.099 | 0.948 |
| umap | M | 0.500 | pool_lgbm | 0.065 | 0.046 | -0.002 | 0.100 | 0.969 |
| umap | M | 0.500 | pool_lgbm | 0.100 | 0.041 | 0.000 | 0.084 | 0.975 |
| umap | M | 0.500 | pool_lgbm | 0.200 | 0.074 | 0.019 | 0.132 | 0.993 |
| umap | M | 0.500 | rank_mean_gates | 0.065 | 0.041 | -0.008 | 0.101 | 0.956 |
| umap | M | 0.500 | rank_mean_gates | 0.100 | 0.035 | -0.001 | 0.072 | 0.973 |
| umap | M | 0.500 | rank_mean_gates | 0.200 | 0.049 | -0.001 | 0.100 | 0.973 |
| umap | M | 0.500 | l4_rank_mean | 0.065 | 0.069 | 0.025 | 0.117 | 0.999 |
| umap | M | 0.500 | l4_rank_mean | 0.100 | 0.072 | 0.034 | 0.115 | 1.000 |
| umap | M | 0.500 | l4_rank_mean | 0.200 | 0.068 | 0.023 | 0.114 | 0.999 |
| umap | M | 0.600 | gate_x_tdi | 0.065 | 0.022 | -0.016 | 0.065 | 0.876 |
| umap | M | 0.600 | gate_x_tdi | 0.100 | 0.013 | -0.016 | 0.050 | 0.785 |
| umap | M | 0.600 | gate_x_tdi | 0.200 | 0.021 | -0.008 | 0.054 | 0.911 |
| umap | M | 0.600 | gate_x_tdi_chemeleon | 0.065 | 0.020 | -0.049 | 0.106 | 0.691 |
| umap | M | 0.600 | gate_x_tdi_chemeleon | 0.100 | 0.021 | -0.040 | 0.087 | 0.722 |
| umap | M | 0.600 | gate_x_tdi_chemeleon | 0.200 | 0.029 | -0.028 | 0.097 | 0.843 |
| umap | M | 0.600 | legacy_sep15_refit | 0.065 | 0.057 | -0.004 | 0.118 | 0.960 |
| umap | M | 0.600 | legacy_sep15_refit | 0.100 | 0.035 | -0.023 | 0.092 | 0.888 |
| umap | M | 0.600 | legacy_sep15_refit | 0.200 | 0.033 | -0.003 | 0.071 | 0.961 |
| umap | M | 0.600 | pool_chemeleon | 0.065 | -0.010 | -0.054 | 0.046 | 0.371 |
| umap | M | 0.600 | pool_chemeleon | 0.100 | -0.004 | -0.044 | 0.042 | 0.454 |
| umap | M | 0.600 | pool_chemeleon | 0.200 | 0.027 | -0.023 | 0.091 | 0.860 |
| umap | M | 0.600 | pool_lgbm | 0.065 | 0.033 | -0.017 | 0.091 | 0.888 |
| umap | M | 0.600 | pool_lgbm | 0.100 | 0.023 | -0.018 | 0.069 | 0.864 |
| umap | M | 0.600 | pool_lgbm | 0.200 | 0.058 | 0.001 | 0.120 | 0.978 |
| umap | M | 0.600 | rank_mean_gates | 0.065 | 0.029 | -0.022 | 0.094 | 0.847 |
| umap | M | 0.600 | rank_mean_gates | 0.100 | 0.018 | -0.020 | 0.058 | 0.818 |
| umap | M | 0.600 | rank_mean_gates | 0.200 | 0.033 | -0.019 | 0.088 | 0.894 |
| umap | M | 0.600 | l4_rank_mean | 0.065 | 0.067 | 0.020 | 0.117 | 0.996 |
| umap | M | 0.600 | l4_rank_mean | 0.100 | 0.066 | 0.026 | 0.113 | 1.000 |
| umap | M | 0.600 | l4_rank_mean | 0.200 | 0.062 | 0.014 | 0.111 | 0.997 |
| umap | M | 0.700 | gate_x_tdi | 0.065 | 0.009 | -0.032 | 0.056 | 0.691 |
| umap | M | 0.700 | gate_x_tdi | 0.100 | -0.004 | -0.038 | 0.036 | 0.407 |
| umap | M | 0.700 | gate_x_tdi | 0.200 | 0.003 | -0.027 | 0.039 | 0.594 |
| umap | M | 0.700 | gate_x_tdi_chemeleon | 0.065 | 0.009 | -0.063 | 0.101 | 0.569 |
| umap | M | 0.700 | gate_x_tdi_chemeleon | 0.100 | 0.006 | -0.059 | 0.075 | 0.552 |
| umap | M | 0.700 | gate_x_tdi_chemeleon | 0.200 | 0.014 | -0.048 | 0.086 | 0.673 |
| umap | M | 0.700 | legacy_sep15_refit | 0.065 | 0.063 | -0.003 | 0.131 | 0.966 |
| umap | M | 0.700 | legacy_sep15_refit | 0.100 | 0.037 | -0.026 | 0.099 | 0.886 |
| umap | M | 0.700 | legacy_sep15_refit | 0.200 | 0.035 | -0.003 | 0.076 | 0.960 |
| umap | M | 0.700 | pool_chemeleon | 0.065 | -0.023 | -0.070 | 0.035 | 0.225 |
| umap | M | 0.700 | pool_chemeleon | 0.100 | -0.022 | -0.064 | 0.029 | 0.190 |
| umap | M | 0.700 | pool_chemeleon | 0.200 | 0.012 | -0.041 | 0.079 | 0.682 |
| umap | M | 0.700 | pool_lgbm | 0.065 | 0.019 | -0.033 | 0.081 | 0.733 |
| umap | M | 0.700 | pool_lgbm | 0.100 | 0.004 | -0.039 | 0.053 | 0.565 |
| umap | M | 0.700 | pool_lgbm | 0.200 | 0.040 | -0.019 | 0.104 | 0.906 |
| umap | M | 0.700 | rank_mean_gates | 0.065 | 0.016 | -0.037 | 0.085 | 0.692 |
| umap | M | 0.700 | rank_mean_gates | 0.100 | -0.001 | -0.042 | 0.041 | 0.468 |
| umap | M | 0.700 | rank_mean_gates | 0.200 | 0.015 | -0.040 | 0.072 | 0.713 |
| umap | M | 0.700 | l4_rank_mean | 0.065 | 0.063 | 0.013 | 0.115 | 0.994 |
| umap | M | 0.700 | l4_rank_mean | 0.100 | 0.058 | 0.016 | 0.107 | 1.000 |
| umap | M | 0.700 | l4_rank_mean | 0.200 | 0.054 | 0.003 | 0.105 | 0.986 |
| umap | M | 0.800 | gate_x_tdi | 0.065 | -0.005 | -0.047 | 0.044 | 0.427 |
| umap | M | 0.800 | gate_x_tdi | 0.100 | -0.024 | -0.062 | 0.017 | 0.124 |
| umap | M | 0.800 | gate_x_tdi | 0.200 | -0.017 | -0.048 | 0.020 | 0.177 |
| umap | M | 0.800 | gate_x_tdi_chemeleon | 0.065 | -0.003 | -0.079 | 0.094 | 0.456 |
| umap | M | 0.800 | gate_x_tdi_chemeleon | 0.100 | -0.011 | -0.079 | 0.061 | 0.380 |
| umap | M | 0.800 | gate_x_tdi_chemeleon | 0.200 | -0.005 | -0.069 | 0.072 | 0.454 |
| umap | M | 0.800 | legacy_sep15_refit | 0.065 | 0.070 | 0.000 | 0.142 | 0.975 |
| umap | M | 0.800 | legacy_sep15_refit | 0.100 | 0.039 | -0.030 | 0.104 | 0.874 |
| umap | M | 0.800 | legacy_sep15_refit | 0.200 | 0.037 | -0.003 | 0.080 | 0.958 |
| umap | M | 0.800 | pool_chemeleon | 0.065 | -0.036 | -0.086 | 0.023 | 0.117 |
| umap | M | 0.800 | pool_chemeleon | 0.100 | -0.041 | -0.086 | 0.013 | 0.066 |
| umap | M | 0.800 | pool_chemeleon | 0.200 | -0.006 | -0.060 | 0.064 | 0.470 |
| umap | M | 0.800 | pool_lgbm | 0.065 | 0.003 | -0.052 | 0.070 | 0.551 |
| umap | M | 0.800 | pool_lgbm | 0.100 | -0.017 | -0.061 | 0.033 | 0.248 |
| umap | M | 0.800 | pool_lgbm | 0.200 | 0.018 | -0.041 | 0.085 | 0.720 |
| umap | M | 0.800 | rank_mean_gates | 0.065 | 0.002 | -0.053 | 0.074 | 0.505 |
| umap | M | 0.800 | rank_mean_gates | 0.100 | -0.022 | -0.067 | 0.024 | 0.167 |
| umap | M | 0.800 | rank_mean_gates | 0.200 | -0.006 | -0.063 | 0.054 | 0.460 |
| umap | M | 0.800 | l4_rank_mean | 0.065 | 0.057 | 0.007 | 0.112 | 0.985 |
| umap | M | 0.800 | l4_rank_mean | 0.100 | 0.048 | 0.004 | 0.100 | 0.988 |
| umap | M | 0.800 | l4_rank_mean | 0.200 | 0.044 | -0.010 | 0.097 | 0.945 |
| umap | M | 0.900 | gate_x_tdi | 0.065 | -0.020 | -0.066 | 0.031 | 0.217 |
| umap | M | 0.900 | gate_x_tdi | 0.100 | -0.045 | -0.087 | -0.001 | 0.022 |
| umap | M | 0.900 | gate_x_tdi | 0.200 | -0.039 | -0.071 | -0.001 | 0.023 |
| umap | M | 0.900 | gate_x_tdi_chemeleon | 0.065 | -0.017 | -0.096 | 0.081 | 0.339 |
| umap | M | 0.900 | gate_x_tdi_chemeleon | 0.100 | -0.030 | -0.103 | 0.045 | 0.215 |
| umap | M | 0.900 | gate_x_tdi_chemeleon | 0.200 | -0.025 | -0.091 | 0.055 | 0.250 |
| umap | M | 0.900 | legacy_sep15_refit | 0.065 | 0.078 | 0.003 | 0.154 | 0.980 |
| umap | M | 0.900 | legacy_sep15_refit | 0.100 | 0.040 | -0.033 | 0.111 | 0.861 |
| umap | M | 0.900 | legacy_sep15_refit | 0.200 | 0.038 | -0.005 | 0.084 | 0.955 |
| umap | M | 0.900 | pool_chemeleon | 0.065 | -0.051 | -0.103 | 0.010 | 0.063 |
| umap | M | 0.900 | pool_chemeleon | 0.100 | -0.061 | -0.110 | -0.006 | 0.019 |
| umap | M | 0.900 | pool_chemeleon | 0.200 | -0.025 | -0.083 | 0.048 | 0.237 |
| umap | M | 0.900 | pool_lgbm | 0.065 | -0.014 | -0.072 | 0.055 | 0.351 |
| umap | M | 0.900 | pool_lgbm | 0.100 | -0.040 | -0.086 | 0.013 | 0.066 |
| umap | M | 0.900 | pool_lgbm | 0.200 | -0.006 | -0.066 | 0.062 | 0.445 |
| umap | M | 0.900 | rank_mean_gates | 0.065 | -0.013 | -0.073 | 0.062 | 0.347 |
| umap | M | 0.900 | rank_mean_gates | 0.100 | -0.045 | -0.093 | 0.002 | 0.037 |
| umap | M | 0.900 | rank_mean_gates | 0.200 | -0.030 | -0.089 | 0.034 | 0.190 |
| umap | M | 0.900 | l4_rank_mean | 0.065 | 0.051 | -0.000 | 0.108 | 0.974 |
| umap | M | 0.900 | l4_rank_mean | 0.100 | 0.037 | -0.008 | 0.090 | 0.946 |
| umap | M | 0.900 | l4_rank_mean | 0.200 | 0.032 | -0.022 | 0.089 | 0.857 |
| umap | M | 1.000 | gate_x_tdi | 0.065 | -0.036 | -0.087 | 0.017 | 0.099 |
| umap | M | 1.000 | gate_x_tdi | 0.100 | -0.068 | -0.115 | -0.021 | 0.001 |
| umap | M | 1.000 | gate_x_tdi | 0.200 | -0.064 | -0.097 | -0.025 | 0.002 |
| umap | M | 1.000 | gate_x_tdi_chemeleon | 0.065 | -0.031 | -0.113 | 0.071 | 0.263 |
| umap | M | 1.000 | gate_x_tdi_chemeleon | 0.100 | -0.050 | -0.127 | 0.026 | 0.106 |
| umap | M | 1.000 | gate_x_tdi_chemeleon | 0.200 | -0.047 | -0.116 | 0.035 | 0.113 |
| umap | M | 1.000 | legacy_sep15_refit | 0.065 | 0.085 | 0.005 | 0.168 | 0.984 |
| umap | M | 1.000 | legacy_sep15_refit | 0.100 | 0.041 | -0.037 | 0.115 | 0.851 |
| umap | M | 1.000 | legacy_sep15_refit | 0.200 | 0.039 | -0.007 | 0.088 | 0.946 |
| umap | M | 1.000 | pool_chemeleon | 0.065 | -0.066 | -0.123 | -0.002 | 0.021 |
| umap | M | 1.000 | pool_chemeleon | 0.100 | -0.083 | -0.137 | -0.024 | 0.003 |
| umap | M | 1.000 | pool_chemeleon | 0.200 | -0.047 | -0.104 | 0.027 | 0.112 |
| umap | M | 1.000 | pool_lgbm | 0.065 | -0.031 | -0.093 | 0.039 | 0.175 |
| umap | M | 1.000 | pool_lgbm | 0.100 | -0.065 | -0.112 | -0.010 | 0.013 |
| umap | M | 1.000 | pool_lgbm | 0.200 | -0.032 | -0.092 | 0.037 | 0.176 |
| umap | M | 1.000 | rank_mean_gates | 0.065 | -0.029 | -0.093 | 0.050 | 0.215 |
| umap | M | 1.000 | rank_mean_gates | 0.100 | -0.069 | -0.124 | -0.019 | 0.007 |
| umap | M | 1.000 | rank_mean_gates | 0.200 | -0.055 | -0.116 | 0.011 | 0.050 |
| umap | M | 1.000 | l4_rank_mean | 0.065 | 0.044 | -0.010 | 0.103 | 0.948 |
| umap | M | 1.000 | l4_rank_mean | 0.100 | 0.024 | -0.021 | 0.079 | 0.832 |
| umap | M | 1.000 | l4_rank_mean | 0.200 | 0.018 | -0.037 | 0.078 | 0.735 |
| scaffold | O | 0.050 | gate_x_tdi | 0.065 | -0.006 | -0.040 | 0.022 | 0.360 |
| scaffold | O | 0.050 | gate_x_tdi | 0.100 | -0.010 | -0.040 | 0.021 | 0.257 |
| scaffold | O | 0.050 | gate_x_tdi | 0.200 | 0.000 | -0.026 | 0.027 | 0.536 |
| scaffold | O | 0.050 | gate_x_tdi_chemeleon | 0.065 | -0.012 | -0.054 | 0.030 | 0.272 |
| scaffold | O | 0.050 | gate_x_tdi_chemeleon | 0.100 | -0.012 | -0.051 | 0.026 | 0.266 |
| scaffold | O | 0.050 | gate_x_tdi_chemeleon | 0.200 | 0.001 | -0.035 | 0.040 | 0.528 |
| scaffold | O | 0.050 | legacy_sep15_refit | 0.065 | 0.041 | -0.000 | 0.084 | 0.973 |
| scaffold | O | 0.050 | legacy_sep15_refit | 0.100 | 0.043 | 0.007 | 0.082 | 0.991 |
| scaffold | O | 0.050 | legacy_sep15_refit | 0.200 | 0.053 | 0.020 | 0.087 | 0.998 |
| scaffold | O | 0.050 | pool_chemeleon | 0.065 | -0.018 | -0.060 | 0.024 | 0.220 |
| scaffold | O | 0.050 | pool_chemeleon | 0.100 | -0.021 | -0.062 | 0.018 | 0.152 |
| scaffold | O | 0.050 | pool_chemeleon | 0.200 | -0.001 | -0.038 | 0.039 | 0.483 |
| scaffold | O | 0.050 | pool_lgbm | 0.065 | 0.019 | -0.019 | 0.057 | 0.831 |
| scaffold | O | 0.050 | pool_lgbm | 0.100 | 0.007 | -0.028 | 0.043 | 0.661 |
| scaffold | O | 0.050 | pool_lgbm | 0.200 | -0.004 | -0.034 | 0.030 | 0.436 |
| scaffold | O | 0.050 | rank_mean_gates | 0.065 | -0.005 | -0.042 | 0.031 | 0.415 |
| scaffold | O | 0.050 | rank_mean_gates | 0.100 | -0.008 | -0.041 | 0.026 | 0.343 |
| scaffold | O | 0.050 | rank_mean_gates | 0.200 | 0.015 | -0.015 | 0.048 | 0.835 |
| scaffold | O | 0.050 | l4_rank_mean | 0.065 | 0.047 | 0.012 | 0.083 | 0.995 |
| scaffold | O | 0.050 | l4_rank_mean | 0.100 | 0.038 | 0.006 | 0.070 | 0.991 |
| scaffold | O | 0.050 | l4_rank_mean | 0.200 | 0.031 | 0.001 | 0.062 | 0.979 |
| scaffold | O | 0.071 | gate_x_tdi | 0.065 | -0.007 | -0.046 | 0.026 | 0.362 |
| scaffold | O | 0.071 | gate_x_tdi | 0.100 | -0.012 | -0.045 | 0.025 | 0.262 |
| scaffold | O | 0.071 | gate_x_tdi | 0.200 | 0.000 | -0.030 | 0.032 | 0.541 |
| scaffold | O | 0.071 | gate_x_tdi_chemeleon | 0.065 | -0.014 | -0.062 | 0.035 | 0.273 |
| scaffold | O | 0.071 | gate_x_tdi_chemeleon | 0.100 | -0.014 | -0.059 | 0.030 | 0.268 |
| scaffold | O | 0.071 | gate_x_tdi_chemeleon | 0.200 | 0.001 | -0.041 | 0.046 | 0.532 |
| scaffold | O | 0.071 | legacy_sep15_refit | 0.065 | 0.047 | -0.000 | 0.096 | 0.974 |
| scaffold | O | 0.071 | legacy_sep15_refit | 0.100 | 0.049 | 0.008 | 0.093 | 0.991 |
| scaffold | O | 0.071 | legacy_sep15_refit | 0.200 | 0.061 | 0.023 | 0.101 | 0.998 |
| scaffold | O | 0.071 | pool_chemeleon | 0.065 | -0.021 | -0.070 | 0.027 | 0.223 |
| scaffold | O | 0.071 | pool_chemeleon | 0.100 | -0.024 | -0.071 | 0.021 | 0.152 |
| scaffold | O | 0.071 | pool_chemeleon | 0.200 | -0.001 | -0.044 | 0.045 | 0.489 |
| scaffold | O | 0.071 | pool_lgbm | 0.065 | 0.022 | -0.021 | 0.066 | 0.834 |
| scaffold | O | 0.071 | pool_lgbm | 0.100 | 0.008 | -0.032 | 0.050 | 0.665 |
| scaffold | O | 0.071 | pool_lgbm | 0.200 | -0.004 | -0.040 | 0.035 | 0.442 |
| scaffold | O | 0.071 | rank_mean_gates | 0.065 | -0.006 | -0.048 | 0.037 | 0.417 |
| scaffold | O | 0.071 | rank_mean_gates | 0.100 | -0.010 | -0.047 | 0.030 | 0.345 |
| scaffold | O | 0.071 | rank_mean_gates | 0.200 | 0.018 | -0.017 | 0.056 | 0.835 |
| scaffold | O | 0.071 | l4_rank_mean | 0.065 | 0.054 | 0.014 | 0.095 | 0.995 |
| scaffold | O | 0.071 | l4_rank_mean | 0.100 | 0.044 | 0.007 | 0.081 | 0.991 |
| scaffold | O | 0.071 | l4_rank_mean | 0.200 | 0.036 | 0.001 | 0.072 | 0.979 |
| scaffold | O | 0.090 | gate_x_tdi | 0.065 | -0.008 | -0.050 | 0.029 | 0.367 |
| scaffold | O | 0.090 | gate_x_tdi | 0.100 | -0.013 | -0.050 | 0.027 | 0.262 |
| scaffold | O | 0.090 | gate_x_tdi | 0.200 | 0.000 | -0.034 | 0.036 | 0.546 |
| scaffold | O | 0.090 | gate_x_tdi_chemeleon | 0.065 | -0.015 | -0.068 | 0.039 | 0.274 |
| scaffold | O | 0.090 | gate_x_tdi_chemeleon | 0.100 | -0.016 | -0.064 | 0.034 | 0.270 |
| scaffold | O | 0.090 | gate_x_tdi_chemeleon | 0.200 | 0.001 | -0.045 | 0.052 | 0.536 |
| scaffold | O | 0.090 | legacy_sep15_refit | 0.065 | 0.051 | -0.000 | 0.105 | 0.974 |
| scaffold | O | 0.090 | legacy_sep15_refit | 0.100 | 0.054 | 0.009 | 0.101 | 0.991 |
| scaffold | O | 0.090 | legacy_sep15_refit | 0.200 | 0.067 | 0.025 | 0.111 | 0.998 |
| scaffold | O | 0.090 | pool_chemeleon | 0.065 | -0.023 | -0.076 | 0.030 | 0.224 |
| scaffold | O | 0.090 | pool_chemeleon | 0.100 | -0.027 | -0.079 | 0.024 | 0.155 |
| scaffold | O | 0.090 | pool_chemeleon | 0.200 | -0.000 | -0.048 | 0.050 | 0.492 |
| scaffold | O | 0.090 | pool_lgbm | 0.065 | 0.025 | -0.023 | 0.072 | 0.835 |
| scaffold | O | 0.090 | pool_lgbm | 0.100 | 0.009 | -0.034 | 0.055 | 0.670 |
| scaffold | O | 0.090 | pool_lgbm | 0.200 | -0.004 | -0.044 | 0.039 | 0.447 |
| scaffold | O | 0.090 | rank_mean_gates | 0.065 | -0.006 | -0.053 | 0.041 | 0.420 |
| scaffold | O | 0.090 | rank_mean_gates | 0.100 | -0.011 | -0.052 | 0.033 | 0.347 |
| scaffold | O | 0.090 | rank_mean_gates | 0.200 | 0.020 | -0.018 | 0.063 | 0.837 |
| scaffold | O | 0.090 | l4_rank_mean | 0.065 | 0.059 | 0.016 | 0.104 | 0.995 |
| scaffold | O | 0.090 | l4_rank_mean | 0.100 | 0.048 | 0.007 | 0.089 | 0.991 |
| scaffold | O | 0.090 | l4_rank_mean | 0.200 | 0.040 | 0.002 | 0.080 | 0.980 |
| scaffold | M | 0.000 | gate_x_tdi | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | gate_x_tdi | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | gate_x_tdi | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | gate_x_tdi_chemeleon | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | gate_x_tdi_chemeleon | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | gate_x_tdi_chemeleon | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | legacy_sep15_refit | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | legacy_sep15_refit | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | legacy_sep15_refit | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | pool_chemeleon | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | pool_chemeleon | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | pool_chemeleon | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | pool_lgbm | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | pool_lgbm | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | pool_lgbm | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | rank_mean_gates | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | rank_mean_gates | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | rank_mean_gates | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | l4_rank_mean | 0.065 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | l4_rank_mean | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.000 | l4_rank_mean | 0.200 | 0.000 | 0.000 | 0.000 | 0.000 |
| scaffold | M | 0.100 | gate_x_tdi | 0.065 | 0.071 | 0.047 | 0.100 | 1.000 |
| scaffold | M | 0.100 | gate_x_tdi | 0.100 | 0.066 | 0.042 | 0.091 | 1.000 |
| scaffold | M | 0.100 | gate_x_tdi | 0.200 | 0.083 | 0.063 | 0.104 | 1.000 |
| scaffold | M | 0.100 | gate_x_tdi_chemeleon | 0.065 | 0.066 | 0.034 | 0.101 | 1.000 |
| scaffold | M | 0.100 | gate_x_tdi_chemeleon | 0.100 | 0.066 | 0.037 | 0.095 | 1.000 |
| scaffold | M | 0.100 | gate_x_tdi_chemeleon | 0.200 | 0.081 | 0.057 | 0.105 | 1.000 |
| scaffold | M | 0.100 | legacy_sep15_refit | 0.065 | 0.025 | 0.001 | 0.049 | 0.977 |
| scaffold | M | 0.100 | legacy_sep15_refit | 0.100 | 0.022 | 0.001 | 0.042 | 0.976 |
| scaffold | M | 0.100 | legacy_sep15_refit | 0.200 | 0.024 | 0.006 | 0.041 | 0.997 |
| scaffold | M | 0.100 | pool_chemeleon | 0.065 | 0.050 | 0.019 | 0.084 | 1.000 |
| scaffold | M | 0.100 | pool_chemeleon | 0.100 | 0.048 | 0.020 | 0.077 | 0.999 |
| scaffold | M | 0.100 | pool_chemeleon | 0.200 | 0.077 | 0.051 | 0.103 | 1.000 |
| scaffold | M | 0.100 | pool_lgbm | 0.065 | 0.095 | 0.064 | 0.130 | 1.000 |
| scaffold | M | 0.100 | pool_lgbm | 0.100 | 0.092 | 0.061 | 0.123 | 1.000 |
| scaffold | M | 0.100 | pool_lgbm | 0.200 | 0.095 | 0.069 | 0.122 | 1.000 |
| scaffold | M | 0.100 | rank_mean_gates | 0.065 | 0.072 | 0.043 | 0.103 | 1.000 |
| scaffold | M | 0.100 | rank_mean_gates | 0.100 | 0.075 | 0.049 | 0.103 | 1.000 |
| scaffold | M | 0.100 | rank_mean_gates | 0.200 | 0.102 | 0.079 | 0.126 | 1.000 |
| scaffold | M | 0.100 | l4_rank_mean | 0.065 | 0.069 | 0.044 | 0.096 | 1.000 |
| scaffold | M | 0.100 | l4_rank_mean | 0.100 | 0.063 | 0.040 | 0.085 | 1.000 |
| scaffold | M | 0.100 | l4_rank_mean | 0.200 | 0.062 | 0.043 | 0.080 | 1.000 |
| scaffold | M | 0.200 | gate_x_tdi | 0.065 | 0.081 | 0.050 | 0.118 | 1.000 |
| scaffold | M | 0.200 | gate_x_tdi | 0.100 | 0.077 | 0.047 | 0.110 | 1.000 |
| scaffold | M | 0.200 | gate_x_tdi | 0.200 | 0.101 | 0.074 | 0.128 | 1.000 |
| scaffold | M | 0.200 | gate_x_tdi_chemeleon | 0.065 | 0.075 | 0.034 | 0.120 | 1.000 |
| scaffold | M | 0.200 | gate_x_tdi_chemeleon | 0.100 | 0.077 | 0.039 | 0.114 | 1.000 |
| scaffold | M | 0.200 | gate_x_tdi_chemeleon | 0.200 | 0.099 | 0.066 | 0.132 | 1.000 |
| scaffold | M | 0.200 | legacy_sep15_refit | 0.065 | 0.035 | 0.002 | 0.067 | 0.977 |
| scaffold | M | 0.200 | legacy_sep15_refit | 0.100 | 0.031 | 0.001 | 0.060 | 0.977 |
| scaffold | M | 0.200 | legacy_sep15_refit | 0.200 | 0.034 | 0.008 | 0.058 | 0.998 |
| scaffold | M | 0.200 | pool_chemeleon | 0.065 | 0.056 | 0.015 | 0.100 | 0.995 |
| scaffold | M | 0.200 | pool_chemeleon | 0.100 | 0.055 | 0.017 | 0.094 | 0.999 |
| scaffold | M | 0.200 | pool_chemeleon | 0.200 | 0.094 | 0.059 | 0.129 | 1.000 |
| scaffold | M | 0.200 | pool_lgbm | 0.065 | 0.112 | 0.072 | 0.156 | 1.000 |
| scaffold | M | 0.200 | pool_lgbm | 0.100 | 0.108 | 0.068 | 0.146 | 1.000 |
| scaffold | M | 0.200 | pool_lgbm | 0.200 | 0.115 | 0.081 | 0.148 | 1.000 |
| scaffold | M | 0.200 | rank_mean_gates | 0.065 | 0.082 | 0.046 | 0.123 | 1.000 |
| scaffold | M | 0.200 | rank_mean_gates | 0.100 | 0.087 | 0.054 | 0.121 | 1.000 |
| scaffold | M | 0.200 | rank_mean_gates | 0.200 | 0.125 | 0.095 | 0.156 | 1.000 |
| scaffold | M | 0.200 | l4_rank_mean | 0.065 | 0.087 | 0.056 | 0.123 | 1.000 |
| scaffold | M | 0.200 | l4_rank_mean | 0.100 | 0.081 | 0.049 | 0.109 | 1.000 |
| scaffold | M | 0.200 | l4_rank_mean | 0.200 | 0.081 | 0.054 | 0.105 | 1.000 |
| scaffold | M | 0.210 | gate_x_tdi | 0.065 | 0.081 | 0.050 | 0.120 | 1.000 |
| scaffold | M | 0.210 | gate_x_tdi | 0.100 | 0.077 | 0.047 | 0.111 | 1.000 |
| scaffold | M | 0.210 | gate_x_tdi | 0.200 | 0.102 | 0.075 | 0.129 | 1.000 |
| scaffold | M | 0.210 | gate_x_tdi_chemeleon | 0.065 | 0.075 | 0.033 | 0.121 | 1.000 |
| scaffold | M | 0.210 | gate_x_tdi_chemeleon | 0.100 | 0.077 | 0.039 | 0.115 | 1.000 |
| scaffold | M | 0.210 | gate_x_tdi_chemeleon | 0.200 | 0.100 | 0.066 | 0.133 | 1.000 |
| scaffold | M | 0.210 | legacy_sep15_refit | 0.065 | 0.036 | 0.002 | 0.069 | 0.977 |
| scaffold | M | 0.210 | legacy_sep15_refit | 0.100 | 0.032 | 0.001 | 0.061 | 0.977 |
| scaffold | M | 0.210 | legacy_sep15_refit | 0.200 | 0.035 | 0.009 | 0.060 | 0.998 |
| scaffold | M | 0.210 | pool_chemeleon | 0.065 | 0.056 | 0.014 | 0.100 | 0.995 |
| scaffold | M | 0.210 | pool_chemeleon | 0.100 | 0.055 | 0.016 | 0.095 | 0.999 |
| scaffold | M | 0.210 | pool_chemeleon | 0.200 | 0.095 | 0.059 | 0.130 | 1.000 |
| scaffold | M | 0.210 | pool_lgbm | 0.065 | 0.112 | 0.073 | 0.157 | 1.000 |
| scaffold | M | 0.210 | pool_lgbm | 0.100 | 0.108 | 0.067 | 0.147 | 1.000 |
| scaffold | M | 0.210 | pool_lgbm | 0.200 | 0.115 | 0.081 | 0.149 | 1.000 |
| scaffold | M | 0.210 | rank_mean_gates | 0.065 | 0.083 | 0.045 | 0.124 | 1.000 |
| scaffold | M | 0.210 | rank_mean_gates | 0.100 | 0.087 | 0.054 | 0.122 | 1.000 |
| scaffold | M | 0.210 | rank_mean_gates | 0.200 | 0.126 | 0.096 | 0.158 | 1.000 |
| scaffold | M | 0.210 | l4_rank_mean | 0.065 | 0.089 | 0.057 | 0.125 | 1.000 |
| scaffold | M | 0.210 | l4_rank_mean | 0.100 | 0.082 | 0.050 | 0.111 | 1.000 |
| scaffold | M | 0.210 | l4_rank_mean | 0.200 | 0.082 | 0.055 | 0.107 | 1.000 |
| scaffold | M | 0.214 | gate_x_tdi | 0.065 | 0.081 | 0.049 | 0.120 | 1.000 |
| scaffold | M | 0.214 | gate_x_tdi | 0.100 | 0.077 | 0.047 | 0.111 | 1.000 |
| scaffold | M | 0.214 | gate_x_tdi | 0.200 | 0.102 | 0.075 | 0.130 | 1.000 |
| scaffold | M | 0.214 | gate_x_tdi_chemeleon | 0.065 | 0.075 | 0.033 | 0.121 | 1.000 |
| scaffold | M | 0.214 | gate_x_tdi_chemeleon | 0.100 | 0.078 | 0.039 | 0.116 | 1.000 |
| scaffold | M | 0.214 | gate_x_tdi_chemeleon | 0.200 | 0.100 | 0.066 | 0.134 | 1.000 |
| scaffold | M | 0.214 | legacy_sep15_refit | 0.065 | 0.036 | 0.002 | 0.069 | 0.977 |
| scaffold | M | 0.214 | legacy_sep15_refit | 0.100 | 0.032 | 0.001 | 0.062 | 0.977 |
| scaffold | M | 0.214 | legacy_sep15_refit | 0.200 | 0.035 | 0.009 | 0.060 | 0.998 |
| scaffold | M | 0.214 | pool_chemeleon | 0.065 | 0.056 | 0.014 | 0.100 | 0.994 |
| scaffold | M | 0.214 | pool_chemeleon | 0.100 | 0.055 | 0.016 | 0.096 | 0.999 |
| scaffold | M | 0.214 | pool_chemeleon | 0.200 | 0.095 | 0.059 | 0.130 | 1.000 |
| scaffold | M | 0.214 | pool_lgbm | 0.065 | 0.113 | 0.073 | 0.157 | 1.000 |
| scaffold | M | 0.214 | pool_lgbm | 0.100 | 0.109 | 0.067 | 0.148 | 1.000 |
| scaffold | M | 0.214 | pool_lgbm | 0.200 | 0.116 | 0.081 | 0.150 | 1.000 |
| scaffold | M | 0.214 | rank_mean_gates | 0.065 | 0.083 | 0.045 | 0.124 | 1.000 |
| scaffold | M | 0.214 | rank_mean_gates | 0.100 | 0.087 | 0.054 | 0.122 | 1.000 |
| scaffold | M | 0.214 | rank_mean_gates | 0.200 | 0.127 | 0.096 | 0.158 | 1.000 |
| scaffold | M | 0.214 | l4_rank_mean | 0.065 | 0.089 | 0.057 | 0.125 | 1.000 |
| scaffold | M | 0.214 | l4_rank_mean | 0.100 | 0.082 | 0.050 | 0.111 | 1.000 |
| scaffold | M | 0.214 | l4_rank_mean | 0.200 | 0.083 | 0.055 | 0.107 | 1.000 |
| scaffold | M | 0.231 | gate_x_tdi | 0.065 | 0.081 | 0.048 | 0.121 | 1.000 |
| scaffold | M | 0.231 | gate_x_tdi | 0.100 | 0.078 | 0.046 | 0.112 | 1.000 |
| scaffold | M | 0.231 | gate_x_tdi | 0.200 | 0.103 | 0.075 | 0.131 | 1.000 |
| scaffold | M | 0.231 | gate_x_tdi_chemeleon | 0.065 | 0.075 | 0.031 | 0.122 | 0.999 |
| scaffold | M | 0.231 | gate_x_tdi_chemeleon | 0.100 | 0.078 | 0.038 | 0.117 | 1.000 |
| scaffold | M | 0.231 | gate_x_tdi_chemeleon | 0.200 | 0.102 | 0.066 | 0.136 | 1.000 |
| scaffold | M | 0.231 | legacy_sep15_refit | 0.065 | 0.037 | 0.002 | 0.072 | 0.977 |
| scaffold | M | 0.231 | legacy_sep15_refit | 0.100 | 0.033 | 0.001 | 0.064 | 0.978 |
| scaffold | M | 0.231 | legacy_sep15_refit | 0.200 | 0.037 | 0.010 | 0.063 | 0.998 |
| scaffold | M | 0.231 | pool_chemeleon | 0.065 | 0.056 | 0.012 | 0.102 | 0.993 |
| scaffold | M | 0.231 | pool_chemeleon | 0.100 | 0.056 | 0.015 | 0.097 | 0.999 |
| scaffold | M | 0.231 | pool_chemeleon | 0.200 | 0.097 | 0.059 | 0.133 | 1.000 |
| scaffold | M | 0.231 | pool_lgbm | 0.065 | 0.113 | 0.073 | 0.159 | 1.000 |
| scaffold | M | 0.231 | pool_lgbm | 0.100 | 0.109 | 0.067 | 0.149 | 1.000 |
| scaffold | M | 0.231 | pool_lgbm | 0.200 | 0.117 | 0.081 | 0.152 | 1.000 |
| scaffold | M | 0.231 | rank_mean_gates | 0.065 | 0.083 | 0.045 | 0.125 | 1.000 |
| scaffold | M | 0.231 | rank_mean_gates | 0.100 | 0.087 | 0.054 | 0.123 | 1.000 |
| scaffold | M | 0.231 | rank_mean_gates | 0.200 | 0.128 | 0.096 | 0.161 | 1.000 |
| scaffold | M | 0.231 | l4_rank_mean | 0.065 | 0.091 | 0.058 | 0.128 | 1.000 |
| scaffold | M | 0.231 | l4_rank_mean | 0.100 | 0.084 | 0.051 | 0.114 | 1.000 |
| scaffold | M | 0.231 | l4_rank_mean | 0.200 | 0.084 | 0.056 | 0.110 | 1.000 |
| scaffold | M | 0.250 | gate_x_tdi | 0.065 | 0.081 | 0.048 | 0.122 | 1.000 |
| scaffold | M | 0.250 | gate_x_tdi | 0.100 | 0.078 | 0.046 | 0.113 | 1.000 |
| scaffold | M | 0.250 | gate_x_tdi | 0.200 | 0.104 | 0.075 | 0.133 | 1.000 |
| scaffold | M | 0.250 | gate_x_tdi_chemeleon | 0.065 | 0.075 | 0.030 | 0.122 | 0.999 |
| scaffold | M | 0.250 | gate_x_tdi_chemeleon | 0.100 | 0.078 | 0.037 | 0.118 | 1.000 |
| scaffold | M | 0.250 | gate_x_tdi_chemeleon | 0.200 | 0.103 | 0.067 | 0.138 | 1.000 |
| scaffold | M | 0.250 | legacy_sep15_refit | 0.065 | 0.038 | 0.002 | 0.074 | 0.978 |
| scaffold | M | 0.250 | legacy_sep15_refit | 0.100 | 0.034 | 0.002 | 0.066 | 0.978 |
| scaffold | M | 0.250 | legacy_sep15_refit | 0.200 | 0.039 | 0.010 | 0.065 | 0.999 |
| scaffold | M | 0.250 | pool_chemeleon | 0.065 | 0.056 | 0.011 | 0.103 | 0.993 |
| scaffold | M | 0.250 | pool_chemeleon | 0.100 | 0.055 | 0.014 | 0.098 | 0.997 |
| scaffold | M | 0.250 | pool_chemeleon | 0.200 | 0.098 | 0.059 | 0.135 | 1.000 |
| scaffold | M | 0.250 | pool_lgbm | 0.065 | 0.114 | 0.072 | 0.161 | 1.000 |
| scaffold | M | 0.250 | pool_lgbm | 0.100 | 0.110 | 0.067 | 0.151 | 1.000 |
| scaffold | M | 0.250 | pool_lgbm | 0.200 | 0.118 | 0.081 | 0.154 | 1.000 |
| scaffold | M | 0.250 | rank_mean_gates | 0.065 | 0.083 | 0.044 | 0.125 | 1.000 |
| scaffold | M | 0.250 | rank_mean_gates | 0.100 | 0.088 | 0.053 | 0.124 | 1.000 |
| scaffold | M | 0.250 | rank_mean_gates | 0.200 | 0.130 | 0.097 | 0.164 | 1.000 |
| scaffold | M | 0.250 | l4_rank_mean | 0.065 | 0.093 | 0.059 | 0.131 | 1.000 |
| scaffold | M | 0.250 | l4_rank_mean | 0.100 | 0.086 | 0.052 | 0.116 | 1.000 |
| scaffold | M | 0.250 | l4_rank_mean | 0.200 | 0.086 | 0.057 | 0.113 | 1.000 |
| scaffold | M | 0.300 | gate_x_tdi | 0.065 | 0.080 | 0.044 | 0.122 | 1.000 |
| scaffold | M | 0.300 | gate_x_tdi | 0.100 | 0.077 | 0.042 | 0.113 | 1.000 |
| scaffold | M | 0.300 | gate_x_tdi | 0.200 | 0.105 | 0.075 | 0.136 | 1.000 |
| scaffold | M | 0.300 | gate_x_tdi_chemeleon | 0.065 | 0.073 | 0.026 | 0.121 | 0.998 |
| scaffold | M | 0.300 | gate_x_tdi_chemeleon | 0.100 | 0.077 | 0.034 | 0.120 | 0.999 |
| scaffold | M | 0.300 | gate_x_tdi_chemeleon | 0.200 | 0.104 | 0.065 | 0.142 | 1.000 |
| scaffold | M | 0.300 | legacy_sep15_refit | 0.065 | 0.042 | 0.002 | 0.080 | 0.980 |
| scaffold | M | 0.300 | legacy_sep15_refit | 0.100 | 0.038 | 0.002 | 0.072 | 0.981 |
| scaffold | M | 0.300 | legacy_sep15_refit | 0.200 | 0.043 | 0.012 | 0.072 | 0.999 |
| scaffold | M | 0.300 | pool_chemeleon | 0.065 | 0.054 | 0.006 | 0.105 | 0.987 |
| scaffold | M | 0.300 | pool_chemeleon | 0.100 | 0.054 | 0.008 | 0.100 | 0.995 |
| scaffold | M | 0.300 | pool_chemeleon | 0.200 | 0.099 | 0.057 | 0.140 | 1.000 |
| scaffold | M | 0.300 | pool_lgbm | 0.065 | 0.114 | 0.071 | 0.163 | 1.000 |
| scaffold | M | 0.300 | pool_lgbm | 0.100 | 0.110 | 0.063 | 0.153 | 1.000 |
| scaffold | M | 0.300 | pool_lgbm | 0.200 | 0.118 | 0.079 | 0.157 | 1.000 |
| scaffold | M | 0.300 | rank_mean_gates | 0.065 | 0.082 | 0.040 | 0.127 | 0.999 |
| scaffold | M | 0.300 | rank_mean_gates | 0.100 | 0.086 | 0.050 | 0.125 | 1.000 |
| scaffold | M | 0.300 | rank_mean_gates | 0.200 | 0.132 | 0.097 | 0.168 | 1.000 |
| scaffold | M | 0.300 | l4_rank_mean | 0.065 | 0.097 | 0.061 | 0.136 | 1.000 |
| scaffold | M | 0.300 | l4_rank_mean | 0.100 | 0.089 | 0.053 | 0.122 | 1.000 |
| scaffold | M | 0.300 | l4_rank_mean | 0.200 | 0.090 | 0.058 | 0.119 | 1.000 |
| scaffold | M | 0.308 | gate_x_tdi | 0.065 | 0.079 | 0.044 | 0.122 | 1.000 |
| scaffold | M | 0.308 | gate_x_tdi | 0.100 | 0.076 | 0.042 | 0.113 | 1.000 |
| scaffold | M | 0.308 | gate_x_tdi | 0.200 | 0.105 | 0.074 | 0.136 | 1.000 |
| scaffold | M | 0.308 | gate_x_tdi_chemeleon | 0.065 | 0.073 | 0.025 | 0.121 | 0.998 |
| scaffold | M | 0.308 | gate_x_tdi_chemeleon | 0.100 | 0.077 | 0.033 | 0.120 | 0.999 |
| scaffold | M | 0.308 | gate_x_tdi_chemeleon | 0.200 | 0.104 | 0.065 | 0.143 | 1.000 |
| scaffold | M | 0.308 | legacy_sep15_refit | 0.065 | 0.042 | 0.002 | 0.081 | 0.980 |
| scaffold | M | 0.308 | legacy_sep15_refit | 0.100 | 0.038 | 0.002 | 0.073 | 0.981 |
| scaffold | M | 0.308 | legacy_sep15_refit | 0.200 | 0.043 | 0.012 | 0.073 | 0.999 |
| scaffold | M | 0.308 | pool_chemeleon | 0.065 | 0.054 | 0.006 | 0.105 | 0.986 |
| scaffold | M | 0.308 | pool_chemeleon | 0.100 | 0.054 | 0.008 | 0.100 | 0.993 |
| scaffold | M | 0.308 | pool_chemeleon | 0.200 | 0.099 | 0.057 | 0.140 | 1.000 |
| scaffold | M | 0.308 | pool_lgbm | 0.065 | 0.114 | 0.071 | 0.163 | 1.000 |
| scaffold | M | 0.308 | pool_lgbm | 0.100 | 0.109 | 0.063 | 0.153 | 1.000 |
| scaffold | M | 0.308 | pool_lgbm | 0.200 | 0.118 | 0.078 | 0.157 | 1.000 |
| scaffold | M | 0.308 | rank_mean_gates | 0.065 | 0.081 | 0.040 | 0.127 | 0.999 |
| scaffold | M | 0.308 | rank_mean_gates | 0.100 | 0.086 | 0.049 | 0.125 | 1.000 |
| scaffold | M | 0.308 | rank_mean_gates | 0.200 | 0.132 | 0.096 | 0.168 | 1.000 |
| scaffold | M | 0.308 | l4_rank_mean | 0.065 | 0.097 | 0.061 | 0.137 | 1.000 |
| scaffold | M | 0.308 | l4_rank_mean | 0.100 | 0.090 | 0.053 | 0.123 | 1.000 |
| scaffold | M | 0.308 | l4_rank_mean | 0.200 | 0.091 | 0.059 | 0.120 | 1.000 |
| scaffold | M | 0.311 | gate_x_tdi | 0.065 | 0.079 | 0.044 | 0.122 | 1.000 |
| scaffold | M | 0.311 | gate_x_tdi | 0.100 | 0.076 | 0.041 | 0.113 | 1.000 |
| scaffold | M | 0.311 | gate_x_tdi | 0.200 | 0.105 | 0.074 | 0.136 | 1.000 |
| scaffold | M | 0.311 | gate_x_tdi_chemeleon | 0.065 | 0.073 | 0.024 | 0.121 | 0.998 |
| scaffold | M | 0.311 | gate_x_tdi_chemeleon | 0.100 | 0.077 | 0.033 | 0.120 | 0.999 |
| scaffold | M | 0.311 | gate_x_tdi_chemeleon | 0.200 | 0.104 | 0.064 | 0.143 | 1.000 |
| scaffold | M | 0.311 | legacy_sep15_refit | 0.065 | 0.042 | 0.002 | 0.081 | 0.980 |
| scaffold | M | 0.311 | legacy_sep15_refit | 0.100 | 0.038 | 0.002 | 0.073 | 0.981 |
| scaffold | M | 0.311 | legacy_sep15_refit | 0.200 | 0.043 | 0.012 | 0.073 | 0.999 |
| scaffold | M | 0.311 | pool_chemeleon | 0.065 | 0.054 | 0.005 | 0.105 | 0.986 |
| scaffold | M | 0.311 | pool_chemeleon | 0.100 | 0.054 | 0.007 | 0.100 | 0.993 |
| scaffold | M | 0.311 | pool_chemeleon | 0.200 | 0.099 | 0.057 | 0.140 | 1.000 |
| scaffold | M | 0.311 | pool_lgbm | 0.065 | 0.114 | 0.071 | 0.163 | 1.000 |
| scaffold | M | 0.311 | pool_lgbm | 0.100 | 0.109 | 0.062 | 0.153 | 1.000 |
| scaffold | M | 0.311 | pool_lgbm | 0.200 | 0.118 | 0.078 | 0.157 | 1.000 |
| scaffold | M | 0.311 | rank_mean_gates | 0.065 | 0.081 | 0.039 | 0.127 | 0.999 |
| scaffold | M | 0.311 | rank_mean_gates | 0.100 | 0.086 | 0.049 | 0.125 | 1.000 |
| scaffold | M | 0.311 | rank_mean_gates | 0.200 | 0.132 | 0.096 | 0.168 | 1.000 |
| scaffold | M | 0.311 | l4_rank_mean | 0.065 | 0.097 | 0.061 | 0.137 | 1.000 |
| scaffold | M | 0.311 | l4_rank_mean | 0.100 | 0.090 | 0.053 | 0.124 | 1.000 |
| scaffold | M | 0.311 | l4_rank_mean | 0.200 | 0.091 | 0.059 | 0.120 | 1.000 |
| scaffold | M | 0.326 | gate_x_tdi | 0.065 | 0.078 | 0.042 | 0.122 | 1.000 |
| scaffold | M | 0.326 | gate_x_tdi | 0.100 | 0.076 | 0.040 | 0.113 | 1.000 |
| scaffold | M | 0.326 | gate_x_tdi | 0.200 | 0.105 | 0.073 | 0.137 | 1.000 |
| scaffold | M | 0.326 | gate_x_tdi_chemeleon | 0.065 | 0.072 | 0.023 | 0.120 | 0.997 |
| scaffold | M | 0.326 | gate_x_tdi_chemeleon | 0.100 | 0.076 | 0.031 | 0.120 | 0.999 |
| scaffold | M | 0.326 | gate_x_tdi_chemeleon | 0.200 | 0.104 | 0.064 | 0.144 | 1.000 |
| scaffold | M | 0.326 | legacy_sep15_refit | 0.065 | 0.043 | 0.002 | 0.083 | 0.980 |
| scaffold | M | 0.326 | legacy_sep15_refit | 0.100 | 0.039 | 0.002 | 0.075 | 0.983 |
| scaffold | M | 0.326 | legacy_sep15_refit | 0.200 | 0.045 | 0.013 | 0.075 | 0.999 |
| scaffold | M | 0.326 | pool_chemeleon | 0.065 | 0.053 | 0.004 | 0.105 | 0.983 |
| scaffold | M | 0.326 | pool_chemeleon | 0.100 | 0.053 | 0.006 | 0.100 | 0.991 |
| scaffold | M | 0.326 | pool_chemeleon | 0.200 | 0.099 | 0.056 | 0.141 | 1.000 |
| scaffold | M | 0.326 | pool_lgbm | 0.065 | 0.114 | 0.070 | 0.163 | 1.000 |
| scaffold | M | 0.326 | pool_lgbm | 0.100 | 0.109 | 0.062 | 0.153 | 1.000 |
| scaffold | M | 0.326 | pool_lgbm | 0.200 | 0.118 | 0.077 | 0.157 | 1.000 |
| scaffold | M | 0.326 | rank_mean_gates | 0.065 | 0.080 | 0.038 | 0.127 | 0.999 |
| scaffold | M | 0.326 | rank_mean_gates | 0.100 | 0.085 | 0.047 | 0.125 | 1.000 |
| scaffold | M | 0.326 | rank_mean_gates | 0.200 | 0.132 | 0.096 | 0.169 | 1.000 |
| scaffold | M | 0.326 | l4_rank_mean | 0.065 | 0.098 | 0.062 | 0.139 | 1.000 |
| scaffold | M | 0.326 | l4_rank_mean | 0.100 | 0.091 | 0.054 | 0.125 | 1.000 |
| scaffold | M | 0.326 | l4_rank_mean | 0.200 | 0.092 | 0.059 | 0.121 | 1.000 |
| scaffold | M | 0.400 | gate_x_tdi | 0.065 | 0.073 | 0.034 | 0.118 | 1.000 |
| scaffold | M | 0.400 | gate_x_tdi | 0.100 | 0.071 | 0.034 | 0.111 | 0.998 |
| scaffold | M | 0.400 | gate_x_tdi | 0.200 | 0.102 | 0.069 | 0.136 | 1.000 |
| scaffold | M | 0.400 | gate_x_tdi_chemeleon | 0.065 | 0.067 | 0.015 | 0.118 | 0.993 |
| scaffold | M | 0.400 | gate_x_tdi_chemeleon | 0.100 | 0.071 | 0.023 | 0.118 | 0.998 |
| scaffold | M | 0.400 | gate_x_tdi_chemeleon | 0.200 | 0.101 | 0.058 | 0.145 | 1.000 |
| scaffold | M | 0.400 | legacy_sep15_refit | 0.065 | 0.047 | 0.002 | 0.089 | 0.979 |
| scaffold | M | 0.400 | legacy_sep15_refit | 0.100 | 0.043 | 0.003 | 0.082 | 0.984 |
| scaffold | M | 0.400 | legacy_sep15_refit | 0.200 | 0.050 | 0.016 | 0.084 | 0.999 |
| scaffold | M | 0.400 | pool_chemeleon | 0.065 | 0.048 | -0.004 | 0.104 | 0.960 |
| scaffold | M | 0.400 | pool_chemeleon | 0.100 | 0.048 | -0.003 | 0.098 | 0.973 |
| scaffold | M | 0.400 | pool_chemeleon | 0.200 | 0.096 | 0.050 | 0.142 | 1.000 |
| scaffold | M | 0.400 | pool_lgbm | 0.065 | 0.110 | 0.065 | 0.161 | 1.000 |
| scaffold | M | 0.400 | pool_lgbm | 0.100 | 0.105 | 0.056 | 0.151 | 1.000 |
| scaffold | M | 0.400 | pool_lgbm | 0.200 | 0.113 | 0.069 | 0.156 | 1.000 |
| scaffold | M | 0.400 | rank_mean_gates | 0.065 | 0.075 | 0.030 | 0.124 | 0.999 |
| scaffold | M | 0.400 | rank_mean_gates | 0.100 | 0.080 | 0.040 | 0.122 | 1.000 |
| scaffold | M | 0.400 | rank_mean_gates | 0.200 | 0.130 | 0.092 | 0.170 | 1.000 |
| scaffold | M | 0.400 | l4_rank_mean | 0.065 | 0.101 | 0.062 | 0.145 | 1.000 |
| scaffold | M | 0.400 | l4_rank_mean | 0.100 | 0.093 | 0.053 | 0.129 | 1.000 |
| scaffold | M | 0.400 | l4_rank_mean | 0.200 | 0.094 | 0.059 | 0.126 | 1.000 |
| scaffold | M | 0.400 | gate_x_tdi | 0.065 | 0.073 | 0.034 | 0.118 | 1.000 |
| scaffold | M | 0.400 | gate_x_tdi | 0.100 | 0.071 | 0.034 | 0.111 | 0.998 |
| scaffold | M | 0.400 | gate_x_tdi | 0.200 | 0.102 | 0.069 | 0.136 | 1.000 |
| scaffold | M | 0.400 | gate_x_tdi_chemeleon | 0.065 | 0.067 | 0.015 | 0.118 | 0.993 |
| scaffold | M | 0.400 | gate_x_tdi_chemeleon | 0.100 | 0.071 | 0.023 | 0.118 | 0.998 |
| scaffold | M | 0.400 | gate_x_tdi_chemeleon | 0.200 | 0.101 | 0.058 | 0.145 | 1.000 |
| scaffold | M | 0.400 | legacy_sep15_refit | 0.065 | 0.047 | 0.002 | 0.089 | 0.979 |
| scaffold | M | 0.400 | legacy_sep15_refit | 0.100 | 0.043 | 0.003 | 0.082 | 0.984 |
| scaffold | M | 0.400 | legacy_sep15_refit | 0.200 | 0.050 | 0.016 | 0.084 | 0.999 |
| scaffold | M | 0.400 | pool_chemeleon | 0.065 | 0.048 | -0.004 | 0.104 | 0.960 |
| scaffold | M | 0.400 | pool_chemeleon | 0.100 | 0.048 | -0.003 | 0.098 | 0.973 |
| scaffold | M | 0.400 | pool_chemeleon | 0.200 | 0.096 | 0.050 | 0.142 | 1.000 |
| scaffold | M | 0.400 | pool_lgbm | 0.065 | 0.110 | 0.065 | 0.161 | 1.000 |
| scaffold | M | 0.400 | pool_lgbm | 0.100 | 0.105 | 0.056 | 0.151 | 1.000 |
| scaffold | M | 0.400 | pool_lgbm | 0.200 | 0.113 | 0.069 | 0.156 | 1.000 |
| scaffold | M | 0.400 | rank_mean_gates | 0.065 | 0.075 | 0.030 | 0.124 | 0.999 |
| scaffold | M | 0.400 | rank_mean_gates | 0.100 | 0.080 | 0.040 | 0.122 | 1.000 |
| scaffold | M | 0.400 | rank_mean_gates | 0.200 | 0.130 | 0.092 | 0.170 | 1.000 |
| scaffold | M | 0.400 | l4_rank_mean | 0.065 | 0.101 | 0.062 | 0.145 | 1.000 |
| scaffold | M | 0.400 | l4_rank_mean | 0.100 | 0.093 | 0.053 | 0.129 | 1.000 |
| scaffold | M | 0.400 | l4_rank_mean | 0.200 | 0.094 | 0.059 | 0.126 | 1.000 |
| scaffold | M | 0.403 | gate_x_tdi | 0.065 | 0.073 | 0.034 | 0.118 | 1.000 |
| scaffold | M | 0.403 | gate_x_tdi | 0.100 | 0.071 | 0.034 | 0.111 | 0.998 |
| scaffold | M | 0.403 | gate_x_tdi | 0.200 | 0.102 | 0.068 | 0.136 | 1.000 |
| scaffold | M | 0.403 | gate_x_tdi_chemeleon | 0.065 | 0.066 | 0.015 | 0.118 | 0.993 |
| scaffold | M | 0.403 | gate_x_tdi_chemeleon | 0.100 | 0.071 | 0.023 | 0.118 | 0.998 |
| scaffold | M | 0.403 | gate_x_tdi_chemeleon | 0.200 | 0.101 | 0.058 | 0.145 | 1.000 |
| scaffold | M | 0.403 | legacy_sep15_refit | 0.065 | 0.047 | 0.002 | 0.089 | 0.979 |
| scaffold | M | 0.403 | legacy_sep15_refit | 0.100 | 0.044 | 0.003 | 0.082 | 0.984 |
| scaffold | M | 0.403 | legacy_sep15_refit | 0.200 | 0.050 | 0.016 | 0.084 | 0.999 |
| scaffold | M | 0.403 | pool_chemeleon | 0.065 | 0.048 | -0.004 | 0.103 | 0.958 |
| scaffold | M | 0.403 | pool_chemeleon | 0.100 | 0.048 | -0.003 | 0.098 | 0.973 |
| scaffold | M | 0.403 | pool_chemeleon | 0.200 | 0.096 | 0.050 | 0.142 | 1.000 |
| scaffold | M | 0.403 | pool_lgbm | 0.065 | 0.110 | 0.064 | 0.161 | 1.000 |
| scaffold | M | 0.403 | pool_lgbm | 0.100 | 0.104 | 0.056 | 0.151 | 1.000 |
| scaffold | M | 0.403 | pool_lgbm | 0.200 | 0.113 | 0.069 | 0.156 | 1.000 |
| scaffold | M | 0.403 | rank_mean_gates | 0.065 | 0.075 | 0.029 | 0.124 | 0.999 |
| scaffold | M | 0.403 | rank_mean_gates | 0.100 | 0.080 | 0.039 | 0.122 | 1.000 |
| scaffold | M | 0.403 | rank_mean_gates | 0.200 | 0.130 | 0.092 | 0.170 | 1.000 |
| scaffold | M | 0.403 | l4_rank_mean | 0.065 | 0.101 | 0.062 | 0.145 | 1.000 |
| scaffold | M | 0.403 | l4_rank_mean | 0.100 | 0.093 | 0.053 | 0.129 | 1.000 |
| scaffold | M | 0.403 | l4_rank_mean | 0.200 | 0.094 | 0.059 | 0.126 | 1.000 |
| scaffold | M | 0.416 | gate_x_tdi | 0.065 | 0.072 | 0.033 | 0.118 | 1.000 |
| scaffold | M | 0.416 | gate_x_tdi | 0.100 | 0.070 | 0.033 | 0.110 | 0.998 |
| scaffold | M | 0.416 | gate_x_tdi | 0.200 | 0.101 | 0.067 | 0.135 | 1.000 |
| scaffold | M | 0.416 | gate_x_tdi_chemeleon | 0.065 | 0.065 | 0.014 | 0.117 | 0.993 |
| scaffold | M | 0.416 | gate_x_tdi_chemeleon | 0.100 | 0.070 | 0.021 | 0.118 | 0.997 |
| scaffold | M | 0.416 | gate_x_tdi_chemeleon | 0.200 | 0.101 | 0.057 | 0.144 | 1.000 |
| scaffold | M | 0.416 | legacy_sep15_refit | 0.065 | 0.048 | 0.002 | 0.091 | 0.979 |
| scaffold | M | 0.416 | legacy_sep15_refit | 0.100 | 0.044 | 0.003 | 0.083 | 0.984 |
| scaffold | M | 0.416 | legacy_sep15_refit | 0.200 | 0.051 | 0.016 | 0.086 | 0.999 |
| scaffold | M | 0.416 | pool_chemeleon | 0.065 | 0.047 | -0.006 | 0.103 | 0.951 |
| scaffold | M | 0.416 | pool_chemeleon | 0.100 | 0.047 | -0.005 | 0.097 | 0.971 |
| scaffold | M | 0.416 | pool_chemeleon | 0.200 | 0.095 | 0.048 | 0.142 | 1.000 |
| scaffold | M | 0.416 | pool_lgbm | 0.065 | 0.109 | 0.063 | 0.161 | 1.000 |
| scaffold | M | 0.416 | pool_lgbm | 0.100 | 0.103 | 0.054 | 0.150 | 1.000 |
| scaffold | M | 0.416 | pool_lgbm | 0.200 | 0.112 | 0.067 | 0.155 | 1.000 |
| scaffold | M | 0.416 | rank_mean_gates | 0.065 | 0.074 | 0.028 | 0.124 | 0.999 |
| scaffold | M | 0.416 | rank_mean_gates | 0.100 | 0.079 | 0.038 | 0.121 | 1.000 |
| scaffold | M | 0.416 | rank_mean_gates | 0.200 | 0.130 | 0.091 | 0.170 | 1.000 |
| scaffold | M | 0.416 | l4_rank_mean | 0.065 | 0.101 | 0.062 | 0.145 | 1.000 |
| scaffold | M | 0.416 | l4_rank_mean | 0.100 | 0.093 | 0.053 | 0.130 | 1.000 |
| scaffold | M | 0.416 | l4_rank_mean | 0.200 | 0.094 | 0.059 | 0.127 | 1.000 |
| scaffold | M | 0.500 | gate_x_tdi | 0.065 | 0.064 | 0.023 | 0.112 | 0.999 |
| scaffold | M | 0.500 | gate_x_tdi | 0.100 | 0.061 | 0.022 | 0.105 | 0.997 |
| scaffold | M | 0.500 | gate_x_tdi | 0.200 | 0.094 | 0.056 | 0.129 | 1.000 |
| scaffold | M | 0.500 | gate_x_tdi_chemeleon | 0.065 | 0.057 | 0.002 | 0.113 | 0.979 |
| scaffold | M | 0.500 | gate_x_tdi_chemeleon | 0.100 | 0.061 | 0.009 | 0.113 | 0.991 |
| scaffold | M | 0.500 | gate_x_tdi_chemeleon | 0.200 | 0.094 | 0.046 | 0.141 | 1.000 |
| scaffold | M | 0.500 | legacy_sep15_refit | 0.065 | 0.051 | 0.002 | 0.096 | 0.979 |
| scaffold | M | 0.500 | legacy_sep15_refit | 0.100 | 0.048 | 0.003 | 0.090 | 0.982 |
| scaffold | M | 0.500 | legacy_sep15_refit | 0.200 | 0.057 | 0.019 | 0.095 | 0.999 |
| scaffold | M | 0.500 | pool_chemeleon | 0.065 | 0.039 | -0.017 | 0.099 | 0.905 |
| scaffold | M | 0.500 | pool_chemeleon | 0.100 | 0.039 | -0.016 | 0.093 | 0.928 |
| scaffold | M | 0.500 | pool_chemeleon | 0.200 | 0.088 | 0.039 | 0.140 | 1.000 |
| scaffold | M | 0.500 | pool_lgbm | 0.065 | 0.103 | 0.054 | 0.157 | 1.000 |
| scaffold | M | 0.500 | pool_lgbm | 0.100 | 0.095 | 0.043 | 0.146 | 1.000 |
| scaffold | M | 0.500 | pool_lgbm | 0.200 | 0.103 | 0.057 | 0.149 | 1.000 |
| scaffold | M | 0.500 | rank_mean_gates | 0.065 | 0.066 | 0.016 | 0.118 | 0.995 |
| scaffold | M | 0.500 | rank_mean_gates | 0.100 | 0.070 | 0.026 | 0.115 | 0.997 |
| scaffold | M | 0.500 | rank_mean_gates | 0.200 | 0.123 | 0.082 | 0.165 | 1.000 |
| scaffold | M | 0.500 | l4_rank_mean | 0.065 | 0.101 | 0.060 | 0.148 | 1.000 |
| scaffold | M | 0.500 | l4_rank_mean | 0.100 | 0.093 | 0.050 | 0.132 | 1.000 |
| scaffold | M | 0.500 | l4_rank_mean | 0.200 | 0.094 | 0.056 | 0.130 | 1.000 |
| scaffold | M | 0.600 | gate_x_tdi | 0.065 | 0.052 | 0.008 | 0.105 | 0.991 |
| scaffold | M | 0.600 | gate_x_tdi | 0.100 | 0.049 | 0.006 | 0.095 | 0.989 |
| scaffold | M | 0.600 | gate_x_tdi | 0.200 | 0.081 | 0.041 | 0.119 | 1.000 |
| scaffold | M | 0.600 | gate_x_tdi_chemeleon | 0.065 | 0.044 | -0.014 | 0.106 | 0.924 |
| scaffold | M | 0.600 | gate_x_tdi_chemeleon | 0.100 | 0.049 | -0.008 | 0.104 | 0.955 |
| scaffold | M | 0.600 | gate_x_tdi_chemeleon | 0.200 | 0.082 | 0.030 | 0.133 | 1.000 |
| scaffold | M | 0.600 | legacy_sep15_refit | 0.065 | 0.055 | 0.001 | 0.103 | 0.978 |
| scaffold | M | 0.600 | legacy_sep15_refit | 0.100 | 0.053 | 0.003 | 0.098 | 0.982 |
| scaffold | M | 0.600 | legacy_sep15_refit | 0.200 | 0.064 | 0.022 | 0.104 | 1.000 |
| scaffold | M | 0.600 | pool_chemeleon | 0.065 | 0.028 | -0.031 | 0.092 | 0.819 |
| scaffold | M | 0.600 | pool_chemeleon | 0.100 | 0.027 | -0.031 | 0.084 | 0.832 |
| scaffold | M | 0.600 | pool_chemeleon | 0.200 | 0.077 | 0.023 | 0.133 | 0.997 |
| scaffold | M | 0.600 | pool_lgbm | 0.065 | 0.092 | 0.042 | 0.149 | 0.999 |
| scaffold | M | 0.600 | pool_lgbm | 0.100 | 0.083 | 0.027 | 0.137 | 1.000 |
| scaffold | M | 0.600 | pool_lgbm | 0.200 | 0.088 | 0.038 | 0.137 | 1.000 |
| scaffold | M | 0.600 | rank_mean_gates | 0.065 | 0.054 | 0.002 | 0.109 | 0.981 |
| scaffold | M | 0.600 | rank_mean_gates | 0.100 | 0.057 | 0.010 | 0.105 | 0.993 |
| scaffold | M | 0.600 | rank_mean_gates | 0.200 | 0.111 | 0.067 | 0.156 | 1.000 |
| scaffold | M | 0.600 | l4_rank_mean | 0.065 | 0.099 | 0.055 | 0.149 | 1.000 |
| scaffold | M | 0.600 | l4_rank_mean | 0.100 | 0.091 | 0.045 | 0.132 | 1.000 |
| scaffold | M | 0.600 | l4_rank_mean | 0.200 | 0.091 | 0.050 | 0.130 | 1.000 |
| scaffold | M | 0.700 | gate_x_tdi | 0.065 | 0.038 | -0.008 | 0.094 | 0.945 |
| scaffold | M | 0.700 | gate_x_tdi | 0.100 | 0.035 | -0.009 | 0.084 | 0.943 |
| scaffold | M | 0.700 | gate_x_tdi | 0.200 | 0.066 | 0.023 | 0.105 | 1.000 |
| scaffold | M | 0.700 | gate_x_tdi_chemeleon | 0.065 | 0.030 | -0.031 | 0.094 | 0.830 |
| scaffold | M | 0.700 | gate_x_tdi_chemeleon | 0.100 | 0.034 | -0.027 | 0.094 | 0.880 |
| scaffold | M | 0.700 | gate_x_tdi_chemeleon | 0.200 | 0.066 | 0.011 | 0.121 | 0.986 |
| scaffold | M | 0.700 | legacy_sep15_refit | 0.065 | 0.058 | 0.000 | 0.108 | 0.975 |
| scaffold | M | 0.700 | legacy_sep15_refit | 0.100 | 0.057 | 0.003 | 0.104 | 0.983 |
| scaffold | M | 0.700 | legacy_sep15_refit | 0.200 | 0.070 | 0.025 | 0.114 | 1.000 |
| scaffold | M | 0.700 | pool_chemeleon | 0.065 | 0.015 | -0.047 | 0.083 | 0.681 |
| scaffold | M | 0.700 | pool_chemeleon | 0.100 | 0.014 | -0.047 | 0.073 | 0.678 |
| scaffold | M | 0.700 | pool_chemeleon | 0.200 | 0.062 | 0.005 | 0.121 | 0.981 |
| scaffold | M | 0.700 | pool_lgbm | 0.065 | 0.080 | 0.027 | 0.140 | 0.998 |
| scaffold | M | 0.700 | pool_lgbm | 0.100 | 0.068 | 0.009 | 0.126 | 0.988 |
| scaffold | M | 0.700 | pool_lgbm | 0.200 | 0.069 | 0.017 | 0.121 | 0.996 |
| scaffold | M | 0.700 | rank_mean_gates | 0.065 | 0.041 | -0.015 | 0.098 | 0.934 |
| scaffold | M | 0.700 | rank_mean_gates | 0.100 | 0.042 | -0.006 | 0.091 | 0.946 |
| scaffold | M | 0.700 | rank_mean_gates | 0.200 | 0.095 | 0.049 | 0.144 | 1.000 |
| scaffold | M | 0.700 | l4_rank_mean | 0.065 | 0.095 | 0.048 | 0.146 | 0.999 |
| scaffold | M | 0.700 | l4_rank_mean | 0.100 | 0.086 | 0.038 | 0.130 | 1.000 |
| scaffold | M | 0.700 | l4_rank_mean | 0.200 | 0.086 | 0.042 | 0.129 | 1.000 |
| scaffold | M | 0.800 | gate_x_tdi | 0.065 | 0.024 | -0.025 | 0.081 | 0.832 |
| scaffold | M | 0.800 | gate_x_tdi | 0.100 | 0.019 | -0.027 | 0.069 | 0.799 |
| scaffold | M | 0.800 | gate_x_tdi | 0.200 | 0.047 | 0.002 | 0.088 | 0.982 |
| scaffold | M | 0.800 | gate_x_tdi_chemeleon | 0.065 | 0.015 | -0.051 | 0.080 | 0.692 |
| scaffold | M | 0.800 | gate_x_tdi_chemeleon | 0.100 | 0.018 | -0.046 | 0.080 | 0.743 |
| scaffold | M | 0.800 | gate_x_tdi_chemeleon | 0.200 | 0.048 | -0.010 | 0.108 | 0.936 |
| scaffold | M | 0.800 | legacy_sep15_refit | 0.065 | 0.061 | -0.001 | 0.114 | 0.971 |
| scaffold | M | 0.800 | legacy_sep15_refit | 0.100 | 0.061 | 0.004 | 0.111 | 0.983 |
| scaffold | M | 0.800 | legacy_sep15_refit | 0.200 | 0.077 | 0.029 | 0.122 | 1.000 |
| scaffold | M | 0.800 | pool_chemeleon | 0.065 | 0.001 | -0.065 | 0.072 | 0.526 |
| scaffold | M | 0.800 | pool_chemeleon | 0.100 | -0.001 | -0.065 | 0.061 | 0.478 |
| scaffold | M | 0.800 | pool_chemeleon | 0.200 | 0.044 | -0.017 | 0.106 | 0.927 |
| scaffold | M | 0.800 | pool_lgbm | 0.065 | 0.066 | 0.011 | 0.131 | 0.989 |
| scaffold | M | 0.800 | pool_lgbm | 0.100 | 0.052 | -0.010 | 0.113 | 0.957 |
| scaffold | M | 0.800 | pool_lgbm | 0.200 | 0.047 | -0.007 | 0.100 | 0.956 |
| scaffold | M | 0.800 | rank_mean_gates | 0.065 | 0.026 | -0.033 | 0.085 | 0.817 |
| scaffold | M | 0.800 | rank_mean_gates | 0.100 | 0.025 | -0.026 | 0.075 | 0.835 |
| scaffold | M | 0.800 | rank_mean_gates | 0.200 | 0.076 | 0.026 | 0.127 | 0.996 |
| scaffold | M | 0.800 | l4_rank_mean | 0.065 | 0.090 | 0.041 | 0.143 | 0.999 |
| scaffold | M | 0.800 | l4_rank_mean | 0.100 | 0.080 | 0.029 | 0.126 | 1.000 |
| scaffold | M | 0.800 | l4_rank_mean | 0.200 | 0.078 | 0.031 | 0.124 | 0.999 |
| scaffold | M | 0.900 | gate_x_tdi | 0.065 | 0.008 | -0.043 | 0.065 | 0.629 |
| scaffold | M | 0.900 | gate_x_tdi | 0.100 | 0.002 | -0.047 | 0.055 | 0.544 |
| scaffold | M | 0.900 | gate_x_tdi | 0.200 | 0.026 | -0.021 | 0.071 | 0.859 |
| scaffold | M | 0.900 | gate_x_tdi_chemeleon | 0.065 | -0.001 | -0.072 | 0.066 | 0.491 |
| scaffold | M | 0.900 | gate_x_tdi_chemeleon | 0.100 | -0.000 | -0.068 | 0.064 | 0.500 |
| scaffold | M | 0.900 | gate_x_tdi_chemeleon | 0.200 | 0.027 | -0.036 | 0.089 | 0.815 |
| scaffold | M | 0.900 | legacy_sep15_refit | 0.065 | 0.063 | -0.003 | 0.119 | 0.968 |
| scaffold | M | 0.900 | legacy_sep15_refit | 0.100 | 0.065 | 0.003 | 0.118 | 0.981 |
| scaffold | M | 0.900 | legacy_sep15_refit | 0.200 | 0.083 | 0.031 | 0.131 | 1.000 |
| scaffold | M | 0.900 | pool_chemeleon | 0.065 | -0.014 | -0.084 | 0.060 | 0.370 |
| scaffold | M | 0.900 | pool_chemeleon | 0.100 | -0.017 | -0.085 | 0.049 | 0.314 |
| scaffold | M | 0.900 | pool_chemeleon | 0.200 | 0.024 | -0.042 | 0.090 | 0.782 |
| scaffold | M | 0.900 | pool_lgbm | 0.065 | 0.051 | -0.007 | 0.118 | 0.957 |
| scaffold | M | 0.900 | pool_lgbm | 0.100 | 0.034 | -0.030 | 0.097 | 0.870 |
| scaffold | M | 0.900 | pool_lgbm | 0.200 | 0.023 | -0.036 | 0.076 | 0.784 |
| scaffold | M | 0.900 | rank_mean_gates | 0.065 | 0.010 | -0.053 | 0.071 | 0.642 |
| scaffold | M | 0.900 | rank_mean_gates | 0.100 | 0.007 | -0.047 | 0.060 | 0.613 |
| scaffold | M | 0.900 | rank_mean_gates | 0.200 | 0.054 | -0.000 | 0.108 | 0.973 |
| scaffold | M | 0.900 | l4_rank_mean | 0.065 | 0.084 | 0.032 | 0.139 | 0.998 |
| scaffold | M | 0.900 | l4_rank_mean | 0.100 | 0.073 | 0.019 | 0.122 | 0.998 |
| scaffold | M | 0.900 | l4_rank_mean | 0.200 | 0.068 | 0.018 | 0.117 | 0.997 |
| scaffold | M | 1.000 | gate_x_tdi | 0.065 | -0.009 | -0.063 | 0.040 | 0.390 |
| scaffold | M | 1.000 | gate_x_tdi | 0.100 | -0.016 | -0.067 | 0.038 | 0.283 |
| scaffold | M | 1.000 | gate_x_tdi | 0.200 | 0.002 | -0.045 | 0.051 | 0.564 |
| scaffold | M | 1.000 | gate_x_tdi_chemeleon | 0.065 | -0.019 | -0.089 | 0.053 | 0.290 |
| scaffold | M | 1.000 | gate_x_tdi_chemeleon | 0.100 | -0.020 | -0.084 | 0.047 | 0.287 |
| scaffold | M | 1.000 | gate_x_tdi_chemeleon | 0.200 | 0.003 | -0.061 | 0.073 | 0.553 |
| scaffold | M | 1.000 | legacy_sep15_refit | 0.065 | 0.064 | 0.000 | 0.130 | 0.975 |
| scaffold | M | 1.000 | legacy_sep15_refit | 0.100 | 0.068 | 0.012 | 0.127 | 0.991 |
| scaffold | M | 1.000 | legacy_sep15_refit | 0.200 | 0.089 | 0.034 | 0.145 | 0.998 |
| scaffold | M | 1.000 | pool_chemeleon | 0.065 | -0.030 | -0.101 | 0.042 | 0.231 |
| scaffold | M | 1.000 | pool_chemeleon | 0.100 | -0.035 | -0.104 | 0.033 | 0.161 |
| scaffold | M | 1.000 | pool_chemeleon | 0.200 | 0.001 | -0.066 | 0.070 | 0.500 |
| scaffold | M | 1.000 | pool_lgbm | 0.065 | 0.035 | -0.025 | 0.097 | 0.854 |
| scaffold | M | 1.000 | pool_lgbm | 0.100 | 0.014 | -0.044 | 0.076 | 0.690 |
| scaffold | M | 1.000 | pool_lgbm | 0.200 | -0.005 | -0.060 | 0.055 | 0.466 |
| scaffold | M | 1.000 | rank_mean_gates | 0.065 | -0.007 | -0.067 | 0.056 | 0.432 |
| scaffold | M | 1.000 | rank_mean_gates | 0.100 | -0.013 | -0.068 | 0.046 | 0.364 |
| scaffold | M | 1.000 | rank_mean_gates | 0.200 | 0.029 | -0.023 | 0.088 | 0.854 |
| scaffold | M | 1.000 | l4_rank_mean | 0.065 | 0.076 | 0.022 | 0.132 | 0.995 |
| scaffold | M | 1.000 | l4_rank_mean | 0.100 | 0.064 | 0.011 | 0.116 | 0.993 |
| scaffold | M | 1.000 | l4_rank_mean | 0.200 | 0.056 | 0.004 | 0.108 | 0.981 |

## Crossover mixture fraction

| split_id | model | rate | delta_at_identified_phi | delta_at_phi_1 | crossings |
|---|---|---|---|---|---|
| umap | gate_x_tdi | 0.065 | 0.049 | -0.036 | 0.750 |
| umap | gate_x_tdi | 0.100 | 0.050 | -0.068 | 0.650 |
| umap | gate_x_tdi | 0.200 | 0.054 | -0.064 | 0.750 |
| umap | gate_x_tdi_chemeleon | 0.065 | 0.044 | -0.031 | 0.750 |
| umap | gate_x_tdi_chemeleon | 0.100 | 0.051 | -0.050 | 0.750 |
| umap | gate_x_tdi_chemeleon | 0.200 | 0.058 | -0.047 | 0.750 |
| umap | legacy_sep15_refit | 0.065 | 0.037 | 0.085 | none |
| umap | legacy_sep15_refit | 0.100 | 0.028 | 0.041 | none |
| umap | legacy_sep15_refit | 0.200 | 0.026 | 0.039 | none |
| umap | pool_chemeleon | 0.065 | 0.018 | -0.066 | 0.550 |
| umap | pool_chemeleon | 0.100 | 0.032 | -0.083 | 0.550 |
| umap | pool_chemeleon | 0.200 | 0.054 | -0.047 | 0.750 |
| umap | pool_lgbm | 0.065 | 0.064 | -0.031 | 0.850 |
| umap | pool_lgbm | 0.100 | 0.064 | -0.065 | 0.750 |
| umap | pool_lgbm | 0.200 | 0.091 | -0.032 | 0.850 |
| umap | rank_mean_gates | 0.065 | 0.056 | -0.029 | 0.850 |
| umap | rank_mean_gates | 0.100 | 0.058 | -0.069 | 0.650 |
| umap | rank_mean_gates | 0.200 | 0.067 | -0.055 | 0.750 |
| umap | l4_rank_mean | 0.065 | 0.068 | 0.044 | none |
| umap | l4_rank_mean | 0.100 | 0.075 | 0.024 | none |
| umap | l4_rank_mean | 0.200 | 0.071 | 0.018 | none |
| scaffold | gate_x_tdi | 0.065 | 0.078 | -0.009 | 0.950 |
| scaffold | gate_x_tdi | 0.100 | 0.076 | -0.016 | 0.950 |
| scaffold | gate_x_tdi | 0.200 | 0.105 | 0.002 | none |
| scaffold | gate_x_tdi_chemeleon | 0.065 | 0.072 | -0.019 | 0.850 |
| scaffold | gate_x_tdi_chemeleon | 0.100 | 0.076 | -0.020 | 0.850 |
| scaffold | gate_x_tdi_chemeleon | 0.200 | 0.104 | 0.003 | none |
| scaffold | legacy_sep15_refit | 0.065 | 0.043 | 0.064 | none |
| scaffold | legacy_sep15_refit | 0.100 | 0.039 | 0.068 | none |
| scaffold | legacy_sep15_refit | 0.200 | 0.045 | 0.089 | none |
| scaffold | pool_chemeleon | 0.065 | 0.053 | -0.030 | 0.850 |
| scaffold | pool_chemeleon | 0.100 | 0.053 | -0.035 | 0.750 |
| scaffold | pool_chemeleon | 0.200 | 0.099 | 0.001 | none |
| scaffold | pool_lgbm | 0.065 | 0.114 | 0.035 | none |
| scaffold | pool_lgbm | 0.100 | 0.109 | 0.014 | none |
| scaffold | pool_lgbm | 0.200 | 0.118 | -0.005 | 0.950 |
| scaffold | rank_mean_gates | 0.065 | 0.080 | -0.007 | 0.950 |
| scaffold | rank_mean_gates | 0.100 | 0.085 | -0.013 | 0.950 |
| scaffold | rank_mean_gates | 0.200 | 0.132 | 0.029 | none |
| scaffold | l4_rank_mean | 0.065 | 0.098 | 0.076 | none |
| scaffold | l4_rank_mean | 0.100 | 0.091 | 0.064 | none |
| scaffold | l4_rank_mean | 0.200 | 0.092 | 0.056 | none |
