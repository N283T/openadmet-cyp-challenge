# A4 pool and activity-gate variants

## Pool counts

```json
{
  "counts": {
    "blind": 750,
    "borderline": 192,
    "official": 3584,
    "pseudo_negative": 161,
    "unlabelled_hit": 2207,
    "unscreened": 2
  },
  "development": 6146,
  "official_positive": 764,
  "official_positive_rate": 0.21316964285714285,
  "official_with_log2fc": 1816,
  "official_weak_primary": 525,
  "official_weak_primary_positive": 1,
  "official_positive_log2fc_max": -0.4828097204,
  "official_positive_log2fc_ge_gate": 28,
  "pool_n": 3745,
  "pool_positive_rate": 0.2040053404539386,
  "screened_development": 4376,
  "gate_active_development": 2911,
  "primary_source_rows": 17504
}
```

## Ranking (pooled held-fold, official labels)

| split_id | model | auroc | ap | n_official | n_pool |
|---|---|---|---|---|---|
| umap | chemeleon_lr | 0.698 | 0.370 | 3584 | 6146 |
| umap | equal7_filter_0.1 | 0.807 | 0.522 | 3584 | 6146 |
| umap | equal7_filter_0.2 | 0.805 | 0.521 | 3584 | 6146 |
| umap | equal7_filter_0.3 | 0.794 | 0.516 | 3584 | 6146 |
| umap | equal7_x_gate | 0.781 | 0.483 | 3584 | 6146 |
| umap | equal7_x_sqrt_gate | 0.798 | 0.504 | 3584 | 6146 |
| umap | gate_x_tdi | 0.769 | 0.457 | 3584 | 6146 |
| umap | official_lgbm | 0.767 | 0.449 | 3584 | 6146 |
| umap | pool_lgbm | 0.769 | 0.454 | 3584 | 6146 |
| umap | rank_mean_equal7_chemeleon | 0.782 | 0.494 | 3584 | 6146 |
| umap | rank_mean_equal7_official | 0.802 | 0.508 | 3584 | 6146 |
| umap | submitted_equal7 | 0.806 | 0.524 | 3584 | 6146 |
| scaffold | chemeleon_lr | 0.765 | 0.430 | 3584 | 6146 |
| scaffold | gate_x_tdi | 0.809 | 0.513 | 3584 | 6146 |
| scaffold | official_lgbm | 0.808 | 0.509 | 3584 | 6146 |
| scaffold | pool_lgbm | 0.806 | 0.507 | 3584 | 6146 |

## Decision

| split_id | model | nested_mcc | nested_calls | reference_mcc | reference_precision | reference_recall | reference_calls_official | reference_calls_unscreened |
|---|---|---|---|---|---|---|---|---|
| umap | chemeleon_lr | 0.233 | 2124 | 0.237 | 0.365 | 0.499 | 1045 | 2 |
| umap | equal7_filter_0.1 | 0.429 | 2274 | 0.406 | 0.483 | 0.623 | 985 | 2 |
| umap | equal7_filter_0.2 | 0.426 | 2274 | 0.405 | 0.483 | 0.622 | 983 | 2 |
| umap | equal7_filter_0.3 | 0.404 | 2300 | 0.410 | 0.486 | 0.626 | 983 | 2 |
| umap | equal7_x_gate | 0.354 | 1831 | 0.361 | 0.451 | 0.590 | 1001 | 1 |
| umap | equal7_x_sqrt_gate | 0.372 | 2186 | 0.385 | 0.468 | 0.607 | 991 | 2 |
| umap | gate_x_tdi | 0.342 | 2458 | 0.330 | 0.431 | 0.563 | 998 | 2 |
| umap | official_lgbm | 0.312 | 1904 | 0.330 | 0.429 | 0.565 | 1006 | 2 |
| umap | pool_lgbm | 0.322 | 2080 | 0.343 | 0.438 | 0.576 | 1004 | 2 |
| umap | rank_mean_equal7_chemeleon | 0.376 | 2040 | 0.376 | 0.461 | 0.601 | 995 | 2 |
| umap | rank_mean_equal7_official | 0.397 | 2347 | 0.389 | 0.473 | 0.606 | 978 | 2 |
| umap | submitted_equal7 | 0.429 | 2274 | 0.406 | 0.483 | 0.623 | 985 | 2 |
| scaffold | chemeleon_lr | 0.338 | 2460 | 0.320 | 0.419 | 0.565 | 1030 | 2 |
| scaffold | gate_x_tdi | 0.360 | 2117 | 0.389 | 0.473 | 0.607 | 982 | 2 |
| scaffold | official_lgbm | 0.389 | 2423 | 0.387 | 0.469 | 0.611 | 996 | 2 |
| scaffold | pool_lgbm | 0.383 | 2265 | 0.389 | 0.471 | 0.610 | 989 | 2 |

## Official-label MCC by call rate (pooled)

| split_id | model | 25.0% | 29.0% | 33.2% | 37.0% | 40.0% |
|---|---|---|---|---|---|---|
| umap | chemeleon_lr | 0.221 | 0.228 | 0.237 | 0.231 | 0.236 |
| umap | equal7_filter_0.1 | 0.372 | 0.387 | 0.406 | 0.429 | 0.423 |
| umap | equal7_filter_0.2 | 0.370 | 0.385 | 0.405 | 0.426 | 0.420 |
| umap | equal7_filter_0.3 | 0.368 | 0.383 | 0.410 | 0.418 | 0.418 |
| umap | equal7_x_gate | 0.368 | 0.359 | 0.361 | 0.367 | 0.360 |
| umap | equal7_x_sqrt_gate | 0.381 | 0.389 | 0.385 | 0.393 | 0.396 |
| umap | gate_x_tdi | 0.317 | 0.324 | 0.330 | 0.332 | 0.342 |
| umap | official_lgbm | 0.332 | 0.329 | 0.330 | 0.329 | 0.330 |
| umap | pool_lgbm | 0.337 | 0.338 | 0.343 | 0.344 | 0.336 |
| umap | rank_mean_equal7_chemeleon | 0.365 | 0.367 | 0.376 | 0.365 | 0.362 |
| umap | rank_mean_equal7_official | 0.380 | 0.390 | 0.389 | 0.404 | 0.403 |
| umap | submitted_equal7 | 0.370 | 0.387 | 0.406 | 0.429 | 0.423 |
| scaffold | chemeleon_lr | 0.285 | 0.306 | 0.320 | 0.330 | 0.338 |
| scaffold | gate_x_tdi | 0.377 | 0.377 | 0.389 | 0.390 | 0.387 |
| scaffold | official_lgbm | 0.383 | 0.383 | 0.387 | 0.391 | 0.395 |
| scaffold | pool_lgbm | 0.372 | 0.384 | 0.389 | 0.391 | 0.390 |

## Paired group bootstrap

| split_id | model | control | operating_point | delta | ci_low | ci_high | p_delta_gt_0 |
|---|---|---|---|---|---|---|---|
| umap | chemeleon_lr | submitted_equal7 | nested | -0.197 | -0.236 | -0.159 | 0.000 |
| umap | chemeleon_lr | submitted_equal7 | reference | -0.169 | -0.216 | -0.127 | 0.000 |
| umap | equal7_filter_0.1 | submitted_equal7 | nested | 0.000 | -0.001 | 0.001 | 0.333 |
| umap | equal7_filter_0.1 | submitted_equal7 | reference | 0.000 | 0.000 | 0.000 | 0.000 |
| umap | equal7_filter_0.2 | submitted_equal7 | nested | -0.003 | -0.008 | 0.000 | 0.037 |
| umap | equal7_filter_0.2 | submitted_equal7 | reference | -0.001 | -0.004 | 0.003 | 0.364 |
| umap | equal7_filter_0.3 | submitted_equal7 | nested | -0.026 | -0.043 | -0.010 | 0.001 |
| umap | equal7_filter_0.3 | submitted_equal7 | reference | 0.004 | -0.005 | 0.012 | 0.802 |
| umap | equal7_x_gate | submitted_equal7 | nested | -0.075 | -0.105 | -0.043 | 0.000 |
| umap | equal7_x_gate | submitted_equal7 | reference | -0.045 | -0.079 | -0.008 | 0.011 |
| umap | equal7_x_sqrt_gate | submitted_equal7 | nested | -0.057 | -0.083 | -0.032 | 0.000 |
| umap | equal7_x_sqrt_gate | submitted_equal7 | reference | -0.021 | -0.049 | 0.010 | 0.071 |
| umap | gate_x_tdi | submitted_equal7 | nested | -0.087 | -0.120 | -0.056 | 0.000 |
| umap | gate_x_tdi | submitted_equal7 | reference | -0.076 | -0.119 | -0.043 | 0.000 |
| umap | official_lgbm | submitted_equal7 | nested | -0.117 | -0.151 | -0.081 | 0.000 |
| umap | official_lgbm | submitted_equal7 | reference | -0.076 | -0.117 | -0.042 | 0.000 |
| umap | pool_lgbm | submitted_equal7 | nested | -0.108 | -0.144 | -0.075 | 0.000 |
| umap | pool_lgbm | submitted_equal7 | reference | -0.063 | -0.100 | -0.029 | 0.000 |
| umap | rank_mean_equal7_chemeleon | submitted_equal7 | nested | -0.054 | -0.089 | -0.019 | 0.000 |
| umap | rank_mean_equal7_chemeleon | submitted_equal7 | reference | -0.030 | -0.065 | 0.005 | 0.047 |
| umap | rank_mean_equal7_official | submitted_equal7 | nested | -0.032 | -0.055 | -0.010 | 0.001 |
| umap | rank_mean_equal7_official | submitted_equal7 | reference | -0.017 | -0.039 | 0.004 | 0.064 |
| umap | chemeleon_lr | official_lgbm | nested | -0.080 | -0.134 | -0.032 | 0.000 |
| umap | chemeleon_lr | official_lgbm | reference | -0.093 | -0.145 | -0.043 | 0.000 |
| umap | equal7_filter_0.1 | official_lgbm | nested | 0.117 | 0.081 | 0.151 | 1.000 |
| umap | equal7_filter_0.1 | official_lgbm | reference | 0.076 | 0.042 | 0.117 | 1.000 |
| umap | equal7_filter_0.2 | official_lgbm | nested | 0.114 | 0.077 | 0.149 | 1.000 |
| umap | equal7_filter_0.2 | official_lgbm | reference | 0.076 | 0.041 | 0.117 | 1.000 |
| umap | equal7_filter_0.3 | official_lgbm | nested | 0.091 | 0.052 | 0.129 | 1.000 |
| umap | equal7_filter_0.3 | official_lgbm | reference | 0.080 | 0.043 | 0.120 | 1.000 |
| umap | equal7_x_gate | official_lgbm | nested | 0.042 | 0.001 | 0.082 | 0.979 |
| umap | equal7_x_gate | official_lgbm | reference | 0.031 | -0.017 | 0.075 | 0.900 |
| umap | equal7_x_sqrt_gate | official_lgbm | nested | 0.060 | 0.019 | 0.100 | 0.996 |
| umap | equal7_x_sqrt_gate | official_lgbm | reference | 0.055 | 0.011 | 0.099 | 0.996 |
| umap | gate_x_tdi | official_lgbm | nested | 0.030 | 0.002 | 0.055 | 0.984 |
| umap | gate_x_tdi | official_lgbm | reference | 0.000 | -0.020 | 0.019 | 0.500 |
| umap | pool_lgbm | official_lgbm | nested | 0.009 | -0.013 | 0.032 | 0.776 |
| umap | pool_lgbm | official_lgbm | reference | 0.013 | -0.004 | 0.031 | 0.924 |
| umap | rank_mean_equal7_chemeleon | official_lgbm | nested | 0.063 | 0.012 | 0.111 | 0.996 |
| umap | rank_mean_equal7_chemeleon | official_lgbm | reference | 0.046 | -0.003 | 0.094 | 0.961 |
| umap | rank_mean_equal7_official | official_lgbm | nested | 0.085 | 0.056 | 0.114 | 1.000 |
| umap | rank_mean_equal7_official | official_lgbm | reference | 0.059 | 0.033 | 0.088 | 1.000 |
| umap | submitted_equal7 | official_lgbm | nested | 0.117 | 0.081 | 0.151 | 1.000 |
| umap | submitted_equal7 | official_lgbm | reference | 0.076 | 0.042 | 0.117 | 1.000 |
| scaffold | chemeleon_lr | official_lgbm | nested | -0.051 | -0.087 | -0.014 | 0.003 |
| scaffold | chemeleon_lr | official_lgbm | reference | -0.068 | -0.106 | -0.025 | 0.000 |
| scaffold | gate_x_tdi | official_lgbm | nested | -0.030 | -0.053 | -0.008 | 0.003 |
| scaffold | gate_x_tdi | official_lgbm | reference | 0.002 | -0.021 | 0.022 | 0.550 |
| scaffold | pool_lgbm | official_lgbm | nested | -0.006 | -0.024 | 0.012 | 0.235 |
| scaffold | pool_lgbm | official_lgbm | reference | 0.002 | -0.014 | 0.019 | 0.579 |

## Blind diagnostics

| model | overlap_at_249 | added_vs_submitted | mean_p_active | p_active_gt_0_5 | umap_nested_rates | scaffold_nested_rates |
|---|---|---|---|---|---|---|
| chemeleon_lr | 178 | 71 | – | – | 0.332,0.332,0.332,0.400,0.332 | 0.400,0.400,0.400,0.400,0.400 |
| equal7_filter_0.1 | 249 | 0 | 0.727 | 572 | 0.370,0.370,0.370,0.370,0.370 |  |
| equal7_filter_0.2 | 247 | 2 | 0.727 | 572 | 0.370,0.370,0.370,0.370,0.370 |  |
| equal7_filter_0.3 | 243 | 6 | 0.727 | 572 | 0.332,0.370,0.370,0.400,0.400 |  |
| equal7_x_gate | 207 | 42 | 0.727 | 572 | 0.250,0.250,0.250,0.370,0.370 |  |
| equal7_x_sqrt_gate | 216 | 33 | 0.727 | 572 | 0.290,0.400,0.290,0.400,0.400 |  |
| gate_x_tdi | 186 | 63 | 0.727 | 572 | 0.400,0.400,0.400,0.400,0.400 | 0.370,0.370,0.400,0.250,0.332 |
| official_lgbm | 194 | 55 | – | – | 0.250,0.250,0.250,0.400,0.400 | 0.400,0.400,0.400,0.400,0.370 |
| pool_lgbm | 198 | 51 | – | – | 0.370,0.332,0.250,0.370,0.370 | 0.370,0.370,0.370,0.400,0.332 |
| rank_mean_equal7_chemeleon | 212 | 37 | 0.727 | 572 | 0.332,0.332,0.332,0.332,0.332 |  |
| rank_mean_equal7_official | 217 | 32 | 0.727 | 572 | 0.400,0.370,0.370,0.370,0.400 |  |
| submitted_equal7 | 249 | 0 | – | – | 0.370,0.370,0.370,0.370,0.370 |  |

```json
{
  "n": 750,
  "mean": 0.7270091373759097,
  "quantiles": {
    "0.1": 0.21897753672010295,
    "0.25": 0.5251819739194877,
    "0.5": 0.8607023844982238,
    "0.75": 0.9689719406076939,
    "0.9": 0.9893691025500355
  },
  "gt_0_5": 572,
  "lt_0.1": 13,
  "lt_0.2": 68,
  "lt_0.3": 106,
  "development_active_rate": 0.6652193784277879
}
```

## Historical A4 transition concordance

| transition | lb_delta | column | n_removed | n_added | mean_removed | mean_added | difference | permutation_p | sign_agrees |
|---|---|---|---|---|---|---|---|---|---|
| sep12->sep13 | 0.060 | chemeleon_lr | 131 | 2 | 0.425 | 0.647 | 0.222 | 0.243 | yes |
| sep12->sep13 | 0.060 | equal7_filter_0.1 | 131 | 2 | 0.551 | 0.755 | 0.205 | 0.004 | yes |
| sep12->sep13 | 0.060 | equal7_filter_0.2 | 131 | 2 | 0.467 | 0.755 | 0.288 | 0.156 | yes |
| sep12->sep13 | 0.060 | equal7_filter_0.3 | 131 | 2 | 0.383 | 0.755 | 0.372 | 0.117 | yes |
| sep12->sep13 | 0.060 | equal7_x_gate | 131 | 2 | 0.362 | 0.656 | 0.294 | 0.007 | yes |
| sep12->sep13 | 0.060 | equal7_x_sqrt_gate | 131 | 2 | 0.433 | 0.704 | 0.271 | 0.001 | yes |
| sep12->sep13 | 0.060 | gate_x_tdi | 131 | 2 | 0.204 | 0.262 | 0.058 | 0.697 | yes |
| sep12->sep13 | 0.060 | official_lgbm | 131 | 2 | 0.278 | 0.306 | 0.029 | 0.849 | yes |
| sep12->sep13 | 0.060 | pool_lgbm | 131 | 2 | 0.279 | 0.280 | 0.000 | 0.999 | yes |
| sep12->sep13 | 0.060 | rank_mean_equal7_chemeleon | 131 | 2 | 0.542 | 0.743 | 0.201 | 0.034 | yes |
| sep12->sep13 | 0.060 | rank_mean_equal7_official | 131 | 2 | 0.543 | 0.680 | 0.136 | 0.125 | yes |
| sep12->sep13 | 0.060 | submitted_equal7 | 131 | 2 | 0.541 | 0.725 | 0.184 | 0.004 | yes |
| sep12->sep13 | 0.060 | gate_p_active | 131 | 2 | 0.660 | 0.869 | 0.209 | 0.354 | yes |
| sep14->sep15 | 0.080 | chemeleon_lr | 18 | 18 | 0.496 | 0.441 | -0.055 | 0.523 | no |
| sep14->sep15 | 0.080 | equal7_filter_0.1 | 18 | 18 | 0.627 | 0.694 | 0.067 | 0.000 | yes |
| sep14->sep15 | 0.080 | equal7_filter_0.2 | 18 | 18 | 0.627 | 0.638 | 0.011 | 0.999 | yes |
| sep14->sep15 | 0.080 | equal7_filter_0.3 | 18 | 18 | 0.627 | 0.638 | 0.011 | 0.999 | yes |
| sep14->sep15 | 0.080 | equal7_x_gate | 18 | 18 | 0.474 | 0.520 | 0.047 | 0.414 | yes |
| sep14->sep15 | 0.080 | equal7_x_sqrt_gate | 18 | 18 | 0.538 | 0.588 | 0.050 | 0.173 | yes |
| sep14->sep15 | 0.080 | gate_x_tdi | 18 | 18 | 0.336 | 0.244 | -0.091 | 0.127 | no |
| sep14->sep15 | 0.080 | official_lgbm | 18 | 18 | 0.442 | 0.295 | -0.147 | 0.011 | no |
| sep14->sep15 | 0.080 | pool_lgbm | 18 | 18 | 0.448 | 0.284 | -0.164 | 0.008 | no |
| sep14->sep15 | 0.080 | rank_mean_equal7_chemeleon | 18 | 18 | 0.610 | 0.618 | 0.008 | 0.835 | yes |
| sep14->sep15 | 0.080 | rank_mean_equal7_official | 18 | 18 | 0.666 | 0.624 | -0.041 | 0.154 | no |
| sep14->sep15 | 0.080 | submitted_equal7 | 18 | 18 | 0.608 | 0.655 | 0.047 | 0.000 | yes |
| sep14->sep15 | 0.080 | gate_p_active | 18 | 18 | 0.760 | 0.749 | -0.011 | 0.895 | no |
