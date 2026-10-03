# A4 equal-seven OOF stratification and blind reweighting

## Pool

```json
{
  "molecules": 6896,
  "development": 6146,
  "blind": 750,
  "labelled": 3584,
  "labelled_positive": 764,
  "labelled_positive_rate": 0.21316964285714285,
  "labelled_with_direct_pic50": 2335,
  "labelled_assigned_negative": 1249,
  "labelled_with_primary": 1816,
  "labelled_positive_max_log2fc": -0.4828097204,
  "screened_development": 4376,
  "unlabelled_primary_hits": 2207,
  "unlabelled_primary_non_hits": 353,
  "nearest_blind_quartile_edges": [
    0.2625,
    0.29577464788732394,
    0.3389830508474576
  ]
}
```

## Strata (labelled rows, 249/750 call rate on the complete pool)

| stratification | stratum | n | positive_rate | auroc | ap | calls_at_249 | tp_at_249 | precision_at_249 | mean_nearest_blind |
|---|---|---|---|---|---|---|---|---|---|
| log2fc_bin | <= -3 | 399 | 0.599 | 0.625 | 0.728 | 252 | 167 | 0.663 | 0.335 |
| log2fc_bin | (-3, -2] | 325 | 0.532 | 0.650 | 0.656 | 177 | 114 | 0.644 | 0.316 |
| log2fc_bin | (-2, -1] | 329 | 0.505 | 0.595 | 0.586 | 150 | 85 | 0.567 | 0.308 |
| log2fc_bin | (-1, -0.5) | 238 | 0.113 | 0.647 | 0.201 | 74 | 11 | 0.149 | 0.307 |
| log2fc_bin | >= -0.5 | 525 | 0.002 | 0.947 | 0.034 | 49 | 1 | 0.020 | 0.305 |
| log2fc_bin | none | 1768 | 0.089 | 0.853 | 0.415 | 283 | 98 | 0.346 | 0.311 |
| pic50_provenance | direct_pic50 | 2335 | 0.327 | 0.741 | 0.564 | 901 | 476 | 0.528 | 0.319 |
| pic50_provenance | assigned_negative | 1249 | 0.000 | – | – | 84 | 0 | 0.000 | 0.300 |
| nearest_blind_quartile | Q1 | 904 | 0.179 | 0.772 | 0.426 | 231 | 92 | 0.398 | 0.236 |
| nearest_blind_quartile | Q2 | 902 | 0.192 | 0.832 | 0.536 | 236 | 111 | 0.470 | 0.280 |
| nearest_blind_quartile | Q3 | 888 | 0.221 | 0.801 | 0.524 | 262 | 128 | 0.489 | 0.316 |
| nearest_blind_quartile | Q4 | 890 | 0.262 | 0.817 | 0.601 | 256 | 145 | 0.566 | 0.420 |
| predicted_log2fc_bin | <= -3 | 814 | 0.154 | 0.857 | 0.501 | 164 | 78 | 0.476 | 0.319 |
| predicted_log2fc_bin | (-3, -2] | 1015 | 0.294 | 0.816 | 0.616 | 431 | 234 | 0.543 | 0.318 |
| predicted_log2fc_bin | (-2, -1] | 1068 | 0.258 | 0.722 | 0.475 | 362 | 155 | 0.428 | 0.308 |
| predicted_log2fc_bin | (-1, -0.5) | 411 | 0.134 | 0.664 | 0.269 | 27 | 8 | 0.296 | 0.303 |
| predicted_log2fc_bin | >= -0.5 | 276 | 0.036 | 0.591 | 0.187 | 1 | 1 | 1.000 | 0.305 |
| gate_-1 | active | 2558 | 0.249 | 0.808 | 0.557 | 879 | 437 | 0.497 | 0.316 |
| gate_-1 | inactive | 1026 | 0.123 | 0.727 | 0.306 | 106 | 39 | 0.368 | 0.304 |
| gate_-2 | active | 1747 | 0.226 | 0.842 | 0.583 | 553 | 287 | 0.519 | 0.319 |
| gate_-2 | inactive | 1837 | 0.201 | 0.765 | 0.452 | 432 | 189 | 0.438 | 0.306 |

## Reweighted population metrics

| weighting | usable | effective_n | projected_prevalence | auroc | ap | mcc_200 | mcc_225 | mcc_249 | mcc_275 | mcc_300 | mcc_325 | mcc_350 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| none | True | 3584.000 | 0.213 | 0.806 | 0.524 | 0.375 | 0.386 | 0.406 | 0.427 | 0.423 | 0.420 | 0.401 |
| log2fc_bin:nn1 | True | 1152.009 | 0.471 | 0.708 | 0.674 | 0.274 | 0.276 | 0.307 | 0.327 | 0.328 | 0.329 | 0.321 |
| log2fc_bin:nn5 | True | 1988.281 | 0.386 | 0.740 | 0.626 | 0.311 | 0.318 | 0.347 | 0.370 | 0.371 | 0.374 | 0.360 |
| pic50_provenance:nn1 | True | 2367.034 | 0.325 | 0.742 | 0.563 | 0.307 | 0.317 | 0.341 | 0.368 | 0.364 | 0.365 | 0.347 |
| pic50_provenance:nn5 | True | 2841.090 | 0.293 | 0.763 | 0.554 | 0.330 | 0.341 | 0.365 | 0.390 | 0.386 | 0.387 | 0.369 |
| nearest_blind_quartile:nn1 | True | 890.000 | 0.262 | 0.817 | 0.601 | 0.423 | 0.439 | 0.440 | 0.457 | 0.455 | 0.443 | 0.419 |
| nearest_blind_quartile:nn5 | True | 978.444 | 0.260 | 0.816 | 0.597 | 0.420 | 0.436 | 0.439 | 0.455 | 0.454 | 0.442 | 0.419 |
| predicted_log2fc_bin:direct | True | 3528.598 | 0.216 | 0.810 | 0.531 | 0.383 | 0.394 | 0.414 | 0.432 | 0.428 | 0.426 | 0.406 |
| predicted_log2fc_bin:nn1 | True | 2835.578 | 0.231 | 0.819 | 0.553 | 0.400 | 0.413 | 0.433 | 0.445 | 0.441 | 0.439 | 0.419 |
| predicted_log2fc_bin:nn5 | True | 3357.021 | 0.229 | 0.808 | 0.541 | 0.383 | 0.394 | 0.414 | 0.432 | 0.427 | 0.424 | 0.404 |
| gate_-1:direct | True | 3542.475 | 0.219 | 0.807 | 0.530 | 0.377 | 0.389 | 0.409 | 0.428 | 0.425 | 0.423 | 0.404 |
| gate_-1:nn1 | True | 2946.405 | 0.240 | 0.809 | 0.550 | 0.380 | 0.395 | 0.414 | 0.432 | 0.431 | 0.429 | 0.410 |
| gate_-1:nn5 | True | 3376.792 | 0.227 | 0.808 | 0.539 | 0.379 | 0.392 | 0.411 | 0.430 | 0.428 | 0.426 | 0.407 |
| gate_-2:direct | True | 3538.697 | 0.215 | 0.811 | 0.531 | 0.382 | 0.394 | 0.414 | 0.432 | 0.429 | 0.427 | 0.407 |
| gate_-2:nn1 | True | 2897.852 | 0.219 | 0.824 | 0.554 | 0.407 | 0.421 | 0.441 | 0.451 | 0.449 | 0.446 | 0.425 |
| gate_-2:nn5 | True | 3448.199 | 0.216 | 0.814 | 0.536 | 0.388 | 0.400 | 0.420 | 0.437 | 0.434 | 0.431 | 0.411 |

## Best count per weighting

| weighting | count_per_750 | projected_mcc | projected_prevalence |
|---|---|---|---|
| nearest_blind_quartile:nn1 | 275 | 0.457 | 0.262 |
| nearest_blind_quartile:nn5 | 275 | 0.455 | 0.260 |
| gate_-2:nn1 | 275 | 0.451 | 0.219 |
| predicted_log2fc_bin:nn1 | 275 | 0.445 | 0.231 |
| gate_-2:nn5 | 275 | 0.437 | 0.216 |
| gate_-2:direct | 275 | 0.432 | 0.215 |
| predicted_log2fc_bin:direct | 275 | 0.432 | 0.216 |
| predicted_log2fc_bin:nn5 | 275 | 0.432 | 0.229 |
| gate_-1:nn1 | 275 | 0.432 | 0.240 |
| gate_-1:nn5 | 275 | 0.430 | 0.227 |
| gate_-1:direct | 275 | 0.428 | 0.219 |
| none | 275 | 0.427 | 0.213 |
| pic50_provenance:nn5 | 275 | 0.390 | 0.293 |
| log2fc_bin:nn5 | 325 | 0.374 | 0.386 |
| pic50_provenance:nn1 | 275 | 0.368 | 0.325 |
| log2fc_bin:nn1 | 325 | 0.329 | 0.471 |

## Blind gate by submitted call

| submitted_call | n | mean_predicted_log2fc | mean_p_active_-1 | mean_p_active_-2 | mean_nearest_dev_tanimoto |
|---|---|---|---|---|---|
| True | 249 | -2.509 | 0.827 | 0.642 | 0.596 |
| False | 501 | -2.070 | 0.677 | 0.492 | 0.599 |

Submitted-vector check: `{"submitted_positives": 249, "top_reference_matches_submission": true, "overlap": 249}`
