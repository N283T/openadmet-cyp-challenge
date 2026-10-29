# Induced D6 coverage defect and the pseudo-negative fix

## Label census

| enzyme | labelled | positives | arms | strong_n | strong_positives | borderline_n | borderline_positives | weak_n | weak_positives | missing_n | missing_positives | weak_any_unlabelled_n | weak_any_unlabelled_positives | device_dose | hidden_from_defect | hidden_positives |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| CYP1A2 | 1413 | 23 | defect,fixed,fixed_full,ceiling_screened,ceiling | 1161 | 22 | 122 | 1 | 130 | 0 | 0 | 0 | 1659 | 0 | 1.541 | 252 | 1 |
| CYP2C9 | 1285 | 38 | defect,fixed,fixed_full,ceiling_screened,ceiling | 875 | 35 | 249 | 3 | 161 | 0 | 0 | 0 | 911 | 0 | 1.225 | 410 | 3 |
| CYP3A4 | 3584 | 764 | defect,fixed,fixed_full,ceiling_screened,ceiling,fixed_missing | 1291 | 605 | 178 | 1 | 347 | 0 | 1768 | 158 | 161 | 0 | 0.393 | 2293 | 159 |

## Nested-rate MCC on real complete labels (pooled held folds)

| enzyme | split_id | arm | mcc_defect | mcc_arm | mcc_ceiling_screened | mcc_ceiling_complete | recovered | recovered_fraction_screened |
|---|---|---|---|---|---|---|---|---|
| CYP1A2 | umap | fixed | 0.113 | 0.110 | 0.164 | 0.164 | -0.003 | -0.058 |
| CYP1A2 | umap | fixed_full | 0.113 | 0.060 | 0.164 | 0.164 | -0.053 | -1.039 |
| CYP1A2 | umap | ceiling_screened | 0.113 | 0.164 | 0.164 | 0.164 | 0.051 | 1.000 |
| CYP1A2 | umap | ceiling | 0.113 | 0.164 | 0.164 | 0.164 | 0.051 | 1.000 |
| CYP1A2 | scaffold | fixed | 0.052 | 0.059 | 0.065 | 0.065 | 0.008 | 0.598 |
| CYP1A2 | scaffold | fixed_full | 0.052 | 0.055 | 0.065 | 0.065 | 0.003 | 0.269 |
| CYP1A2 | scaffold | ceiling_screened | 0.052 | 0.065 | 0.065 | 0.065 | 0.013 | 1.000 |
| CYP1A2 | scaffold | ceiling | 0.052 | 0.065 | 0.065 | 0.065 | 0.013 | 1.000 |
| CYP2C9 | umap | fixed | 0.074 | 0.025 | 0.103 | 0.103 | -0.049 | -1.709 |
| CYP2C9 | umap | fixed_full | 0.074 | 0.039 | 0.103 | 0.103 | -0.035 | -1.206 |
| CYP2C9 | umap | ceiling_screened | 0.074 | 0.103 | 0.103 | 0.103 | 0.029 | 1.000 |
| CYP2C9 | umap | ceiling | 0.074 | 0.103 | 0.103 | 0.103 | 0.029 | 1.000 |
| CYP2C9 | scaffold | fixed | 0.052 | 0.067 | 0.048 | 0.048 | 0.015 | – |
| CYP2C9 | scaffold | fixed_full | 0.052 | 0.017 | 0.048 | 0.048 | -0.035 | – |
| CYP2C9 | scaffold | ceiling_screened | 0.052 | 0.048 | 0.048 | 0.048 | -0.004 | – |
| CYP2C9 | scaffold | ceiling | 0.052 | 0.048 | 0.048 | 0.048 | -0.004 | – |
| CYP3A4 | umap | fixed | 0.089 | 0.125 | 0.107 | 0.312 | 0.036 | 2.019 |
| CYP3A4 | umap | fixed_full | 0.089 | 0.132 | 0.107 | 0.312 | 0.043 | 2.388 |
| CYP3A4 | umap | ceiling_screened | 0.089 | 0.107 | 0.107 | 0.312 | 0.018 | 1.000 |
| CYP3A4 | umap | ceiling | 0.089 | 0.312 | 0.107 | 0.312 | 0.223 | 12.539 |
| CYP3A4 | umap | fixed_missing | 0.089 | 0.330 | 0.107 | 0.312 | 0.241 | 13.541 |
| CYP3A4 | scaffold | fixed | 0.095 | 0.125 | 0.124 | 0.395 | 0.031 | 1.042 |
| CYP3A4 | scaffold | fixed_full | 0.095 | 0.131 | 0.124 | 0.395 | 0.037 | 1.246 |
| CYP3A4 | scaffold | ceiling_screened | 0.095 | 0.124 | 0.124 | 0.395 | 0.029 | 1.000 |
| CYP3A4 | scaffold | ceiling | 0.095 | 0.395 | 0.124 | 0.395 | 0.300 | 10.252 |
| CYP3A4 | scaffold | fixed_missing | 0.095 | 0.340 | 0.124 | 0.395 | 0.245 | 8.365 |

## Ranking on real labels

| enzyme | split_id | arm | population | n | positives | auroc | ap |
|---|---|---|---|---|---|---|---|
| CYP1A2 | umap | defect | broad | 1413 | 23 | 0.661 | 0.095 |
| CYP1A2 | umap | defect | strong | 1161 | 22 | 0.681 | 0.107 |
| CYP1A2 | umap | defect | weak | 130 | 0 | – | – |
| CYP1A2 | umap | fixed | broad | 1413 | 23 | 0.630 | 0.090 |
| CYP1A2 | umap | fixed | strong | 1161 | 22 | 0.640 | 0.100 |
| CYP1A2 | umap | fixed | weak | 130 | 0 | – | – |
| CYP1A2 | umap | fixed_full | broad | 1413 | 23 | 0.582 | 0.088 |
| CYP1A2 | umap | fixed_full | strong | 1161 | 22 | 0.568 | 0.095 |
| CYP1A2 | umap | fixed_full | weak | 130 | 0 | – | – |
| CYP1A2 | umap | ceiling_screened | broad | 1413 | 23 | 0.639 | 0.061 |
| CYP1A2 | umap | ceiling_screened | strong | 1161 | 22 | 0.650 | 0.068 |
| CYP1A2 | umap | ceiling_screened | weak | 130 | 0 | – | – |
| CYP1A2 | umap | ceiling | broad | 1413 | 23 | 0.639 | 0.061 |
| CYP1A2 | umap | ceiling | strong | 1161 | 22 | 0.650 | 0.068 |
| CYP1A2 | umap | ceiling | weak | 130 | 0 | – | – |
| CYP1A2 | scaffold | defect | broad | 1413 | 23 | 0.744 | 0.079 |
| CYP1A2 | scaffold | defect | strong | 1161 | 22 | 0.758 | 0.092 |
| CYP1A2 | scaffold | defect | weak | 130 | 0 | – | – |
| CYP1A2 | scaffold | fixed | broad | 1413 | 23 | 0.718 | 0.077 |
| CYP1A2 | scaffold | fixed | strong | 1161 | 22 | 0.728 | 0.087 |
| CYP1A2 | scaffold | fixed | weak | 130 | 0 | – | – |
| CYP1A2 | scaffold | fixed_full | broad | 1413 | 23 | 0.666 | 0.086 |
| CYP1A2 | scaffold | fixed_full | strong | 1161 | 22 | 0.648 | 0.092 |
| CYP1A2 | scaffold | fixed_full | weak | 130 | 0 | – | – |
| CYP1A2 | scaffold | ceiling_screened | broad | 1413 | 23 | 0.697 | 0.082 |
| CYP1A2 | scaffold | ceiling_screened | strong | 1161 | 22 | 0.706 | 0.092 |
| CYP1A2 | scaffold | ceiling_screened | weak | 130 | 0 | – | – |
| CYP1A2 | scaffold | ceiling | broad | 1413 | 23 | 0.697 | 0.082 |
| CYP1A2 | scaffold | ceiling | strong | 1161 | 22 | 0.706 | 0.092 |
| CYP1A2 | scaffold | ceiling | weak | 130 | 0 | – | – |
| CYP2C9 | umap | defect | broad | 1285 | 38 | 0.571 | 0.067 |
| CYP2C9 | umap | defect | strong | 875 | 35 | 0.607 | 0.101 |
| CYP2C9 | umap | defect | weak | 161 | 0 | – | – |
| CYP2C9 | umap | fixed | broad | 1285 | 38 | 0.563 | 0.056 |
| CYP2C9 | umap | fixed | strong | 875 | 35 | 0.578 | 0.094 |
| CYP2C9 | umap | fixed | weak | 161 | 0 | – | – |
| CYP2C9 | umap | fixed_full | broad | 1285 | 38 | 0.581 | 0.058 |
| CYP2C9 | umap | fixed_full | strong | 875 | 35 | 0.559 | 0.074 |
| CYP2C9 | umap | fixed_full | weak | 161 | 0 | – | – |
| CYP2C9 | umap | ceiling_screened | broad | 1285 | 38 | 0.524 | 0.068 |
| CYP2C9 | umap | ceiling_screened | strong | 875 | 35 | 0.520 | 0.092 |
| CYP2C9 | umap | ceiling_screened | weak | 161 | 0 | – | – |
| CYP2C9 | umap | ceiling | broad | 1285 | 38 | 0.524 | 0.068 |
| CYP2C9 | umap | ceiling | strong | 875 | 35 | 0.520 | 0.092 |
| CYP2C9 | umap | ceiling | weak | 161 | 0 | – | – |
| CYP2C9 | scaffold | defect | broad | 1285 | 38 | 0.655 | 0.059 |
| CYP2C9 | scaffold | defect | strong | 875 | 35 | 0.707 | 0.101 |
| CYP2C9 | scaffold | defect | weak | 161 | 0 | – | – |
| CYP2C9 | scaffold | fixed | broad | 1285 | 38 | 0.639 | 0.049 |
| CYP2C9 | scaffold | fixed | strong | 875 | 35 | 0.667 | 0.072 |
| CYP2C9 | scaffold | fixed | weak | 161 | 0 | – | – |
| CYP2C9 | scaffold | fixed_full | broad | 1285 | 38 | 0.603 | 0.051 |
| CYP2C9 | scaffold | fixed_full | strong | 875 | 35 | 0.586 | 0.068 |
| CYP2C9 | scaffold | fixed_full | weak | 161 | 0 | – | – |
| CYP2C9 | scaffold | ceiling_screened | broad | 1285 | 38 | 0.620 | 0.046 |
| CYP2C9 | scaffold | ceiling_screened | strong | 875 | 35 | 0.631 | 0.066 |
| CYP2C9 | scaffold | ceiling_screened | weak | 161 | 0 | – | – |
| CYP2C9 | scaffold | ceiling | broad | 1285 | 38 | 0.620 | 0.046 |
| CYP2C9 | scaffold | ceiling | strong | 875 | 35 | 0.631 | 0.066 |
| CYP2C9 | scaffold | ceiling | weak | 161 | 0 | – | – |
| CYP3A4 | umap | defect | broad | 3584 | 764 | 0.570 | 0.248 |
| CYP3A4 | umap | defect | strong | 1291 | 605 | 0.604 | 0.540 |
| CYP3A4 | umap | defect | weak | 347 | 0 | – | – |
| CYP3A4 | umap | defect | missing | 1768 | 158 | 0.502 | 0.087 |
| CYP3A4 | umap | fixed | broad | 3584 | 764 | 0.596 | 0.258 |
| CYP3A4 | umap | fixed | strong | 1291 | 605 | 0.612 | 0.560 |
| CYP3A4 | umap | fixed | weak | 347 | 0 | – | – |
| CYP3A4 | umap | fixed | missing | 1768 | 158 | 0.490 | 0.088 |
| CYP3A4 | umap | fixed_full | broad | 3584 | 764 | 0.602 | 0.259 |
| CYP3A4 | umap | fixed_full | strong | 1291 | 605 | 0.620 | 0.569 |
| CYP3A4 | umap | fixed_full | weak | 347 | 0 | – | – |
| CYP3A4 | umap | fixed_full | missing | 1768 | 158 | 0.498 | 0.084 |
| CYP3A4 | umap | ceiling_screened | broad | 3584 | 764 | 0.593 | 0.251 |
| CYP3A4 | umap | ceiling_screened | strong | 1291 | 605 | 0.615 | 0.565 |
| CYP3A4 | umap | ceiling_screened | weak | 347 | 0 | – | – |
| CYP3A4 | umap | ceiling_screened | missing | 1768 | 158 | 0.485 | 0.082 |
| CYP3A4 | umap | ceiling | broad | 3584 | 764 | 0.767 | 0.449 |
| CYP3A4 | umap | ceiling | strong | 1291 | 605 | 0.613 | 0.568 |
| CYP3A4 | umap | ceiling | weak | 347 | 0 | – | – |
| CYP3A4 | umap | ceiling | missing | 1768 | 158 | 0.809 | 0.315 |
| CYP3A4 | umap | fixed_missing | broad | 3584 | 764 | 0.759 | 0.443 |
| CYP3A4 | umap | fixed_missing | strong | 1291 | 605 | 0.618 | 0.570 |
| CYP3A4 | umap | fixed_missing | weak | 347 | 0 | – | – |
| CYP3A4 | umap | fixed_missing | missing | 1768 | 158 | 0.781 | 0.296 |
| CYP3A4 | scaffold | defect | broad | 3584 | 764 | 0.594 | 0.270 |
| CYP3A4 | scaffold | defect | strong | 1291 | 605 | 0.645 | 0.600 |
| CYP3A4 | scaffold | defect | weak | 347 | 0 | – | – |
| CYP3A4 | scaffold | defect | missing | 1768 | 158 | 0.504 | 0.086 |
| CYP3A4 | scaffold | fixed | broad | 3584 | 764 | 0.623 | 0.276 |
| CYP3A4 | scaffold | fixed | strong | 1291 | 605 | 0.645 | 0.599 |
| CYP3A4 | scaffold | fixed | weak | 347 | 0 | – | – |
| CYP3A4 | scaffold | fixed | missing | 1768 | 158 | 0.510 | 0.087 |
| CYP3A4 | scaffold | fixed_full | broad | 3584 | 764 | 0.622 | 0.270 |
| CYP3A4 | scaffold | fixed_full | strong | 1291 | 605 | 0.646 | 0.594 |
| CYP3A4 | scaffold | fixed_full | weak | 347 | 0 | – | – |
| CYP3A4 | scaffold | fixed_full | missing | 1768 | 158 | 0.506 | 0.085 |
| CYP3A4 | scaffold | ceiling_screened | broad | 3584 | 764 | 0.616 | 0.268 |
| CYP3A4 | scaffold | ceiling_screened | strong | 1291 | 605 | 0.643 | 0.604 |
| CYP3A4 | scaffold | ceiling_screened | weak | 347 | 0 | – | – |
| CYP3A4 | scaffold | ceiling_screened | missing | 1768 | 158 | 0.493 | 0.082 |
| CYP3A4 | scaffold | ceiling | broad | 3584 | 764 | 0.808 | 0.509 |
| CYP3A4 | scaffold | ceiling | strong | 1291 | 605 | 0.641 | 0.601 |
| CYP3A4 | scaffold | ceiling | weak | 347 | 0 | – | – |
| CYP3A4 | scaffold | ceiling | missing | 1768 | 158 | 0.890 | 0.452 |
| CYP3A4 | scaffold | fixed_missing | broad | 3584 | 764 | 0.781 | 0.470 |
| CYP3A4 | scaffold | fixed_missing | strong | 1291 | 605 | 0.635 | 0.590 |
| CYP3A4 | scaffold | fixed_missing | weak | 347 | 0 | – | – |
| CYP3A4 | scaffold | fixed_missing | missing | 1768 | 158 | 0.785 | 0.252 |

## Paired group bootstrap versus the induced-defect arm (real labels)

| enzyme | split_id | arm | operating_point | population | delta | ci_low | ci_high | p_delta_gt_0 |
|---|---|---|---|---|---|---|---|---|
| CYP1A2 | umap | fixed | nested | broad | -0.003 | -0.059 | 0.073 | 0.572 |
| CYP1A2 | umap | fixed | nested | strong | -0.008 | -0.063 | 0.067 | 0.549 |
| CYP1A2 | umap | fixed | nested | mixture_phi_0.33 | 0.067 | -0.050 | 0.184 | 0.887 |
| CYP1A2 | umap | fixed | nested | mixture_phi_1.00 | -0.008 | -0.063 | 0.067 | 0.549 |
| CYP1A2 | umap | fixed | rate_0.0653333333333333 | broad | -0.006 | -0.017 | 0.003 | 0.090 |
| CYP1A2 | umap | fixed | rate_0.0653333333333333 | strong | -0.007 | -0.022 | 0.004 | 0.117 |
| CYP1A2 | umap | fixed | rate_0.0653333333333333 | mixture_phi_0.33 | -0.001 | -0.019 | 0.015 | 0.434 |
| CYP1A2 | umap | fixed | rate_0.0653333333333333 | mixture_phi_1.00 | -0.007 | -0.022 | 0.004 | 0.125 |
| CYP1A2 | umap | fixed | rate_0.2 | broad | -0.019 | -0.052 | -0.001 | 0.007 |
| CYP1A2 | umap | fixed | rate_0.2 | strong | -0.025 | -0.062 | -0.005 | 0.000 |
| CYP1A2 | umap | fixed | rate_0.2 | mixture_phi_0.33 | -0.003 | -0.025 | 0.010 | 0.410 |
| CYP1A2 | umap | fixed | rate_0.2 | mixture_phi_1.00 | -0.025 | -0.062 | -0.005 | 0.000 |
| CYP1A2 | umap | fixed_full | nested | broad | -0.053 | -0.145 | 0.079 | 0.175 |
| CYP1A2 | umap | fixed_full | nested | strong | -0.069 | -0.176 | 0.081 | 0.155 |
| CYP1A2 | umap | fixed_full | nested | mixture_phi_0.33 | -0.009 | -0.100 | 0.083 | 0.401 |
| CYP1A2 | umap | fixed_full | nested | mixture_phi_1.00 | -0.069 | -0.176 | 0.081 | 0.155 |
| CYP1A2 | umap | fixed_full | rate_0.0653333333333333 | broad | -0.014 | -0.036 | 0.012 | 0.127 |
| CYP1A2 | umap | fixed_full | rate_0.0653333333333333 | strong | -0.024 | -0.053 | 0.004 | 0.044 |
| CYP1A2 | umap | fixed_full | rate_0.0653333333333333 | mixture_phi_0.33 | 0.004 | -0.020 | 0.031 | 0.538 |
| CYP1A2 | umap | fixed_full | rate_0.0653333333333333 | mixture_phi_1.00 | -0.024 | -0.053 | 0.004 | 0.044 |
| CYP1A2 | umap | fixed_full | rate_0.2 | broad | -0.029 | -0.069 | 0.010 | 0.071 |
| CYP1A2 | umap | fixed_full | rate_0.2 | strong | -0.041 | -0.088 | 0.006 | 0.038 |
| CYP1A2 | umap | fixed_full | rate_0.2 | mixture_phi_0.33 | 0.001 | -0.025 | 0.025 | 0.547 |
| CYP1A2 | umap | fixed_full | rate_0.2 | mixture_phi_1.00 | -0.041 | -0.088 | 0.006 | 0.038 |
| CYP1A2 | umap | ceiling_screened | nested | broad | 0.051 | -0.001 | 0.130 | 0.971 |
| CYP1A2 | umap | ceiling_screened | nested | strong | 0.052 | -0.004 | 0.122 | 0.969 |
| CYP1A2 | umap | ceiling_screened | nested | mixture_phi_0.33 | 0.124 | 0.002 | 0.243 | 0.993 |
| CYP1A2 | umap | ceiling_screened | nested | mixture_phi_1.00 | 0.052 | -0.004 | 0.122 | 0.969 |
| CYP1A2 | umap | ceiling_screened | rate_0.0653333333333333 | broad | -0.003 | -0.012 | 0.010 | 0.286 |
| CYP1A2 | umap | ceiling_screened | rate_0.0653333333333333 | strong | -0.007 | -0.020 | 0.006 | 0.118 |
| CYP1A2 | umap | ceiling_screened | rate_0.0653333333333333 | mixture_phi_0.33 | 0.002 | -0.011 | 0.020 | 0.610 |
| CYP1A2 | umap | ceiling_screened | rate_0.0653333333333333 | mixture_phi_1.00 | -0.007 | -0.020 | 0.006 | 0.127 |
| CYP1A2 | umap | ceiling_screened | rate_0.2 | broad | 0.014 | -0.030 | 0.064 | 0.738 |
| CYP1A2 | umap | ceiling_screened | rate_0.2 | strong | 0.010 | -0.040 | 0.067 | 0.643 |
| CYP1A2 | umap | ceiling_screened | rate_0.2 | mixture_phi_0.33 | 0.020 | -0.010 | 0.052 | 0.898 |
| CYP1A2 | umap | ceiling_screened | rate_0.2 | mixture_phi_1.00 | 0.010 | -0.040 | 0.067 | 0.643 |
| CYP1A2 | umap | ceiling | nested | broad | 0.051 | -0.001 | 0.130 | 0.971 |
| CYP1A2 | umap | ceiling | nested | strong | 0.052 | -0.004 | 0.122 | 0.969 |
| CYP1A2 | umap | ceiling | nested | mixture_phi_0.33 | 0.124 | 0.002 | 0.243 | 0.993 |
| CYP1A2 | umap | ceiling | nested | mixture_phi_1.00 | 0.052 | -0.004 | 0.122 | 0.969 |
| CYP1A2 | umap | ceiling | rate_0.0653333333333333 | broad | -0.003 | -0.012 | 0.010 | 0.286 |
| CYP1A2 | umap | ceiling | rate_0.0653333333333333 | strong | -0.007 | -0.020 | 0.006 | 0.118 |
| CYP1A2 | umap | ceiling | rate_0.0653333333333333 | mixture_phi_0.33 | 0.002 | -0.011 | 0.020 | 0.610 |
| CYP1A2 | umap | ceiling | rate_0.0653333333333333 | mixture_phi_1.00 | -0.007 | -0.020 | 0.006 | 0.127 |
| CYP1A2 | umap | ceiling | rate_0.2 | broad | 0.014 | -0.030 | 0.064 | 0.738 |
| CYP1A2 | umap | ceiling | rate_0.2 | strong | 0.010 | -0.040 | 0.067 | 0.643 |
| CYP1A2 | umap | ceiling | rate_0.2 | mixture_phi_0.33 | 0.020 | -0.010 | 0.052 | 0.898 |
| CYP1A2 | umap | ceiling | rate_0.2 | mixture_phi_1.00 | 0.010 | -0.040 | 0.067 | 0.643 |
| CYP1A2 | scaffold | fixed | nested | broad | 0.008 | -0.044 | 0.067 | 0.651 |
| CYP1A2 | scaffold | fixed | nested | strong | 0.009 | -0.054 | 0.072 | 0.638 |
| CYP1A2 | scaffold | fixed | nested | mixture_phi_0.33 | 0.002 | -0.034 | 0.034 | 0.552 |
| CYP1A2 | scaffold | fixed | nested | mixture_phi_1.00 | 0.009 | -0.054 | 0.072 | 0.638 |
| CYP1A2 | scaffold | fixed | rate_0.0653333333333333 | broad | -0.010 | -0.025 | 0.002 | 0.042 |
| CYP1A2 | scaffold | fixed | rate_0.0653333333333333 | strong | -0.015 | -0.035 | -0.001 | 0.011 |
| CYP1A2 | scaffold | fixed | rate_0.0653333333333333 | mixture_phi_0.33 | -0.010 | -0.030 | 0.003 | 0.075 |
| CYP1A2 | scaffold | fixed | rate_0.0653333333333333 | mixture_phi_1.00 | -0.015 | -0.035 | -0.001 | 0.012 |
| CYP1A2 | scaffold | fixed | rate_0.2 | broad | -0.035 | -0.084 | 0.014 | 0.085 |
| CYP1A2 | scaffold | fixed | rate_0.2 | strong | -0.047 | -0.108 | 0.011 | 0.055 |
| CYP1A2 | scaffold | fixed | rate_0.2 | mixture_phi_0.33 | -0.009 | -0.044 | 0.026 | 0.308 |
| CYP1A2 | scaffold | fixed | rate_0.2 | mixture_phi_1.00 | -0.047 | -0.108 | 0.011 | 0.055 |
| CYP1A2 | scaffold | fixed_full | nested | broad | 0.003 | -0.049 | 0.065 | 0.571 |
| CYP1A2 | scaffold | fixed_full | nested | strong | -0.001 | -0.061 | 0.057 | 0.461 |
| CYP1A2 | scaffold | fixed_full | nested | mixture_phi_0.33 | 0.004 | -0.031 | 0.040 | 0.606 |
| CYP1A2 | scaffold | fixed_full | nested | mixture_phi_1.00 | -0.001 | -0.061 | 0.057 | 0.464 |
| CYP1A2 | scaffold | fixed_full | rate_0.0653333333333333 | broad | -0.054 | -0.102 | -0.017 | 0.000 |
| CYP1A2 | scaffold | fixed_full | rate_0.0653333333333333 | strong | -0.067 | -0.126 | -0.024 | 0.000 |
| CYP1A2 | scaffold | fixed_full | rate_0.0653333333333333 | mixture_phi_0.33 | -0.028 | -0.078 | 0.001 | 0.027 |
| CYP1A2 | scaffold | fixed_full | rate_0.0653333333333333 | mixture_phi_1.00 | -0.067 | -0.126 | -0.024 | 0.000 |
| CYP1A2 | scaffold | fixed_full | rate_0.2 | broad | -0.020 | -0.086 | 0.050 | 0.298 |
| CYP1A2 | scaffold | fixed_full | rate_0.2 | strong | -0.034 | -0.106 | 0.041 | 0.190 |
| CYP1A2 | scaffold | fixed_full | rate_0.2 | mixture_phi_0.33 | 0.008 | -0.037 | 0.053 | 0.645 |
| CYP1A2 | scaffold | fixed_full | rate_0.2 | mixture_phi_1.00 | -0.034 | -0.106 | 0.041 | 0.190 |
| CYP1A2 | scaffold | ceiling_screened | nested | broad | 0.013 | -0.041 | 0.074 | 0.686 |
| CYP1A2 | scaffold | ceiling_screened | nested | strong | 0.014 | -0.048 | 0.076 | 0.658 |
| CYP1A2 | scaffold | ceiling_screened | nested | mixture_phi_0.33 | 0.002 | -0.033 | 0.036 | 0.565 |
| CYP1A2 | scaffold | ceiling_screened | nested | mixture_phi_1.00 | 0.014 | -0.048 | 0.076 | 0.658 |
| CYP1A2 | scaffold | ceiling_screened | rate_0.0653333333333333 | broad | -0.012 | -0.101 | 0.075 | 0.372 |
| CYP1A2 | scaffold | ceiling_screened | rate_0.0653333333333333 | strong | -0.021 | -0.125 | 0.073 | 0.355 |
| CYP1A2 | scaffold | ceiling_screened | rate_0.0653333333333333 | mixture_phi_0.33 | 0.008 | -0.058 | 0.082 | 0.622 |
| CYP1A2 | scaffold | ceiling_screened | rate_0.0653333333333333 | mixture_phi_1.00 | -0.021 | -0.125 | 0.073 | 0.355 |
| CYP1A2 | scaffold | ceiling_screened | rate_0.2 | broad | -0.021 | -0.075 | 0.034 | 0.231 |
| CYP1A2 | scaffold | ceiling_screened | rate_0.2 | strong | -0.030 | -0.098 | 0.039 | 0.193 |
| CYP1A2 | scaffold | ceiling_screened | rate_0.2 | mixture_phi_0.33 | 0.002 | -0.037 | 0.044 | 0.531 |
| CYP1A2 | scaffold | ceiling_screened | rate_0.2 | mixture_phi_1.00 | -0.030 | -0.098 | 0.039 | 0.193 |
| CYP1A2 | scaffold | ceiling | nested | broad | 0.013 | -0.041 | 0.074 | 0.686 |
| CYP1A2 | scaffold | ceiling | nested | strong | 0.014 | -0.048 | 0.076 | 0.658 |
| CYP1A2 | scaffold | ceiling | nested | mixture_phi_0.33 | 0.002 | -0.033 | 0.036 | 0.565 |
| CYP1A2 | scaffold | ceiling | nested | mixture_phi_1.00 | 0.014 | -0.048 | 0.076 | 0.658 |
| CYP1A2 | scaffold | ceiling | rate_0.0653333333333333 | broad | -0.012 | -0.101 | 0.075 | 0.372 |
| CYP1A2 | scaffold | ceiling | rate_0.0653333333333333 | strong | -0.021 | -0.125 | 0.073 | 0.355 |
| CYP1A2 | scaffold | ceiling | rate_0.0653333333333333 | mixture_phi_0.33 | 0.008 | -0.058 | 0.082 | 0.622 |
| CYP1A2 | scaffold | ceiling | rate_0.0653333333333333 | mixture_phi_1.00 | -0.021 | -0.125 | 0.073 | 0.355 |
| CYP1A2 | scaffold | ceiling | rate_0.2 | broad | -0.021 | -0.075 | 0.034 | 0.231 |
| CYP1A2 | scaffold | ceiling | rate_0.2 | strong | -0.030 | -0.098 | 0.039 | 0.193 |
| CYP1A2 | scaffold | ceiling | rate_0.2 | mixture_phi_0.33 | 0.002 | -0.037 | 0.044 | 0.531 |
| CYP1A2 | scaffold | ceiling | rate_0.2 | mixture_phi_1.00 | -0.030 | -0.098 | 0.039 | 0.193 |
| CYP2C9 | umap | fixed | nested | broad | -0.049 | -0.097 | -0.009 | 0.000 |
| CYP2C9 | umap | fixed | nested | strong | -0.061 | -0.124 | -0.008 | 0.000 |
| CYP2C9 | umap | fixed | nested | mixture_phi_0.33 | -0.022 | -0.073 | -0.005 | 0.000 |
| CYP2C9 | umap | fixed | nested | mixture_phi_1.00 | -0.061 | -0.124 | -0.008 | 0.000 |
| CYP2C9 | umap | fixed | rate_0.0653333333333333 | broad | -0.007 | -0.017 | 0.002 | 0.042 |
| CYP2C9 | umap | fixed | rate_0.0653333333333333 | strong | -0.015 | -0.036 | 0.002 | 0.032 |
| CYP2C9 | umap | fixed | rate_0.0653333333333333 | mixture_phi_0.33 | -0.002 | -0.012 | 0.006 | 0.326 |
| CYP2C9 | umap | fixed | rate_0.0653333333333333 | mixture_phi_1.00 | -0.015 | -0.036 | 0.002 | 0.035 |
| CYP2C9 | umap | fixed | rate_0.2 | broad | -0.025 | -0.066 | 0.020 | 0.119 |
| CYP2C9 | umap | fixed | rate_0.2 | strong | -0.044 | -0.099 | 0.012 | 0.070 |
| CYP2C9 | umap | fixed | rate_0.2 | mixture_phi_0.33 | -0.006 | -0.036 | 0.024 | 0.322 |
| CYP2C9 | umap | fixed | rate_0.2 | mixture_phi_1.00 | -0.044 | -0.099 | 0.012 | 0.070 |
| CYP2C9 | umap | fixed_full | nested | broad | -0.035 | -0.112 | 0.071 | 0.222 |
| CYP2C9 | umap | fixed_full | nested | strong | -0.059 | -0.152 | 0.078 | 0.192 |
| CYP2C9 | umap | fixed_full | nested | mixture_phi_0.33 | 0.001 | -0.060 | 0.070 | 0.459 |
| CYP2C9 | umap | fixed_full | nested | mixture_phi_1.00 | -0.059 | -0.152 | 0.078 | 0.192 |
| CYP2C9 | umap | fixed_full | rate_0.0653333333333333 | broad | 0.029 | -0.007 | 0.074 | 0.906 |
| CYP2C9 | umap | fixed_full | rate_0.0653333333333333 | strong | 0.016 | -0.038 | 0.083 | 0.664 |
| CYP2C9 | umap | fixed_full | rate_0.0653333333333333 | mixture_phi_0.33 | 0.045 | 0.011 | 0.102 | 0.998 |
| CYP2C9 | umap | fixed_full | rate_0.0653333333333333 | mixture_phi_1.00 | 0.016 | -0.038 | 0.083 | 0.664 |
| CYP2C9 | umap | fixed_full | rate_0.2 | broad | 0.000 | -0.066 | 0.075 | 0.514 |
| CYP2C9 | umap | fixed_full | rate_0.2 | strong | -0.038 | -0.126 | 0.056 | 0.205 |
| CYP2C9 | umap | fixed_full | rate_0.2 | mixture_phi_0.33 | 0.018 | -0.027 | 0.077 | 0.765 |
| CYP2C9 | umap | fixed_full | rate_0.2 | mixture_phi_1.00 | -0.038 | -0.126 | 0.056 | 0.205 |
| CYP2C9 | umap | ceiling_screened | nested | broad | 0.029 | -0.016 | 0.094 | 0.838 |
| CYP2C9 | umap | ceiling_screened | nested | strong | 0.036 | -0.026 | 0.124 | 0.801 |
| CYP2C9 | umap | ceiling_screened | nested | mixture_phi_0.33 | 0.023 | -0.001 | 0.086 | 0.951 |
| CYP2C9 | umap | ceiling_screened | nested | mixture_phi_1.00 | 0.036 | -0.026 | 0.124 | 0.809 |
| CYP2C9 | umap | ceiling_screened | rate_0.0653333333333333 | broad | -0.005 | -0.015 | 0.005 | 0.154 |
| CYP2C9 | umap | ceiling_screened | rate_0.0653333333333333 | strong | -0.017 | -0.040 | -0.001 | 0.013 |
| CYP2C9 | umap | ceiling_screened | rate_0.0653333333333333 | mixture_phi_0.33 | 0.004 | -0.009 | 0.020 | 0.708 |
| CYP2C9 | umap | ceiling_screened | rate_0.0653333333333333 | mixture_phi_1.00 | -0.017 | -0.040 | -0.001 | 0.015 |
| CYP2C9 | umap | ceiling_screened | rate_0.2 | broad | -0.004 | -0.054 | 0.054 | 0.419 |
| CYP2C9 | umap | ceiling_screened | rate_0.2 | strong | -0.023 | -0.091 | 0.054 | 0.237 |
| CYP2C9 | umap | ceiling_screened | rate_0.2 | mixture_phi_0.33 | 0.014 | -0.019 | 0.053 | 0.775 |
| CYP2C9 | umap | ceiling_screened | rate_0.2 | mixture_phi_1.00 | -0.023 | -0.091 | 0.054 | 0.237 |
| CYP2C9 | umap | ceiling | nested | broad | 0.029 | -0.016 | 0.094 | 0.838 |
| CYP2C9 | umap | ceiling | nested | strong | 0.036 | -0.026 | 0.124 | 0.801 |
| CYP2C9 | umap | ceiling | nested | mixture_phi_0.33 | 0.023 | -0.001 | 0.086 | 0.951 |
| CYP2C9 | umap | ceiling | nested | mixture_phi_1.00 | 0.036 | -0.026 | 0.124 | 0.809 |
| CYP2C9 | umap | ceiling | rate_0.0653333333333333 | broad | -0.005 | -0.015 | 0.005 | 0.154 |
| CYP2C9 | umap | ceiling | rate_0.0653333333333333 | strong | -0.017 | -0.040 | -0.001 | 0.013 |
| CYP2C9 | umap | ceiling | rate_0.0653333333333333 | mixture_phi_0.33 | 0.004 | -0.009 | 0.020 | 0.708 |
| CYP2C9 | umap | ceiling | rate_0.0653333333333333 | mixture_phi_1.00 | -0.017 | -0.040 | -0.001 | 0.015 |
| CYP2C9 | umap | ceiling | rate_0.2 | broad | -0.004 | -0.054 | 0.054 | 0.419 |
| CYP2C9 | umap | ceiling | rate_0.2 | strong | -0.023 | -0.091 | 0.054 | 0.237 |
| CYP2C9 | umap | ceiling | rate_0.2 | mixture_phi_0.33 | 0.014 | -0.019 | 0.053 | 0.775 |
| CYP2C9 | umap | ceiling | rate_0.2 | mixture_phi_1.00 | -0.023 | -0.091 | 0.054 | 0.237 |
| CYP2C9 | scaffold | fixed | nested | broad | 0.015 | -0.062 | 0.089 | 0.645 |
| CYP2C9 | scaffold | fixed | nested | strong | 0.028 | -0.072 | 0.120 | 0.699 |
| CYP2C9 | scaffold | fixed | nested | mixture_phi_0.33 | 0.026 | -0.031 | 0.080 | 0.823 |
| CYP2C9 | scaffold | fixed | nested | mixture_phi_1.00 | 0.028 | -0.072 | 0.120 | 0.699 |
| CYP2C9 | scaffold | fixed | rate_0.0653333333333333 | broad | -0.051 | -0.112 | -0.009 | 0.000 |
| CYP2C9 | scaffold | fixed | rate_0.0653333333333333 | strong | -0.079 | -0.158 | -0.019 | 0.000 |
| CYP2C9 | scaffold | fixed | rate_0.0653333333333333 | mixture_phi_0.33 | -0.029 | -0.076 | 0.003 | 0.053 |
| CYP2C9 | scaffold | fixed | rate_0.0653333333333333 | mixture_phi_1.00 | -0.079 | -0.158 | -0.019 | 0.000 |
| CYP2C9 | scaffold | fixed | rate_0.2 | broad | -0.014 | -0.056 | 0.022 | 0.221 |
| CYP2C9 | scaffold | fixed | rate_0.2 | strong | -0.029 | -0.081 | 0.023 | 0.128 |
| CYP2C9 | scaffold | fixed | rate_0.2 | mixture_phi_0.33 | -0.001 | -0.032 | 0.030 | 0.506 |
| CYP2C9 | scaffold | fixed | rate_0.2 | mixture_phi_1.00 | -0.029 | -0.081 | 0.023 | 0.128 |
| CYP2C9 | scaffold | fixed_full | nested | broad | -0.035 | -0.119 | 0.056 | 0.202 |
| CYP2C9 | scaffold | fixed_full | nested | strong | -0.054 | -0.165 | 0.050 | 0.152 |
| CYP2C9 | scaffold | fixed_full | nested | mixture_phi_0.33 | -0.012 | -0.071 | 0.063 | 0.365 |
| CYP2C9 | scaffold | fixed_full | nested | mixture_phi_1.00 | -0.054 | -0.165 | 0.050 | 0.152 |
| CYP2C9 | scaffold | fixed_full | rate_0.0653333333333333 | broad | -0.013 | -0.085 | 0.056 | 0.351 |
| CYP2C9 | scaffold | fixed_full | rate_0.0653333333333333 | strong | -0.042 | -0.143 | 0.045 | 0.180 |
| CYP2C9 | scaffold | fixed_full | rate_0.0653333333333333 | mixture_phi_0.33 | 0.012 | -0.049 | 0.067 | 0.649 |
| CYP2C9 | scaffold | fixed_full | rate_0.0653333333333333 | mixture_phi_1.00 | -0.042 | -0.143 | 0.045 | 0.180 |
| CYP2C9 | scaffold | fixed_full | rate_0.2 | broad | -0.001 | -0.058 | 0.059 | 0.494 |
| CYP2C9 | scaffold | fixed_full | rate_0.2 | strong | -0.059 | -0.122 | 0.006 | 0.046 |
| CYP2C9 | scaffold | fixed_full | rate_0.2 | mixture_phi_0.33 | 0.016 | -0.022 | 0.056 | 0.789 |
| CYP2C9 | scaffold | fixed_full | rate_0.2 | mixture_phi_1.00 | -0.059 | -0.122 | 0.006 | 0.046 |
| CYP2C9 | scaffold | ceiling_screened | nested | broad | -0.004 | -0.093 | 0.084 | 0.458 |
| CYP2C9 | scaffold | ceiling_screened | nested | strong | -0.008 | -0.117 | 0.094 | 0.448 |
| CYP2C9 | scaffold | ceiling_screened | nested | mixture_phi_0.33 | 0.024 | -0.039 | 0.083 | 0.751 |
| CYP2C9 | scaffold | ceiling_screened | nested | mixture_phi_1.00 | -0.008 | -0.117 | 0.094 | 0.448 |
| CYP2C9 | scaffold | ceiling_screened | rate_0.0653333333333333 | broad | -0.046 | -0.105 | -0.003 | 0.010 |
| CYP2C9 | scaffold | ceiling_screened | rate_0.0653333333333333 | strong | -0.076 | -0.157 | -0.014 | 0.001 |
| CYP2C9 | scaffold | ceiling_screened | rate_0.0653333333333333 | mixture_phi_0.33 | -0.014 | -0.057 | 0.020 | 0.257 |
| CYP2C9 | scaffold | ceiling_screened | rate_0.0653333333333333 | mixture_phi_1.00 | -0.076 | -0.157 | -0.014 | 0.001 |
| CYP2C9 | scaffold | ceiling_screened | rate_0.2 | broad | -0.027 | -0.071 | 0.014 | 0.104 |
| CYP2C9 | scaffold | ceiling_screened | rate_0.2 | strong | -0.056 | -0.118 | -0.003 | 0.020 |
| CYP2C9 | scaffold | ceiling_screened | rate_0.2 | mixture_phi_0.33 | 0.002 | -0.035 | 0.039 | 0.557 |
| CYP2C9 | scaffold | ceiling_screened | rate_0.2 | mixture_phi_1.00 | -0.056 | -0.118 | -0.003 | 0.020 |
| CYP2C9 | scaffold | ceiling | nested | broad | -0.004 | -0.093 | 0.084 | 0.458 |
| CYP2C9 | scaffold | ceiling | nested | strong | -0.008 | -0.117 | 0.094 | 0.448 |
| CYP2C9 | scaffold | ceiling | nested | mixture_phi_0.33 | 0.024 | -0.039 | 0.083 | 0.751 |
| CYP2C9 | scaffold | ceiling | nested | mixture_phi_1.00 | -0.008 | -0.117 | 0.094 | 0.448 |
| CYP2C9 | scaffold | ceiling | rate_0.0653333333333333 | broad | -0.046 | -0.105 | -0.003 | 0.010 |
| CYP2C9 | scaffold | ceiling | rate_0.0653333333333333 | strong | -0.076 | -0.157 | -0.014 | 0.001 |
| CYP2C9 | scaffold | ceiling | rate_0.0653333333333333 | mixture_phi_0.33 | -0.014 | -0.057 | 0.020 | 0.257 |
| CYP2C9 | scaffold | ceiling | rate_0.0653333333333333 | mixture_phi_1.00 | -0.076 | -0.157 | -0.014 | 0.001 |
| CYP2C9 | scaffold | ceiling | rate_0.2 | broad | -0.027 | -0.071 | 0.014 | 0.104 |
| CYP2C9 | scaffold | ceiling | rate_0.2 | strong | -0.056 | -0.118 | -0.003 | 0.020 |
| CYP2C9 | scaffold | ceiling | rate_0.2 | mixture_phi_0.33 | 0.002 | -0.035 | 0.039 | 0.557 |
| CYP2C9 | scaffold | ceiling | rate_0.2 | mixture_phi_1.00 | -0.056 | -0.118 | -0.003 | 0.020 |
| CYP3A4 | umap | fixed | nested | broad | 0.036 | 0.002 | 0.067 | 0.983 |
| CYP3A4 | umap | fixed | nested | strong | 0.041 | 0.009 | 0.073 | 0.993 |
| CYP3A4 | umap | fixed | nested | mixture_phi_0.33 | 0.200 | 0.161 | 0.236 | 1.000 |
| CYP3A4 | umap | fixed | nested | mixture_phi_1.00 | 0.041 | 0.009 | 0.073 | 0.993 |
| CYP3A4 | umap | fixed | rate_0.0653333333333333 | broad | -0.002 | -0.041 | 0.035 | 0.462 |
| CYP3A4 | umap | fixed | rate_0.0653333333333333 | strong | 0.040 | -0.012 | 0.094 | 0.928 |
| CYP3A4 | umap | fixed | rate_0.0653333333333333 | mixture_phi_0.33 | 0.077 | 0.017 | 0.135 | 0.991 |
| CYP3A4 | umap | fixed | rate_0.0653333333333333 | mixture_phi_1.00 | 0.040 | -0.012 | 0.094 | 0.929 |
| CYP3A4 | umap | fixed | rate_0.2 | broad | 0.005 | -0.033 | 0.041 | 0.597 |
| CYP3A4 | umap | fixed | rate_0.2 | strong | 0.018 | -0.037 | 0.072 | 0.731 |
| CYP3A4 | umap | fixed | rate_0.2 | mixture_phi_0.33 | 0.135 | 0.080 | 0.191 | 1.000 |
| CYP3A4 | umap | fixed | rate_0.2 | mixture_phi_1.00 | 0.018 | -0.037 | 0.072 | 0.731 |
| CYP3A4 | umap | fixed_full | nested | broad | 0.043 | 0.010 | 0.074 | 0.994 |
| CYP3A4 | umap | fixed_full | nested | strong | 0.047 | 0.008 | 0.084 | 0.988 |
| CYP3A4 | umap | fixed_full | nested | mixture_phi_0.33 | 0.220 | 0.172 | 0.262 | 1.000 |
| CYP3A4 | umap | fixed_full | nested | mixture_phi_1.00 | 0.047 | 0.008 | 0.084 | 0.988 |
| CYP3A4 | umap | fixed_full | rate_0.0653333333333333 | broad | -0.001 | -0.042 | 0.042 | 0.464 |
| CYP3A4 | umap | fixed_full | rate_0.0653333333333333 | strong | 0.057 | 0.003 | 0.111 | 0.978 |
| CYP3A4 | umap | fixed_full | rate_0.0653333333333333 | mixture_phi_0.33 | 0.105 | 0.036 | 0.163 | 0.999 |
| CYP3A4 | umap | fixed_full | rate_0.0653333333333333 | mixture_phi_1.00 | 0.057 | 0.003 | 0.111 | 0.979 |
| CYP3A4 | umap | fixed_full | rate_0.2 | broad | 0.019 | -0.027 | 0.065 | 0.762 |
| CYP3A4 | umap | fixed_full | rate_0.2 | strong | 0.044 | -0.010 | 0.102 | 0.938 |
| CYP3A4 | umap | fixed_full | rate_0.2 | mixture_phi_0.33 | 0.167 | 0.107 | 0.225 | 1.000 |
| CYP3A4 | umap | fixed_full | rate_0.2 | mixture_phi_1.00 | 0.044 | -0.010 | 0.102 | 0.938 |
| CYP3A4 | umap | ceiling_screened | nested | broad | 0.018 | -0.017 | 0.054 | 0.821 |
| CYP3A4 | umap | ceiling_screened | nested | strong | 0.036 | -0.008 | 0.080 | 0.946 |
| CYP3A4 | umap | ceiling_screened | nested | mixture_phi_0.33 | 0.202 | 0.153 | 0.247 | 1.000 |
| CYP3A4 | umap | ceiling_screened | nested | mixture_phi_1.00 | 0.036 | -0.008 | 0.080 | 0.946 |
| CYP3A4 | umap | ceiling_screened | rate_0.0653333333333333 | broad | -0.024 | -0.059 | 0.011 | 0.092 |
| CYP3A4 | umap | ceiling_screened | rate_0.0653333333333333 | strong | 0.021 | -0.030 | 0.070 | 0.803 |
| CYP3A4 | umap | ceiling_screened | rate_0.0653333333333333 | mixture_phi_0.33 | 0.069 | 0.009 | 0.122 | 0.988 |
| CYP3A4 | umap | ceiling_screened | rate_0.0653333333333333 | mixture_phi_1.00 | 0.021 | -0.030 | 0.070 | 0.804 |
| CYP3A4 | umap | ceiling_screened | rate_0.2 | broad | 0.009 | -0.032 | 0.056 | 0.654 |
| CYP3A4 | umap | ceiling_screened | rate_0.2 | strong | 0.049 | -0.017 | 0.119 | 0.922 |
| CYP3A4 | umap | ceiling_screened | rate_0.2 | mixture_phi_0.33 | 0.165 | 0.101 | 0.228 | 1.000 |
| CYP3A4 | umap | ceiling_screened | rate_0.2 | mixture_phi_1.00 | 0.049 | -0.017 | 0.119 | 0.923 |
| CYP3A4 | umap | ceiling | nested | broad | 0.223 | 0.156 | 0.286 | 1.000 |
| CYP3A4 | umap | ceiling | nested | strong | -0.026 | -0.092 | 0.038 | 0.219 |
| CYP3A4 | umap | ceiling | nested | mixture_phi_0.33 | 0.152 | 0.104 | 0.194 | 1.000 |
| CYP3A4 | umap | ceiling | nested | mixture_phi_1.00 | -0.026 | -0.092 | 0.038 | 0.219 |
| CYP3A4 | umap | ceiling | rate_0.0653333333333333 | broad | 0.179 | 0.106 | 0.250 | 1.000 |
| CYP3A4 | umap | ceiling | rate_0.0653333333333333 | strong | 0.090 | 0.014 | 0.165 | 0.989 |
| CYP3A4 | umap | ceiling | rate_0.0653333333333333 | mixture_phi_0.33 | 0.153 | 0.072 | 0.228 | 1.000 |
| CYP3A4 | umap | ceiling | rate_0.0653333333333333 | mixture_phi_1.00 | 0.090 | 0.014 | 0.165 | 0.989 |
| CYP3A4 | umap | ceiling | rate_0.2 | broad | 0.247 | 0.188 | 0.300 | 1.000 |
| CYP3A4 | umap | ceiling | rate_0.2 | strong | 0.067 | 0.007 | 0.127 | 0.992 |
| CYP3A4 | umap | ceiling | rate_0.2 | mixture_phi_0.33 | 0.154 | 0.103 | 0.197 | 1.000 |
| CYP3A4 | umap | ceiling | rate_0.2 | mixture_phi_1.00 | 0.067 | 0.007 | 0.127 | 0.992 |
| CYP3A4 | umap | fixed_missing | nested | broad | 0.241 | 0.178 | 0.300 | 1.000 |
| CYP3A4 | umap | fixed_missing | nested | strong | 0.031 | -0.019 | 0.086 | 0.877 |
| CYP3A4 | umap | fixed_missing | nested | mixture_phi_0.33 | 0.139 | 0.099 | 0.177 | 1.000 |
| CYP3A4 | umap | fixed_missing | nested | mixture_phi_1.00 | 0.031 | -0.019 | 0.086 | 0.877 |
| CYP3A4 | umap | fixed_missing | rate_0.0653333333333333 | broad | 0.172 | 0.107 | 0.230 | 1.000 |
| CYP3A4 | umap | fixed_missing | rate_0.0653333333333333 | strong | 0.057 | -0.004 | 0.122 | 0.963 |
| CYP3A4 | umap | fixed_missing | rate_0.0653333333333333 | mixture_phi_0.33 | 0.126 | 0.044 | 0.197 | 1.000 |
| CYP3A4 | umap | fixed_missing | rate_0.0653333333333333 | mixture_phi_1.00 | 0.057 | -0.004 | 0.122 | 0.963 |
| CYP3A4 | umap | fixed_missing | rate_0.2 | broad | 0.252 | 0.191 | 0.306 | 1.000 |
| CYP3A4 | umap | fixed_missing | rate_0.2 | strong | 0.084 | 0.022 | 0.152 | 0.995 |
| CYP3A4 | umap | fixed_missing | rate_0.2 | mixture_phi_0.33 | 0.173 | 0.114 | 0.225 | 1.000 |
| CYP3A4 | umap | fixed_missing | rate_0.2 | mixture_phi_1.00 | 0.084 | 0.022 | 0.152 | 0.995 |
| CYP3A4 | scaffold | fixed | nested | broad | 0.031 | -0.000 | 0.062 | 0.971 |
| CYP3A4 | scaffold | fixed | nested | strong | -0.008 | -0.054 | 0.038 | 0.368 |
| CYP3A4 | scaffold | fixed | nested | mixture_phi_0.33 | 0.138 | 0.091 | 0.184 | 1.000 |
| CYP3A4 | scaffold | fixed | nested | mixture_phi_1.00 | -0.008 | -0.054 | 0.038 | 0.368 |
| CYP3A4 | scaffold | fixed | rate_0.0653333333333333 | broad | -0.011 | -0.046 | 0.024 | 0.256 |
| CYP3A4 | scaffold | fixed | rate_0.0653333333333333 | strong | -0.004 | -0.050 | 0.045 | 0.452 |
| CYP3A4 | scaffold | fixed | rate_0.0653333333333333 | mixture_phi_0.33 | 0.063 | 0.009 | 0.118 | 0.989 |
| CYP3A4 | scaffold | fixed | rate_0.0653333333333333 | mixture_phi_1.00 | -0.004 | -0.050 | 0.045 | 0.452 |
| CYP3A4 | scaffold | fixed | rate_0.2 | broad | 0.002 | -0.027 | 0.032 | 0.538 |
| CYP3A4 | scaffold | fixed | rate_0.2 | strong | -0.016 | -0.060 | 0.031 | 0.234 |
| CYP3A4 | scaffold | fixed | rate_0.2 | mixture_phi_0.33 | 0.081 | 0.030 | 0.129 | 1.000 |
| CYP3A4 | scaffold | fixed | rate_0.2 | mixture_phi_1.00 | -0.016 | -0.060 | 0.031 | 0.234 |
| CYP3A4 | scaffold | fixed_full | nested | broad | 0.037 | 0.003 | 0.069 | 0.985 |
| CYP3A4 | scaffold | fixed_full | nested | strong | 0.012 | -0.038 | 0.063 | 0.681 |
| CYP3A4 | scaffold | fixed_full | nested | mixture_phi_0.33 | 0.185 | 0.136 | 0.232 | 1.000 |
| CYP3A4 | scaffold | fixed_full | nested | mixture_phi_1.00 | 0.012 | -0.038 | 0.063 | 0.681 |
| CYP3A4 | scaffold | fixed_full | rate_0.0653333333333333 | broad | -0.016 | -0.051 | 0.020 | 0.173 |
| CYP3A4 | scaffold | fixed_full | rate_0.0653333333333333 | strong | 0.002 | -0.050 | 0.054 | 0.519 |
| CYP3A4 | scaffold | fixed_full | rate_0.0653333333333333 | mixture_phi_0.33 | 0.078 | 0.016 | 0.132 | 0.993 |
| CYP3A4 | scaffold | fixed_full | rate_0.0653333333333333 | mixture_phi_1.00 | 0.002 | -0.050 | 0.054 | 0.520 |
| CYP3A4 | scaffold | fixed_full | rate_0.2 | broad | -0.015 | -0.045 | 0.019 | 0.180 |
| CYP3A4 | scaffold | fixed_full | rate_0.2 | strong | -0.028 | -0.074 | 0.021 | 0.123 |
| CYP3A4 | scaffold | fixed_full | rate_0.2 | mixture_phi_0.33 | 0.091 | 0.040 | 0.142 | 0.999 |
| CYP3A4 | scaffold | fixed_full | rate_0.2 | mixture_phi_1.00 | -0.028 | -0.074 | 0.021 | 0.123 |
| CYP3A4 | scaffold | ceiling_screened | nested | broad | 0.029 | -0.007 | 0.062 | 0.950 |
| CYP3A4 | scaffold | ceiling_screened | nested | strong | 0.007 | -0.048 | 0.059 | 0.604 |
| CYP3A4 | scaffold | ceiling_screened | nested | mixture_phi_0.33 | 0.178 | 0.130 | 0.227 | 1.000 |
| CYP3A4 | scaffold | ceiling_screened | nested | mixture_phi_1.00 | 0.007 | -0.048 | 0.059 | 0.604 |
| CYP3A4 | scaffold | ceiling_screened | rate_0.0653333333333333 | broad | -0.023 | -0.059 | 0.013 | 0.097 |
| CYP3A4 | scaffold | ceiling_screened | rate_0.0653333333333333 | strong | -0.002 | -0.053 | 0.051 | 0.473 |
| CYP3A4 | scaffold | ceiling_screened | rate_0.0653333333333333 | mixture_phi_0.33 | 0.052 | -0.005 | 0.107 | 0.963 |
| CYP3A4 | scaffold | ceiling_screened | rate_0.0653333333333333 | mixture_phi_1.00 | -0.002 | -0.053 | 0.051 | 0.473 |
| CYP3A4 | scaffold | ceiling_screened | rate_0.2 | broad | -0.023 | -0.055 | 0.009 | 0.075 |
| CYP3A4 | scaffold | ceiling_screened | rate_0.2 | strong | -0.018 | -0.067 | 0.030 | 0.233 |
| CYP3A4 | scaffold | ceiling_screened | rate_0.2 | mixture_phi_0.33 | 0.082 | 0.032 | 0.136 | 0.998 |
| CYP3A4 | scaffold | ceiling_screened | rate_0.2 | mixture_phi_1.00 | -0.018 | -0.067 | 0.030 | 0.233 |
| CYP3A4 | scaffold | ceiling | nested | broad | 0.300 | 0.260 | 0.342 | 1.000 |
| CYP3A4 | scaffold | ceiling | nested | strong | 0.008 | -0.051 | 0.068 | 0.608 |
| CYP3A4 | scaffold | ceiling | nested | mixture_phi_0.33 | 0.159 | 0.111 | 0.209 | 1.000 |
| CYP3A4 | scaffold | ceiling | nested | mixture_phi_1.00 | 0.008 | -0.051 | 0.068 | 0.608 |
| CYP3A4 | scaffold | ceiling | rate_0.0653333333333333 | broad | 0.191 | 0.143 | 0.238 | 1.000 |
| CYP3A4 | scaffold | ceiling | rate_0.0653333333333333 | strong | 0.025 | -0.032 | 0.090 | 0.765 |
| CYP3A4 | scaffold | ceiling | rate_0.0653333333333333 | mixture_phi_0.33 | 0.128 | 0.063 | 0.190 | 1.000 |
| CYP3A4 | scaffold | ceiling | rate_0.0653333333333333 | mixture_phi_1.00 | 0.025 | -0.032 | 0.090 | 0.765 |
| CYP3A4 | scaffold | ceiling | rate_0.2 | broad | 0.249 | 0.206 | 0.294 | 1.000 |
| CYP3A4 | scaffold | ceiling | rate_0.2 | strong | -0.009 | -0.069 | 0.050 | 0.369 |
| CYP3A4 | scaffold | ceiling | rate_0.2 | mixture_phi_0.33 | 0.097 | 0.045 | 0.154 | 1.000 |
| CYP3A4 | scaffold | ceiling | rate_0.2 | mixture_phi_1.00 | -0.009 | -0.069 | 0.050 | 0.369 |
| CYP3A4 | scaffold | fixed_missing | nested | broad | 0.245 | 0.201 | 0.287 | 1.000 |
| CYP3A4 | scaffold | fixed_missing | nested | strong | -0.004 | -0.063 | 0.054 | 0.457 |
| CYP3A4 | scaffold | fixed_missing | nested | mixture_phi_0.33 | 0.092 | 0.040 | 0.146 | 0.999 |
| CYP3A4 | scaffold | fixed_missing | nested | mixture_phi_1.00 | -0.004 | -0.063 | 0.054 | 0.457 |
| CYP3A4 | scaffold | fixed_missing | rate_0.0653333333333333 | broad | 0.147 | 0.101 | 0.195 | 1.000 |
| CYP3A4 | scaffold | fixed_missing | rate_0.0653333333333333 | strong | 0.009 | -0.052 | 0.073 | 0.606 |
| CYP3A4 | scaffold | fixed_missing | rate_0.0653333333333333 | mixture_phi_0.33 | 0.102 | 0.033 | 0.163 | 0.999 |
| CYP3A4 | scaffold | fixed_missing | rate_0.0653333333333333 | mixture_phi_1.00 | 0.009 | -0.052 | 0.073 | 0.606 |
| CYP3A4 | scaffold | fixed_missing | rate_0.2 | broad | 0.199 | 0.156 | 0.240 | 1.000 |
| CYP3A4 | scaffold | fixed_missing | rate_0.2 | strong | -0.019 | -0.076 | 0.038 | 0.245 |
| CYP3A4 | scaffold | fixed_missing | rate_0.2 | mixture_phi_0.33 | 0.079 | 0.026 | 0.135 | 0.997 |
| CYP3A4 | scaffold | fixed_missing | rate_0.2 | mixture_phi_1.00 | -0.019 | -0.076 | 0.038 | 0.245 |
