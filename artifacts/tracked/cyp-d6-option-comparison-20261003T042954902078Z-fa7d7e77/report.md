# Three D6 options in public-metric units

Deployed point MCC 0.36165 against a reported Live 0.364248.

## Expected change and spread

| option | variant | calls | changed_calls | scored_changed | scenario | added_precision | expected_delta | sd_delta | interval_low | interval_high | p_delta_positive |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gate_filter | l4_filter_tau0.2 | 49 | 7 | 4.0000 | base_rate | +0.0706 | -0.0515 | +0.0441 | -0.1198 | +0.0000 | +0.0442 |
| gate_filter | l4_filter_tau0.2 | 49 | 7 | 4.0000 | offline_mixture_ratio | +0.2652 | -0.0204 | +0.0526 | -0.0799 | +0.0399 | +0.2186 |
| gate_filter | l4_filter_tau0.2 | 49 | 7 | 4.0000 | unchanged | +0.3929 | +0.0000 | +0.0552 | -0.0799 | +0.0799 | +0.3590 |
| gate_filter | l4_filter_tau0.2 | 49 | 7 | 4.0000 | neighbour_rate_ratio | +0.2960 | -0.0154 | +0.0534 | -0.0799 | +0.0399 | +0.2518 |
| gate_filter | l4_filter_tau0.3 | 49 | 21 | 12.0000 | base_rate | +0.0706 | -0.1545 | +0.0763 | -0.2397 | -0.0399 | +0.0102 |
| gate_filter | l4_filter_tau0.3 | 49 | 21 | 12.0000 | offline_mixture_ratio | +0.2652 | -0.0612 | +0.0913 | -0.1598 | +0.0399 | +0.1861 |
| gate_filter | l4_filter_tau0.3 | 49 | 21 | 12.0000 | unchanged | +0.3929 | +0.0001 | +0.0956 | -0.1198 | +0.1198 | +0.4177 |
| gate_filter | l4_filter_tau0.3 | 49 | 21 | 12.0000 | neighbour_rate_ratio | +0.2475 | -0.0699 | +0.0901 | -0.1997 | +0.0399 | +0.1584 |
| gate_filter | l4_filter_tau0.4 | 49 | 29 | 16.6000 | base_rate | +0.0706 | -0.2185 | +0.0899 | -0.3196 | -0.1198 | +0.0039 |
| gate_filter | l4_filter_tau0.4 | 49 | 29 | 16.6000 | offline_mixture_ratio | +0.2652 | -0.0861 | +0.1075 | -0.2397 | +0.0399 | +0.1621 |
| gate_filter | l4_filter_tau0.4 | 49 | 29 | 16.6000 | unchanged | +0.3929 | +0.0007 | +0.1130 | -0.1598 | +0.1598 | +0.4317 |
| gate_filter | l4_filter_tau0.4 | 49 | 29 | 16.6000 | neighbour_rate_ratio | +0.2410 | -0.1025 | +0.1060 | -0.2397 | +0.0399 | +0.1251 |
| gate_filter | l4_filter_tau0.5 | 49 | 40 | 22.9000 | base_rate | +0.0706 | -0.2858 | +0.0917 | -0.3995 | -0.1598 | +0.0013 |
| gate_filter | l4_filter_tau0.5 | 49 | 40 | 22.9000 | offline_mixture_ratio | +0.2652 | -0.1070 | +0.1147 | -0.2397 | +0.0399 | +0.1383 |
| gate_filter | l4_filter_tau0.5 | 49 | 40 | 22.9000 | unchanged | +0.3929 | +0.0111 | +0.1215 | -0.1598 | +0.1598 | +0.4570 |
| gate_filter | l4_filter_tau0.5 | 49 | 40 | 22.9000 | neighbour_rate_ratio | +0.1790 | -0.1859 | +0.1070 | -0.3196 | -0.0399 | +0.0352 |
| gate_filter | l4_x_p_active | 49 | 30 | 17.1000 | base_rate | +0.0706 | -0.2184 | +0.0897 | -0.3196 | -0.1198 | +0.0039 |
| gate_filter | l4_x_p_active | 49 | 30 | 17.1000 | offline_mixture_ratio | +0.2652 | -0.0865 | +0.1075 | -0.2397 | +0.0399 | +0.1614 |
| gate_filter | l4_x_p_active | 49 | 30 | 17.1000 | unchanged | +0.3929 | +0.0004 | +0.1130 | -0.1598 | +0.1598 | +0.4305 |
| gate_filter | l4_x_p_active | 49 | 30 | 17.1000 | neighbour_rate_ratio | +0.2003 | -0.1300 | +0.1030 | -0.2796 | +0.0000 | +0.0736 |
| gate_filter | l4_x_sqrt_p_active | 49 | 24 | 13.7000 | base_rate | +0.0706 | -0.1805 | +0.0825 | -0.2796 | -0.0799 | +0.0070 |
| gate_filter | l4_x_sqrt_p_active | 49 | 24 | 13.7000 | offline_mixture_ratio | +0.2652 | -0.0713 | +0.0987 | -0.1997 | +0.0399 | +0.1771 |
| gate_filter | l4_x_sqrt_p_active | 49 | 24 | 13.7000 | unchanged | +0.3929 | +0.0003 | +0.1032 | -0.1198 | +0.1198 | +0.4247 |
| gate_filter | l4_x_sqrt_p_active | 49 | 24 | 13.7000 | neighbour_rate_ratio | +0.2390 | -0.0857 | +0.0968 | -0.1997 | +0.0399 | +0.1376 |
| extension_94 | l4_top94 | 94 | 45 | 22.5000 | uncalibrated_mixture | +0.1330 | -0.0389 | +0.0492 | -0.0983 | +0.0252 | +0.1587 |
| extension_94 | l4_top94 | 94 | 45 | 22.5000 | break_even | +0.1880 | -0.0015 | +0.0567 | -0.0674 | +0.0870 | +0.3989 |
| extension_94 | l4_top94 | 94 | 45 | 22.5000 | calibrated_by_1_99 | +0.2650 | +0.0509 | +0.0639 | -0.0366 | +0.1488 | +0.7316 |
| hold | deployed_oct03 | 49 | 0 | +0.0000 | unchanged | – | +0.0000 | +0.0000 | +0.0000 | +0.0000 | +0.0000 |

## Sensitivity to the scored fraction of the changed ranks

| option | variant | scored_fraction | expected_delta | interval_low | interval_high |
|---|---|---|---|---|---|
| gate_filter | l4_filter_tau0.2 | +0.4990 | -0.0152 | -0.0799 | +0.0399 |
| extension_94 | l4_top94 | +0.4990 | -0.0388 | -0.0983 | +0.0252 |
| gate_filter | l4_filter_tau0.2 | +0.5714 | -0.0203 | -0.0799 | +0.0399 |
| extension_94 | l4_top94 | +0.5714 | -0.0413 | -0.1150 | +0.0345 |
| gate_filter | l4_filter_tau0.2 | +0.6500 | -0.0255 | -0.0799 | +0.0399 |
| extension_94 | l4_top94 | +0.6500 | -0.0430 | -0.0973 | +0.0197 |

## Probability that the pre-registered extension gives a decisive outcome

| scenario | added_precision | p_above_plus_080 | p_below_minus_080 | p_decisive | p_inside |
|---|---|---|---|---|---|
| uncalibrated_mixture | +0.1330 | +0.0203 | +0.1897 | +0.2100 | +0.7900 |
| break_even | +0.1880 | +0.1021 | +0.0619 | +0.1641 | +0.8359 |
| calibrated_by_1_99 | +0.2650 | +0.3594 | +0.0103 | +0.3697 | +0.6303 |

