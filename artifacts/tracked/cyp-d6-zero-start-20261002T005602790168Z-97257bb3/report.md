# Independent CYP2D6 restart: consolidated development evidence

All models below were trained from source molecules/assays with no historical ensemble, fitted encoder, prediction or call budget. Four outer/three inner folds exclude each held-out molecule's complete assay record from fitting. Outer calls use inner-selected thresholds, evaluated against all supplied labels.

| Family | Random AP | Random nested MCC | Scaffold AP | Scaffold nested MCC |
| --- | ---: | ---: | ---: | ---: |
| morgan_logistic | 0.3136 | 0.1159 | 0.2869 | 0.0861 |
| morgan_rdkit_classifier | 0.2612 | 0.0672 | 0.2707 | 0.0636 |
| paired_measurement_margin | 0.2521 | 0.0281 | 0.2647 | 0.0495 |
| confidence_weighted_logistic | 0.3115 | 0.0935 | 0.2941 | 0.1124 |
| regime_multiclass | 0.3119 | 0.1161 | 0.2804 | 0.0560 |
| confidence_weighted_regime | 0.3092 | 0.0962 | 0.2833 | 0.0911 |
| fresh_multiassay_network | 0.3382 | 0.1282 | 0.3243 | 0.1123 |
| fresh_d6_network | 0.3383 | 0.1514 | 0.3174 | 0.1276 |

## Paired development uncertainty

Both sides use the same labelled molecules and the same scaffold-group resamples. Intervals are conditional on saved predictions; they do not include retraining uncertainty, model-selection effects, or transfer to the unknown Live scored population.

| Split | Comparison | Metric | Mean delta | 95% interval |
| --- | --- | --- | ---: | --- |
| random | multiassay minus morgan_logistic | AP | +0.0243 | [+0.0022, +0.0457] |
| random | multiassay minus morgan_logistic | nested_calls_MCC | +0.0127 | [-0.0435, +0.0647] |
| random | multiassay minus fresh_d6_network | AP | -0.0010 | [-0.0170, +0.0165] |
| random | multiassay minus fresh_d6_network | nested_calls_MCC | -0.0221 | [-0.0637, +0.0207] |
| scaffold | multiassay minus morgan_logistic | AP | +0.0370 | [+0.0119, +0.0623] |
| scaffold | multiassay minus morgan_logistic | nested_calls_MCC | +0.0250 | [-0.0261, +0.0769] |
| scaffold | multiassay minus fresh_d6_network | AP | +0.0063 | [-0.0118, +0.0238] |
| scaffold | multiassay minus fresh_d6_network | nested_calls_MCC | -0.0170 | [-0.0647, +0.0312] |

## Unsubmitted diagnostic operating points

| Family | Validation | Threshold | Blind positive calls |
| --- | --- | ---: | ---: |
| morgan_logistic | random | 0.1200 | 300 |
| morgan_rdkit_classifier | random | 0.0800 | 592 |
| paired_measurement_margin | random | 0.1800 | 61 |
| morgan_logistic | scaffold | 0.1700 | 238 |
| morgan_rdkit_classifier | scaffold | 0.2400 | 192 |
| paired_measurement_margin | scaffold | -0.0400 | 382 |
| confidence_weighted_logistic | random | 0.1000 | 339 |
| regime_multiclass | random | 0.0600 | 347 |
| confidence_weighted_regime | random | 0.0700 | 346 |
| confidence_weighted_logistic | scaffold | 0.1600 | 257 |
| regime_multiclass | scaffold | 0.0900 | 303 |
| confidence_weighted_regime | scaffold | 0.0900 | 304 |
| fresh_multiassay_network | random | 0.1400 | 153 |
| fresh_multiassay_network | scaffold | 0.9000 | 58 |
| fresh_d6_network | random | 0.0400 | 203 |
| fresh_d6_network | scaffold | 0.9200 | 94 |

No count was selected in advance. These point thresholds are fitted on pooled development OOF scores and are not separate unbiased tests. No submission payload was prepared and no upload was made.

## Provenance

- v1: `classification/20261002T004016766672Z-e2d420d4`.
- v2: `classification/20261002T004536560534Z-21c1509c`.
- v3: `classification/20261002T005041934540Z-c62b6e78`.
- v4: `classification/20261002T005413496423Z-f92f58e6`.
