# Independent CYP2D6 restart: consolidated development evidence

All models below were trained from source molecules/assays with no historical ensemble, fitted encoder, prediction or call budget. Four outer/three inner folds exclude each held-out molecule's complete assay record from fitting. Outer calls use inner-selected thresholds, evaluated against all supplied labels.

| Family | Validation | AP | Nested MCC | AUROC |
| --- | --- | ---: | ---: | ---: |
| morgan_logistic | random | 0.3136 | 0.1159 | 0.6257 |
| morgan_rdkit_classifier | random | 0.2612 | 0.0672 | 0.5690 |
| paired_measurement_margin | random | 0.2521 | 0.0281 | 0.5342 |
| morgan_logistic | scaffold | 0.2869 | 0.0861 | 0.5924 |
| morgan_rdkit_classifier | scaffold | 0.2707 | 0.0636 | 0.5777 |
| paired_measurement_margin | scaffold | 0.2647 | 0.0495 | 0.5581 |
| confidence_weighted_logistic | random | 0.3115 | 0.0935 | 0.6260 |
| regime_multiclass | random | 0.3119 | 0.1161 | 0.6247 |
| confidence_weighted_regime | random | 0.3092 | 0.0962 | 0.6253 |
| confidence_weighted_logistic | scaffold | 0.2941 | 0.1124 | 0.5964 |
| regime_multiclass | scaffold | 0.2804 | 0.0560 | 0.5852 |
| confidence_weighted_regime | scaffold | 0.2833 | 0.0911 | 0.5873 |
| fresh_multiassay_network | random | 0.3382 | 0.1282 | 0.6323 |
| fresh_multiassay_network | scaffold | 0.3243 | 0.1123 | 0.6093 |
| fresh_d6_network | random | 0.3383 | 0.1514 | 0.6386 |
| fresh_d6_network | scaffold | 0.3174 | 0.1276 | 0.6087 |
| morgan_logistic | similarity | 0.2907 | 0.0482 | 0.6119 |
| fresh_d6_network | similarity | 0.3038 | 0.1184 | 0.6134 |
| fresh_multiassay_network | similarity | 0.3016 | 0.0881 | 0.6124 |

## Paired development uncertainty

Both sides use the same labelled molecules and the same scaffold-group resamples. Intervals are conditional on saved predictions; they do not include retraining uncertainty, model-selection effects, or transfer to the unknown Live scored population.

| Split | Comparison | Metric | Mean delta | 95% interval |
| --- | --- | --- | ---: | --- |
| random | multiassay minus morgan_logistic | AP | +0.0243 | [+0.0022, +0.0457] |
| random | multiassay minus morgan_logistic | nested_calls_MCC | +0.0127 | [-0.0435, +0.0647] |
| random | multiassay minus fresh_d6_network | AP | -0.0010 | [-0.0170, +0.0165] |
| random | multiassay minus fresh_d6_network | nested_calls_MCC | -0.0221 | [-0.0637, +0.0207] |
| random | D6 network minus morgan_logistic | AP | +0.0253 | [+0.0040, +0.0468] |
| random | D6 network minus morgan_logistic | nested_calls_MCC | +0.0348 | [-0.0153, +0.0835] |
| scaffold | multiassay minus morgan_logistic | AP | +0.0370 | [+0.0119, +0.0623] |
| scaffold | multiassay minus morgan_logistic | nested_calls_MCC | +0.0250 | [-0.0261, +0.0769] |
| scaffold | multiassay minus fresh_d6_network | AP | +0.0063 | [-0.0118, +0.0238] |
| scaffold | multiassay minus fresh_d6_network | nested_calls_MCC | -0.0170 | [-0.0647, +0.0312] |
| scaffold | D6 network minus morgan_logistic | AP | +0.0307 | [+0.0108, +0.0534] |
| scaffold | D6 network minus morgan_logistic | nested_calls_MCC | +0.0421 | [-0.0183, +0.1040] |
| similarity | multiassay minus morgan_logistic | AP | +0.0110 | [-0.0081, +0.0316] |
| similarity | multiassay minus morgan_logistic | nested_calls_MCC | +0.0405 | [-0.0217, +0.1006] |
| similarity | multiassay minus fresh_d6_network | AP | -0.0022 | [-0.0176, +0.0132] |
| similarity | multiassay minus fresh_d6_network | nested_calls_MCC | -0.0285 | [-0.0695, +0.0113] |
| similarity | D6 network minus morgan_logistic | AP | +0.0132 | [-0.0066, +0.0351] |
| similarity | D6 network minus morgan_logistic | nested_calls_MCC | +0.0690 | [+0.0174, +0.1239] |

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
| morgan_logistic | similarity | 0.2100 | 209 |
| fresh_d6_network | similarity | 0.0100 | 227 |
| fresh_multiassay_network | similarity | 0.0200 | 219 |

No count was selected in advance. These point thresholds are fitted on pooled development OOF scores and are not separate unbiased tests. No submission payload was prepared and no upload was made.

## Provenance

- v1: `classification/20261002T004016766672Z-e2d420d4`.
- v2: `classification/20261002T004536560534Z-21c1509c`.
- v3: `classification/20261002T005041934540Z-c62b6e78`.
- v4: `classification/20261002T005413496423Z-f92f58e6`.
- v5: `classification/20261002T005927124518Z-ce5df58c`.
