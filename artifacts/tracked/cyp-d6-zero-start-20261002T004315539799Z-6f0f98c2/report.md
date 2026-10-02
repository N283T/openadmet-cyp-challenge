# Independent D6 restart: pilot result and diagnosis

Source model run: `classification/20261002T004016766672Z-e2d420d4`.

Every fitted model starts from source molecular structure or paired assay measurements. No historical ensemble, scores, call count, or prior selected weights are inputs.

| Family | Random AP | Random nested MCC | Scaffold AP | Scaffold nested MCC |
| --- | ---: | ---: | ---: | ---: |
| morgan_logistic | 0.3136 | 0.1159 | 0.2869 | 0.0861 |
| morgan_rdkit_classifier | 0.2612 | 0.0672 | 0.2707 | 0.0636 |
| paired_measurement_margin | 0.2521 | 0.0281 | 0.2647 | 0.0495 |

MCC uses outer calls with inner-selected thresholds. It is not the optimistic maximum MCC after fitting a threshold to all pooled OOF labels.

## Measurement uncertainty

Among the 1,364 supplied labels with direct pIC50 above 4:

- shift_boundary_overlap: 903 molecules.
- shift_negative_interval: 378 molecules.
- shift_positive_interval: 83 molecules.

Shift intervals are formed as TDI lower minus direct upper, and TDI upper minus direct lower. They are conservative combinations of marginal intervals, not joint posterior shift probabilities. Do not discard or relabel official labels from this diagnostic. Inferred positives need separate interpretation because their direct-arm potency is below the assay's reliable range.

## Operating-point transfer

| Family | Validation split | Selected threshold | Blind positive calls |
| --- | --- | ---: | ---: |
| morgan_logistic | random | 0.1200 | 300 |
| morgan_rdkit_classifier | random | 0.0800 | 592 |
| paired_measurement_margin | random | 0.1800 | 61 |
| morgan_logistic | scaffold | 0.1700 | 238 |
| morgan_rdkit_classifier | scaffold | 0.2400 | 192 |
| paired_measurement_margin | scaffold | -0.0400 | 382 |

These are unsubmitted diagnostic predictions. Counts emerge from the model and threshold; none was targeted. Their large split sensitivity is evidence against treating one OOF-optimal point threshold as established.

## Decision

Continue the independent restart. Logistic molecular features provide a weak ranking signal; the point-estimate paired-measurement recipe is a negative result. Neither justifies adopting a new submission. This does not authorize restoring the historical model or returning to a fixed call budget.

Next hypotheses must address paired-shift uncertainty, separation of conventional and inferred positives, and acquisition-population transfer. More trees or another small mixture with historical scores is not the next experiment.
