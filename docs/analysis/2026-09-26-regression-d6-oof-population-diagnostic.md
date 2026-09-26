# CYP2D6 regression: existing OOF population diagnostic

Recorded run `regression/20260926T013148480112Z-72cf4682` aligns the current
two-member ensemble's unplaced parent-series OOF with its keyed labels and the saved
other-CYP pseudo-analog set. The 1,457 OOF predictions cover nonanchor D6
labels. This is a descriptive check of an existing model, not a new split or
candidate selection.

| Population | n | Mean observed pIC50 | Mean predicted | ST-RAE | MAE | Spearman |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| All nonanchors | 1,457 | 4.772 | 4.801 | .8284 | .555 | .475 |
| Other-CYP pseudo-analogs | 163 | 4.938 | 4.956 | .8129 | .600 | .534 |
| Outside that proxy | 1,294 | 4.751 | 4.781 | .8373 | .549 | .459 |

The 163 labelled pseudo-analogs are **more** potent than the full OOF population.
Their slightly better ST-RAE and rank correlation therefore do not validate
transfer to the weaker-activity D6 test distribution described by the organizer.
The proxy and nonproxy ST-RAE denominators use their own population means, so
their ratio values should not be interpreted as a paired improvement.
One of the 164 labelled proxy compounds in the earlier coverage audit is an
anchor and has no OOF prediction.

Outcome-conditioned stress strata expose asymmetric point errors:

| Observed pIC50 | n | Mean prediction minus observation | MAE | Mean soft error |
| --- | ---: | ---: | ---: | ---: |
| Below 4 | 128 | +1.836 | 1.836 | .690 |
| 4 to below 5 | 801 | +.151 | .295 | .163 |
| At least 5 | 528 | −.596 | .638 | .531 |

The model's predictions occupy a narrower range than the observed activity,
and its low-activity tail is strongly overpredicted when stratifying by truth.
This pattern was already present in the September 11 baseline diagnosis. It
does not by itself justify widening or shifting predictions: weak predictors
can correctly shrink toward the conditional mean, and truth-defined strata
cannot be identified at test time. Measurement bounds reduce the soft error, so raw
point MAE does not equal the challenge loss. The below-4 group contributes
88.37 of the 499.40 total soft-error numerator despite being 8.8% of OOF
labels. Only 16 pseudo-analog OOF labels lie below 4. These strata are defined
from the **outcome**, which is unavailable for hidden test molecules; they
cannot support an activity-conditioned deployment rule. Their within-stratum
ST-RAE is deliberately omitted because the constant-mean denominator can
collapse after conditioning on the response. The initial run `ffb007f4`
computed those unstable ratios and is marked superseded; `72cf4682` contains
the corrected output.

## Decision

Keep the submitted D6 vector and parent/scaffold readouts as controls. The
low-activity stress result identifies a useful next question: whether a
fold-safe feature or training objective improves low-activity error **and**
overall ST-RAE and ordering on both prespecified splits. The available proxy
cannot serve as a model-selection gate. A D6-specific acquisition view still
needs enough independent weak-activity labels, as specified in the
[validation design note](2026-09-26-regression-d6-oof-design.md).
The first [weak-activity intervention tests](2026-09-26-regression-d6-weak-activity-tests.md)
did not improve ST-RAE across the two splits.
Compact aggregate results are in
[`artifacts/tracked/regression-d6-oof-population-20260926/`](../../artifacts/tracked/regression-d6-oof-population-20260926/README.md).
