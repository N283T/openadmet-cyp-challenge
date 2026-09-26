# Cross-CYP weak-label transfer under synthetic selection masking

## Question

The D6 direct-inhibition training labels are selected almost entirely from
strong primary-screen responses. Can weak-primary direct-inhibition labels from
other CYPs reveal a missing low-potency regime? The official assay description
states that the library was screened first and a selected set proceeded to
dose-response measurement. The prior [selection audit](2026-09-11-label-selection.md)
found no public D6 direct labels among 2,823 compounds with D6 primary log2FC
at least −0.5. Other CYPs have some direct labels in that range, allowing a
synthetic masking test before changing D6.

## Design

Recorded run `regression/20260926T035054510431Z-7a02201c` holds out the
same parent-series and scaffold folds used in prior work. For each target in
1A2, 2C9 and 3A4, it retains only target direct labels whose measured primary
log2FC is below −0.5 for training. This hides the weak-primary labels and
also 523 3A4 held-out labels with missing public primary values. The
`selected_only` arm fits the retained target labels;
`selected_plus_other_weak` adds weak-primary direct
labels from the other two CYPs at weight 0.25. Donor labels for every held-out
compound are excluded in that fold. An `all_labels_oracle` arm uses target
weak labels and is an unreachable reference for D6. All arms use the same
647 frozen score/embedding features and LightGBM L1 recipe. They are matched
within this experiment, not directly comparable to the submitted D6 ensemble.

The held-out weak-screen cohorts contain 249 1A2, 410 2C9 and 522 3A4
measurements. The generalization check is the held-out cohort, not training
fit. ST-RAE is calculated on all held-out direct labels. Weak-screen values
are mean soft error, because computing the official metric's denominator
separately by a truth-defined stratum would not describe the full challenge.

| Split | Target | All-label ST-RAE, selected only → donor weak | Weak-screen mean soft error, selected only → donor weak | Selected-like mean soft error, selected only → donor weak |
| --- | --- | ---: | ---: | ---: |
| Parent | 1A2 | .829 → .735 | 1.049 → .584 | .242 → .292 |
| Parent | 2C9 | .678 → .591 | .317 → .137 | .138 → .188 |
| Parent | 3A4 | .543 → .493 | .369 → .268 | .244 → .248 |
| Scaffold | 1A2 | .826 → .731 | 1.057 → .563 | .239 → .294 |
| Scaffold | 2C9 | .680 → .585 | .316 → .135 | .140 → .187 |
| Scaffold | 3A4 | .538 → .483 | .377 → .259 | .245 → .248 |

The donor labels address a real selection gap in these **synthetic** settings:
the weak-screen group improves on all targets and both splits. They also harm
the already labelled strong-response group, especially for 1A2 and 2C9.
The all-label oracle remains better on every target. These results do not
show that the D6 test set has the same weak-screen prevalence or response
distribution as the synthetic target cohorts.

## Direct D6 transfer check

Recorded run `regression/20260926T035426767911Z-9b7ef6c0` fits D6 direct
labels with and without all three other-CYP weak-primary labels, excluding
held-out compounds' donor measurements from each fold. The donor weight grid
was 0, 0.1, 0.25 and 0.5; the underlying D6 model and features are identical
within the comparison. Each split has 1,457 nonanchor held-out D6 labels.

| Donor weight | Parent ST-RAE | Scaffold ST-RAE | Parent mean prediction |
| ---: | ---: | ---: | ---: |
| 0 | .836 | .827 | 4.818 |
| 0.1 | .868 | .859 | 4.750 |
| 0.25 | .953 | .939 | 4.681 |
| 0.5 | 1.022 | 1.022 | 4.617 |

Naive pooling systematically lowers D6 predictions and worsens known D6
held-out errors. This is consistent with the strong-label harm in the masking
benchmark, but does not measure error on the unlabelled D6 weak-screen region.
Do not select a candidate or submit this pooling model based on the synthetic
benefit alone.

## Endpoint-conditioned check

Recorded run `regression/20260926T035618807049Z-46298171` repeats the D6
check with three additional donor-endpoint indicators. D6 examples and D6
queries have all-zero indicators; each donor example marks its enzyme. The
zero-weight baseline is numerically identical to the direct D6 control.

| Donor weight | Parent ST-RAE | Scaffold ST-RAE |
| ---: | ---: | ---: |
| 0 | .836 | .827 |
| 0.1 | .841 | .838 |
| 0.25 | .847 | .831 |
| 0.5 | .848 | .837 |

Separating the donor identities prevents most of the damage from naive pooling,
but no tested weight improves either split over the matched D6-only baseline.
This small screen does not rule out a more specialized shared representation
or loss. It supplies no validated D6 replacement.

## Decision and limits

The transferable information in other-CYP weak direct labels is real in
the synthetic target-label masking test. A single pooled direct regression
model is not a useful D6 replacement because it damages the only measurable
D6 validation population. Endpoint indicators reduce that damage but do not
produce a gain. Any follow-up must preserve D6-specific calibration while
learning the shared weak-activity signal. The
critical missing observation is still D6 direct potency on compounds whose
D6 primary response is weak; the hidden challenge labels cannot establish
that error directly.

Compact metrics and inventories are promoted under
[`artifacts/tracked/regression-selection-mask-cross-cyp-20260926/`](../../artifacts/tracked/regression-selection-mask-cross-cyp-20260926/README.md).
Keyed OOF predictions remain in the recorded runs.
