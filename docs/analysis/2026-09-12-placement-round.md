# Broad prediction calibration and placement comparison

The user authorized broad calibration experiments after native model reproduction
failed to materially improve D6. This round keeps predictions fixed and separates
train-label calibration from author-specified blind placement.

## Completed local calibration screen

Affine run `4a9dd6ab` tests 17 fixed methods: point/interval offsets, point/interval
affine maps, partial-strength interval maps, fixed slopes0.5/0.75/1.25/1.5/2,
and neighbor-only or neighbor-weighted interval fits. Nonlinear run `d94acc48`
tests14 methods including isotonic, quantile maps and regularized monotonic
piecewise interval maps. Identity appears in both:30 distinct methods total.
Both native mean and previous global predictions are evaluated on all four CYPs.
Consolidation run `a348113f` retains every method, without constructing an
endpoint-wise winner.

Each calibration fit excludes the assessed fold's labels. The base predictions
are saved development OOF, however, so this is not fully nested end-to-end model
validation. Test maps are fitted on all available OOF labels. Screen winners are
exploratory; fixed variants remain available separately.

| Fixed method | Native macro ST-RAE | Previous-global macro ST-RAE |
|---|---:|---:|
| Identity | 0.622574 | 0.582816 |
| Interval offset | 0.612447 | 0.580793 |
| Interval affine | 0.602363 | 0.581874 |
| Monotonic interval,5knots,no penalty | 0.601682 | 0.581381 |

For D6, native identity0.823983 moves to0.822792 with an interval offset and
0.819776 with the9-knot interval map (penalty0.05). Previous-global identity
0.825436 moves to0.820871 with the unpenalized9-knot interval map. These are
small changes. Isotonic squared-error calibration worsens D6 to0.865352/native
and0.876555/global. Quantile widening and aggressively increasing the slope also
worsen local validation; widening is not a generally valid correction for the
observed training population.

Independent audit `d1889f83` checks all canonical identities/labels/bounds,
reconstructs every affine prediction and reproduces all816 metric entries with
the published official scorer. Paired500-replicate bootstrap over1,185 D6 groups
finds native offset delta−0.001191 with95% interval[−0.004195,+0.001868];
9-knot interval delta−0.004207 with[−0.009076,+0.000995]. Both include zero.
These intervals condition on the saved predictions/maps and do not account for
method screening or refitting uncertainty. Nonlinear map reconstruction and
monotonicity independently pass audit `e43fa269`.

## Author placement and same-submission controls

The frozen author's actual `place()` function is executed for native scenarios
`46754704` and previously submitted-vector scenarios `d6e4f5d9`. Source ST-RAE
placement normalizes each750-compound endpoint vector by its sample mean/SD,
sets prescribed moments and clips below1. For D6 the specified mean/SD are
3.57/0.90. These are author's leaderboard-informed constants, not estimates
learned from our training data or confirmed moments of our scored subset.
Sensitivity variants retain this distinction; no synthetic test scores are
reported. See the detailed source report for the different R2-optimal formula
and the source solver's all750-versus-scored-subset limitation.

Crucially, the previously submitted input is verified numerically identical to
the actually sent file (all3,000 values). D6-only placement preserves the other
2,250 values exactly and preserves D6 ranks. Its mean/SD change from4.638156/
0.405434 to3.570000/0.900000; minimum/maximum become2.16949/6.55948, with no
lower clipping required. This creates a controlled next leaderboard comparison
that changes calibration alone, rather than changing model and calibration
simultaneously.

Packaging `20260912T113446532334Z-31ac9555` freezes three750-row CSV files,
all passing official identity/schema validation:

- `raw.csv`: unchanged previously submitted prediction control.
- `d6_only_author_strae.csv`: only D6 uses source ST-RAE placement.
- `author_solved_strae.csv`: source placement on all four endpoints.

No upload occurred. The large-blind-shift hypothesis has a concrete test artifact,
but its leaderboard benefit has not been measured. Local interval calibration
helps the newly trained native model overall, while its D6 gain remains small.
The same-vector D6-only CSV is the most direct next submission for isolating the
placement hypothesis; do not claim it already improves the leaderboard.

## Supporting reports

- [Nonlinear maps](2026-09-12-scp-nonlinear-calibration.md)
- [Source placement and scenarios](2026-09-12-scp-placement-scenarios.md)
- [Independent reconstruction and uncertainty](2026-09-12-calibration-independent-audit.md)

Full-OOF nonlinear test inference `5d73ca01` additionally saves identity,
5-knot/no-penalty and9-knot/0.05 variants,3,000 keyed predictions each, plus
8 full-data fitted maps. Failed first inference `76791782` is retained; it
encountered float32/float64 concatenation after fitting, corrected in the fresh
run by explicit dtype normalization. No completed output was overwritten.
