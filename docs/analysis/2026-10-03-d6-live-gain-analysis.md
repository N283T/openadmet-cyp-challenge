# Sizing the 2026-10-03 D6 Live gain and attributing it

Recorded run `20261003T012000305555Z-bca66716` (classification, 11 s).
Code: `scripts/evaluate_cyp_d6_live_gain.py`,
`scripts/launch_cyp_d6_live_gain.py`, `scripts/retain_cyp_d6_live_gain.py`,
`configs/classification/d6-live-gain-v1.json`.
Compact outputs with SHA-256 provenance:
`artifacts/tracked/cyp-d6-live-gain-20261003T012000305555Z-bca66716/`.
Nothing was fitted, no candidate or prepared submission was created, no upload
occurred and the challenge API was not contacted.

## Question

The L4 rank-mean payload `20261002T223229553466Z-a08301d2` moved the D6 Live MCC
from the project's previous best .308635 to .364248 at an unchanged 49 calls
([RL0217](../research/2026-10.md)). The reported bootstrap standard deviation is
.087528, so the move is 0.64 reported standard deviations. This analysis asks
four narrow statistical questions and nothing else:

1. How many individual scored row outcomes does the move actually consist of?
2. How surprising is a move that size under explicit nulls?
3. Which of the two score columns in the rank mean produced the new calls?
4. Did the offline estimate for this exact variant predict the realised value,
   and what does the new transition do to the frozen concordance tally?

No biological or chemical claim is made anywhere below. All quantities are
reconstructions from published leaderboard aggregates plus frozen score columns.

## Part 1 - what the published aggregates imply

The leaderboard scores an unknown subset of the 750 blind rows and publishes the
bootstrap mean and standard deviation of Accuracy, Precision, Recall, F1 and
MCC. Three independent routes to the scored-subset size agree.

**Variance inversion.** If the published standard deviations come from
resampling the scored rows, then `P ~ r(1-r)/sd_r^2` on recall estimates the
scored positive count, `C ~ p(1-p)/sd_p^2` on precision estimates the scored
call count, and `N ~ a(1-a)/sd_a^2` on accuracy estimates the scored subset.

**Prevalence identity.** With `r` recall and `p` precision,
`prevalence = (1-a)/(1-2r+r/p)` and `call fraction = prevalence*r/p`. Both are
scale free, so combining the call fraction with the precision inversion gives a
second estimate of `N`.

| State | Submitted | Calls | MCC | P from sd(recall) | C from sd(prec.) | N from sd(acc.) | prevalence | call fraction | N from identity |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| sep13 | 2026-09-13 11:38 UTC | 49 | .2710 | 26.31 | 29.44 | 375.9 | .07061 | .07947 | 370.2 |
| sep14 | 2026-09-14 02:08 UTC | 49 | .2786 | 26.31 | 28.36 | 369.8 | .07064 | .07665 | 369.9 |
| sep15 | 2026-09-15 13:00 UTC | 49 | .3086 | 25.88 | 29.37 | 369.5 | .07061 | .07936 | 370.0 |
| sep18 | 2026-09-18 09:45 UTC | 49 | .1097 | 24.99 | 28.70 | 349.6 | .07064 | .08171 | 351.3 |
| sep23 | 2026-09-23 18:03 UTC | 49 | .2348 | 24.56 | 26.74 | 366.7 | .07058 | .07693 | 347.5 |
| oct01 | 2026-10-01 14:13 UTC | 49 | .1050 | 23.25 | 29.83 | 393.6 | .07034 | .08519 | 350.1 |
| oct02 | 2026-10-02 12:08 UTC | 175 | .1177 | 26.04 | 90.03 | 386.1 | .07062 | .23870 | 377.2 |
| oct03 | 2026-10-03 00:29 UTC | 49 | .3642 | 26.33 | 27.41 | 402.6 | .07060 | .07396 | 370.6 |

The identity returns the same prevalence, .0706, from every one of the eight
states to four decimal places, including the 175-call payload whose call
fraction is three times larger. That is the strongest available evidence that
one fixed scored mask and one fixed label set underlie all eight rows, and it
reproduces the 7.06% figure already used in
[the evaluation protocol](../classification-evaluation.md). The accuracy
inversion is the noisiest of the three routes because accuracy is dominated by
true negatives.

### Joint integer fit

Because all eight states share one scored subset and one label set, the scored
size `N` and the scored positive count `P` are common parameters and only the
true-positive count and the scored call count vary per state. A grid search
over `N` in 300-520 and `P` in 16-44 picks, per state, the integer
(true positives, scored calls) minimising the sum over the five metrics of the
squared deviation from the published mean in units of the published standard
deviation, and sums that loss over the eight states.

The minimum is **N = 374 scored rows of 750, P = 26 scored positives**
(prevalence .0695), total loss .0639. Every one of the ten best (N, P) pairs has
P = 26 and N between 365 and 377, and every one assigns the same
true-positive counts to all eight states. The scored subset is therefore close
to half the blind set, consistent with the A4 figure of about 368 inferred
independently in [the A4 diagnosis](2026-10-03-a4-diagnosis.md).

Reported values against a 40,000-resample bootstrap of the fitted confusion
matrix:

| State | Calls | Scored calls | TP | FP | FN | MCC rep. | MCC sim. | Prec. rep. | Prec. sim. | Recall rep. | Recall sim. |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| sep13 | 49 | 30 | 9 | 21 | 17 | .2710 | .2668 | .3080 | .3003 | .3469 | .3472 |
| sep14 | 49 | 29 | 9 | 20 | 17 | .2786 | .2722 | .3197 | .3096 | .3469 | .3454 |
| sep15 | 49 | 30 | 10 | 20 | 16 | .3086 | .3043 | .3412 | .3329 | .3835 | .3843 |
| sep18 | 49 | 31 | 5 | 26 | 21 | .1097 | .1078 | .1655 | .1613 | .1914 | .1924 |
| sep23 | 49 | 29 | 8 | 21 | 18 | .2348 | .2343 | .2800 | .2763 | .3052 | .3080 |
| oct01 | 49 | 32 | 5 | 27 | 21 | .1050 | .1040 | .1588 | .1565 | .1924 | .1927 |
| oct02 | 175 | 89 | 11 | 78 | 15 | .1177 | .1184 | .1247 | .1236 | .4215 | .4234 |
| oct03 | 49 | 28 | 11 | 17 | 15 | .3642 | .3603 | .4024 | .3931 | .4215 | .4235 |

Two independent consistency checks fall out. The sep13 and sep14 recalls are
identical to sixteen digits, which forces the same true-positive count and the
same positive denominator for those two states, and the fit assigns 9 to both.
The 175-call payload is fitted with 89 scored calls, a scored fraction of .509
against the overall .499, whereas the 49-call payloads are fitted at .57-.65
scored, so the top-ranked calls are modestly enriched in scored rows.

### The change in row outcomes

| Quantity | Sep 15 | Oct 3 | Change |
|---|---:|---:|---:|
| Scored rows | 374 | 374 | 0 |
| Scored positives | 26 | 26 | 0 |
| Total calls | 49 | 49 | 0 |
| Scored calls | 30 | 28 | -2 |
| True positives | 10 | 11 | **+1** |
| False positives | 20 | 17 | **-3** |
| False negatives | 16 | 15 | -1 |

Holding the scale-free identity fixed and varying only the scored-subset size
gives the sensitivity band:

| Scored rows | Sep 15 TP | Oct 3 TP | TP change | Sep 15 FP | Oct 3 FP | FP change |
|---:|---:|---:|---:|---:|---:|---:|
| 330 | 8.94 | 9.82 | +0.88 | 17.25 | 14.58 | -2.67 |
| 350 | 9.48 | 10.42 | +0.94 | 18.30 | 15.47 | -2.83 |
| 366 | 9.91 | 10.89 | +0.98 | 19.14 | 16.18 | -2.96 |
| 374 | 10.13 | 11.13 | +1.00 | 19.55 | 16.53 | -3.02 |
| 390 | 10.56 | 11.61 | +1.05 | 20.39 | 17.24 | -3.15 |
| 410 | 11.10 | 12.20 | +1.10 | 21.44 | 18.12 | -3.32 |
| 430 | 11.64 | 12.80 | +1.15 | 22.48 | 19.00 | -3.48 |

**Stated plainly: the headline .0556 MCC gain is one scored molecule changing
from a false negative to a true positive, together with three scored molecules
changing from false positive to uncalled.** Across every plausible scored-subset
size from 330 to 430 the change is +0.9 to +1.2 true positives and -2.7 to -3.5
false positives. The denominator is 26 scored positives.

At this operating point one true positive is worth .0399 MCC and one false
positive only .0085, so the metric is effectively quantised in units of about
.04 per scored positive. Decomposing the move: raising the true-positive count
from 10 to 11 at 30 scored calls is +.0387, and dropping from 30 to 28 scored
calls at 11 true positives adds +.0166. **About 70% of the headline gain is one
molecule.** The reported bootstrap standard deviation of .0875 corresponds to
2.2 true positives.

## Part 2 - how surprising is a change this size

Two nulls are reported. Both hold the fitted confusion matrix of the earlier
state fixed and resample only which molecules the swap touches. Point MCC
deltas are used; a check on 2,000 draws with an inner 400-resample bootstrap
shows the bootstrap-mean delta differs from the point delta by +.0031 on
average with the same spread (.0591 against .0588), so the two are
interchangeable at the precision that matters here.

### Null A - an unskilled swap, as specified

Remove `k` of the 49 calls uniformly at random and replace them with `k` rows
drawn uniformly from the 701 uncalled rows. The composition of both draws is
multivariate hypergeometric over (scored positive, scored negative, unscored)
at the fitted base rates. 200,000 draws per transition.

| Transition | Swapped | Observed | In reported sd | Null mean | Null sd | Null 5% | Null 50% | Null 95% | P(null >= observed) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| sep15 -> oct03 | 17 | +.0556 | +0.64 | -.1068 | .0574 | -.2012 | -.1045 | -.0141 | **.0018** |
| sep14 -> sep15 | 9 | +.0301 | +0.35 | -.0507 | .0443 | -.1230 | -.0458 | +.0177 | .0322 |
| sep15 -> sep18 | 27 | -.1989 | -2.70 | -.1734 | .0637 | -.2709 | -.1744 | -.0645 | .6489 |

This null is strongly negative-shifted, because replacing calls whose fitted
precision is .33 by rows drawn at a .07 base rate destroys performance. The
observed Sep 15 -> Oct 3 move sits at the 99.8th percentile, so an unskilled
17-call swap essentially never reproduces it. That answers the question as
posed, but the bar is low: it only establishes that the 17 added molecules are
better than 17 random molecules from the bottom 701 of the deployed score. The
Sep 15 -> Sep 18 loss, by contrast, is indistinguishable from an unskilled swap
(p = .65), which is exactly what a real failure looks like under this null.

### Null B - an exchangeable split of the changed calls

The sharper question is whether the recipe ordered the 34 changed molecules
correctly, not whether they beat random molecules. Under Null B the union of
the removed and added calls is fixed and split at random into two halves of
`k`; the retained calls contribute a fixed count. The composition of the
changed set is not observed, so every (retained true positives, retained scored
calls) pair consistent with both fitted confusion matrices is enumerated and
simulated separately, 200,000 draws each.

| Transition | Swapped | Observed | Scenarios | Pooled mean | Pooled sd | Pooled P(>=) | Central P(>=) | Min P(>=) | Median P(>=) | Max P(>=) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| sep15 -> oct03 | 17 | +.0556 | 110 | +.0270 | .0435 | .193 | **.238** | .000 | .224 | .260 |
| sep14 -> sep15 | 9 | +.0301 | 45 | +.0159 | .0317 | .455 | .427 | .376 | .453 | .502 |
| sep15 -> sep18 | 27 | -.1989 | 116 | -.0999 | .0512 | .981 | .997 | .958 | .981 | 1.000 |

The central scenario is the neutral composition in which the retained, removed
and added groups share one scored fraction and one precision among their scored
calls. For Sep 15 -> Oct 3 that is 9 scored positives, 11 scored negatives and
14 unscored among the 34 changed calls, with the removed half holding 4 true
positives and 7 false positives and the added half 5 and 4.

**Under Null B the observed gain sits at about the 76th percentile and is not
distinguishable from a coin-flip assignment of the changed calls.** The
probability of a move at least this large never exceeds .26 under any
consistent composition and is .19-.24 under the plausible ones; the only
scenarios giving a small p are the degenerate ones in which the 34 changed
calls contain a single scored positive. The parity of the arithmetic makes this
intuitive: a difference of exactly one true positive between two halves of a
random split is the median outcome whenever the changed set holds an odd number
of scored positives, and the extra information in the move comes only from the
three-false-positive difference, which is worth .009 MCC each.

The historical comparators calibrate this. The Sep 14 -> Sep 15 gain of +.0300
is at the 55th percentile of Null B, a coin flip. The Sep 15 -> Sep 18 loss of
-.1989 is at the 2nd percentile, clearly worse than an exchangeable reshuffle
and therefore a real failure. **On this scale the current gain resembles the
Sep 14 -> Sep 15 gain far more than it resembles the Sep 18 failure: the two
observed D6 gains are both inside the exchangeable null, while the one large
observed loss is outside it.**

## Part 3 - which score column produced the new calls

The L4 score is the equal percentile-rank mean of the deployed September 15
blind score and `pool_lgbm` from run `20261002T215345688167Z-7920b49a`. This run
recomputes both percentile ranks from the frozen score tables, asserts that the
deployed score reproduces the 49 September 15 calls and that the rank mean
reproduces the 49 October 3 calls, and checks every changed-call identity
against the recorded diagnostic `20261002T223515148996Z-ccb05cf3`:
`call_diagnostics.parquet` joins on 66 rows with 0 group mismatches and a
maximum absolute difference of exactly 0 in the legacy rank, pool rank,
combined rank and gate `p_active`, and `changed_calls.csv` lists exactly the
same 34 changed molecules with 0 group mismatches and 0 rank difference.
Nothing was recomputed from features. The recorded
`combined_rank` is the percentile rank of the rank mean, which orders the 750
rows identically to the rank mean itself.

| Group | n | Legacy rank | Pool rank | Rank mean | p_active | Neighbour positive rate | In `pool_lgbm` top 49 |
|---|---:|---:|---:|---:|---:|---:|---:|
| retained | 32 | .9681 | .9478 | .9579 | .3994 | .4661 | 22 |
| removed | 17 | .9659 | .6098 | .7879 | .1750 | .3483 | **0** |
| added | 17 | .9015 | .9530 | .9273 | .3094 | .3930 | **13** |

Per-molecule ranks for all 66 calls are in `attribution_calls.csv`. The
structure is unambiguous: the removed calls are high under the deployed score
(mean rank .966, essentially the same as the retained calls) and low under
`pool_lgbm` (.610, mean position 293); the added calls are high under
`pool_lgbm` (.953) and only moderately high under the deployed score (.902,
mean position 74).

### Would `pool_lgbm` alone have produced a vector closer to the new one?

| Comparison | Overlap at 49 calls |
|---|---:|
| `pool_lgbm` top 49 against the Sep 15 vector | 22 |
| `pool_lgbm` top 49 against the Oct 3 vector | **35** |
| `pool_lgbm` top 49 against the 32 retained calls | 22 |
| Sep 15 against Oct 3 | 32 |
| Of the 17 additions, ranked top 49 by `pool_lgbm` alone | **13** |
| Of the 17 removals, ranked top 49 by `pool_lgbm` alone | **0** |

`pool_lgbm` alone is substantially closer to the new vector (35 shared) than to
the old one (22 shared). Thirteen of the seventeen additions were already in
its own top 49, and none of the seventeen removals were. The rank mean did not
invent a new ordering; it imported `pool_lgbm`'s top-of-list preferences into
the deployed vector, and the deployed score then filtered them: 13 of
`pool_lgbm`'s own top 49 entered, the other 14 were kept out, 10 calls that
`pool_lgbm` would have dropped were protected, and 4 of the 17 additions came
from outside `pool_lgbm`'s top 49 altogether.

The attribution is therefore clear in composition and silent on value. One
scored true positive moved. Which of the 17 additions it was, and whether it
was one of the 13 that `pool_lgbm` alone would also have called, cannot be
determined without the hidden labels.

## Part 4 - what the offline estimates predicted

The paired group bootstrap of run `20261002T221918697661Z-2d8d20f5`
([concordance analysis](2026-10-03-d6-gate-concordance.md)) gave L4 against the
L0 legacy control at the 6.53% operating point:

| Label set | Split family | Offline delta | 95% interval | Contains realised +.0556 |
|---|---|---:|---|---|
| test-like | umap | +.029 | [-.014, +.073] | yes |
| test-like | scaffold | +.054 | [+.014, +.096] | yes |
| official | umap | -.041 | [-.102, +.031] | no |
| official | scaffold | +.012 | [-.044, +.068] | yes |

The realised +.0556 is almost exactly the test-like scaffold point estimate of
+.054 and sits inside three of the four intervals, the exception being the
official-label umap one. That is a favourable outcome for the test-like
estimator, and it is the first prospective D6 result that the test-like
criterion called in the right direction. It is also, and this is the point,
an outcome that three of the four mutually contradictory offline estimates
would each have claimed as a success.

It does not establish calibration, and three specific reasons say so.

- One observation. A single transition cannot distinguish estimators; nothing
  here is a calibration curve.
- The Live measurement is less precise than the offline estimates it is being
  used to check. Its own uncertainty is .0875 reported bootstrap sd, or .0435
  under the Null B exchangeable distribution, while the whole spread between
  the competing offline predictions is .095 (from -.041 to +.054). The
  observation cannot resolve a difference smaller than the gap between the
  hypotheses.
- The two label sets are not being compared on equal footing. The offline
  test-like delta is an MCC over 4,190 pseudo-labelled molecules at about 400
  calls; the Live delta is an MCC over 26 scored positives at 28 scored calls
  under official labels. The Live board is official-label scoring, so the
  nominally relevant offline row is the official one, which predicted about
  zero. The test-like row that matched the outcome is measuring a different
  quantity and matched it by agreeing on a sign that both label sets could
  easily have produced by chance given Null B.

The honest summary is that the test-like estimator was not refuted and happened
to land on the realised value, while the official-label estimator was also not
refuted on scaffold. Three of four intervals covering the outcome is not
calibration evidence; it is a measurement too coarse to separate them. Treat
the match as encouraging and uninformative.

### The eighth concordance transition

The frozen concordance design of run `2d8d20f5` is applied unchanged to the new
Sep 15 -> Oct 3 transition: mean score of the 17 removed against the 17 added
calls, two-sided permutation over 10,000 relabelings of the union with a fixed
seed, sign agreement against the Live MCC delta of +.0556.

| Column | Mean removed | Mean added | Difference | Perm. p | Sign agrees | Selection component |
|---|---:|---:|---:|---:|---|---|
| gate `p_active` (Morgan) | .1750 | .3094 | +.1344 | .0017 | yes | no |
| gate `p_active` (CheMeleon) | .1752 | .3281 | +.1530 | .0055 | yes | no |
| `pool_lgbm` score | .0653 | .4096 | +.3443 | .0001 | yes | **yes** |
| `pool_chemeleon` score | .5325 | .4448 | -.0877 | .3089 | no | no |
| `rank_mean_gates` score | .8356 | .8258 | -.0099 | .8309 | no | no |
| gate `p_tdi` (official head) | .4749 | .3833 | -.0916 | .2476 | no | no |
| `official_lgbm` score | .4749 | .3833 | -.0916 | .2448 | no | no |

Updated tally over eight transitions:

| Column | Sign agreements | Before this transition | Agrees here | Transitions with perm. p < .05 |
|---|---:|---:|---|---:|
| gate `p_active` (Morgan) | **7 of 8** | 6 of 7 | yes | 4 |
| gate `p_active` (CheMeleon) | **7 of 8** | 6 of 7 | yes | 4 |
| `pool_lgbm` score | **7 of 8** | 6 of 7 | yes | 4 |
| `pool_chemeleon` score | 6 of 8 | 6 of 7 | no | 2 |
| `rank_mean_gates` score | 6 of 8 | 6 of 7 | no | 4 |
| gate `p_tdi` (official head) | 5 of 8 | 5 of 7 | no | 4 |
| `official_lgbm` score | 5 of 8 | 5 of 7 | no | 4 |

**The pseudo-negative-trained column and both gate columns go from 6 of 7 to
7 of 8; the two columns that now agree with the new transition and previously
agreed with six are `gate_x_tdi.p_active`, `gate_x_tdi_chemeleon.p_active` and
`pool_lgbm.score`.** `pool_chemeleon` and `rank_mean_gates` stay at 6 but now
out of 8, and the official-only head stays at 5 of 8.

Two warnings attach to this update and they are not small.

- `pool_lgbm` is one of the two components of the score that selected this
  payload, so its +.3443 separation on this transition is close to
  tautological: the added calls were chosen partly because `pool_lgbm` ranks
  them high and the removed calls partly because it ranks them low. Its tally
  should be read as 6 of 7 with one non-independent eighth.
- The gate columns are not components of L4 and their agreement here is a
  genuine out-of-sample observation, but the eighth transition carries the same
  weight as the others, which is to say very little: eight transitions against
  one hidden set, one scored mask and one bootstrap estimate, with Null B
  placing this particular direction at the 76th percentile of chance.

## Part 5 - what the next submission-worthy experiment should be

### The sample-size reality

Each D6 upload resolves 26 scored positives and 28-32 scored calls. One scored
true positive is worth .0399 MCC and one scored false positive .0085. A
candidate whose true effect is less than roughly one expected scored true
positive cannot be measured by one upload, and the Null B spread of .0435
means that even a one-true-positive effect is not individually resolved.

Two facts make probing cheap anyway. The board reports the latest submission
rather than the best, and the project has already restored a superseded best
once (the Sep 24 re-send of `20260924T050705396203Z-10e1d3ad` recovered .3086),
so a probe that loses costs one opportunity and one re-send rather than the
record. With one opportunity per twelve hours per track and a 2026-11-03
deadline there are on the order of sixty remaining classification
opportunities, so the binding constraint is which questions are worth asking,
not how many uploads are available.

The implication for design is blunt: **spend opportunities only on probes that
change enough calls to move the metric by more than one expected true positive,
and settle everything else offline.** A probe that changes fewer than about ten
calls cannot pay for an opportunity.

### Ranked by information gain per opportunity

**1. `pool_lgbm` alone at 49 calls.** The highest-information probe available.
It changes 14 of the current 49 calls (overlap 35). The two live hypotheses
predict materially different outcomes for it: the test-like estimator puts B
above L0 by +.068 [+.018, +.118] on scaffold and +.022 [-.028, +.077] on umap
against L4's
+.054 and +.029, while the official-label estimator puts B at -.030 and -.116
against L4's +.012 and -.041. Differencing those against L4 gives an expected
effect on the current best of -.007 to +.014 under the test-like hypothesis and
-.042 to -.075 under the official-label one, so the two hypotheses are .05 to
.09 apart, above the resolution. A loss of that size is itself the
informative outcome, because it would say the deployed legacy component is
carrying the vector and the pseudo-negative model is not. This is the only
candidate on the list whose two hypotheses differ by more than the measurement
resolution, and even here the separation is only about 1.4 Null B standard
deviations, so one upload shifts the balance between the two views without
settling it.

**2. The current L4 vector extended to 60 calls.** A pure addition: the L4 top
60 contains all 49 current calls plus 11, so nothing is removed and the
experiment is monotone. It directly measures the marginal precision of ranks
50-60, which is the open question left by the rate-curve work, with no
confounding recipe change. About 6 of the 11 added rows will be scored; the
expected effect is +.027 if the marginal precision matches the current .39,
-.009 at a marginal precision near .17, and -.046 if none of the additions is a
positive. The spread across those scenarios, .073, exceeds the resolution, so
the probe is readable. It also tests the test-like criterion's clearest and
least label-dependent prediction, that the rate should be 7-10%.

**3. Rank-mean weight variation.** Reweighting the mean towards either
component moves few calls, since the two components already agree on the 32
retained ones. Any weight between the extremes produces a vector between the L4
and `pool_lgbm` vectors that experiments 1 and 2 already bracket. Expected
effect below the one-true-positive resolution. Settle offline on the frozen
score tables; do not spend an opportunity.

**4. Tightening or loosening the pseudo-negative log2fc cut.** Already answered
offline and negatively. Part 4 of the concordance run rebuilt the pool with the
cut at log2fc >= 0 and found the variant ordering unchanged (Spearman .983,
adjacent swaps only), nested pooled MCC moving by at most .003 for every
variant but one, and every bootstrap interval keeping its sign and endpoints to
within .005. A blind vector built from a model refitted at a different cut
would change very few of the 49 calls. Expected Live effect far below
resolution. This should never consume an opportunity.

### What to do instead of another upload

The limiting factor is not the number of opportunities but the information per
opportunity, and 26 scored positives will not change before the deadline. The
higher-value work is offline: widen the rate grid past the 10% boundary that
currently truncates every selection, and build an evaluation whose label set is
not the thing under test, because the whole pseudo-negative question is
currently being judged by a criterion that assumes its own answer.

## What this does and does not establish

Established by this run:

- One scored mask and one label set underlie all eight saved D6 Live states;
  the prevalence identity returns .0706 from every one of them, and a joint
  integer fit over all eight gives 374 scored rows of 750 and 26 scored
  positives with every one of the ten best (N, P) pairs agreeing on P = 26 and
  on the per-state true-positive counts.
- The Sep 15 -> Oct 3 move is +1 true positive and -3 false positives, stable at
  +0.9 to +1.2 and -2.7 to -3.5 over scored sizes from 330 to 430. One true
  positive accounts for about 70% of the headline .0556.
- The swap was not unskilled: under a draw from the uncalled pool at the base
  rate, P(move at least this large) = .0018.
- The swap was not demonstrably skilled either: under an exchangeable split of
  the 34 changed calls, P(move at least this large) is .19-.26 under every
  consistent composition, about the 76th percentile. The Sep 14 -> Sep 15 gain
  is at the 55th percentile and the Sep 15 -> Sep 18 loss at the 2nd.
- `pool_lgbm` alone at 49 calls is closer to the new vector (35 shared) than to
  the old one (22 shared); 13 of the 17 additions and 0 of the 17 removals are
  in its own top 49. The additions are `pool_lgbm`'s preferences, filtered by
  the deployed score.
- The gate columns and `pool_lgbm` now agree with 7 of 8 Live directions;
  `pool_lgbm`'s eighth agreement is not independent of the payload that
  produced it.

Not established:

- That the recipe change caused the gain. The gain is inside the exchangeable
  null, so a chance ordering of the 34 changed molecules is a sufficient
  explanation.
- That the offline test-like estimator is calibrated. It matched once, on an
  observation whose own uncertainty exceeds the spread between the competing
  predictions.
- Anything about which molecules are scored or labelled. The 374 and 26 are
  reconstructions from aggregates; the identities remain hidden.
- Any effect size for a future submission. The per-upload resolution is about
  .04 MCC per scored true positive against a .044 null spread.

## Limitations of the procedure itself

- The reconstruction assumes the published standard deviations come from a
  plain bootstrap of the scored rows. A stratified or smoothed bootstrap would
  bias the variance inversions, though not the prevalence identity.
- The integer fit weights the five metrics equally in units of their own
  published standard deviations and does not model the correlation between
  them, so its loss is a fit criterion, not a likelihood.
- Null B enumerates compositions rather than inferring them, and reports a
  range instead of a posterior. The range is narrow (.19 to .26 excluding the
  degenerate single-positive scenarios) so no prior is needed, but the number
  reported is not a p-value from a single well-specified model.
- Point MCC deltas are used throughout in place of bootstrap-mean deltas, on
  the evidence of a 2,000-draw check whose mean discrepancy is +.0031.
- The eight Live states are not eight independent experiments. They share one
  hidden set and consecutive payloads overlap heavily.
- Run `20261003T011015661743Z-f86e4d4d` is the superseded first execution of
  this analysis. Its outputs are byte-identical to the retained run except for
  the four `changed_calls.csv` identity fields added to `summary.json`; it is
  preserved as a run record but not promoted.
