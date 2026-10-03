# Validating the pseudo-negative device without using it as the yardstick

Four recorded classification runs. None fitted a blind vector, created a candidate,
prepared a submission, or contacted the challenge API.

| Part | Run | Seconds | Code |
| --- | --- | ---: | --- |
| 1 Direct rate measurement | `20261003T023544961394Z-95af3eea` | 2 | `scripts/evaluate_cyp_pseudo_negative_rate.py`, `scripts/launch_cyp_pseudo_negative_rate.py`, `configs/classification/pseudo-negative-rate-v1.json` |
| 2 Induced defect and fix | `20261003T024225174508Z-eff829e4` | 300 | `scripts/evaluate_cyp_induced_defect.py`, `scripts/launch_cyp_induced_defect.py`, `configs/classification/induced-defect-v1.json` |
| 3 Non-presupposing D6 metric | `20261003T024837735720Z-008f8427` | 180 | `scripts/evaluate_cyp_d6_official_metric.py`, `scripts/launch_cyp_d6_official_metric.py`, `configs/classification/d6-official-metric-v1.json` |
| 4 Probe pre-registration | `20261003T025354186349Z-621b6ad3` | 1 | `scripts/evaluate_cyp_d6_probe_preregistration.py`, `scripts/launch_cyp_d6_probe_preregistration.py`, `configs/classification/d6-probe-preregistration-v1.json` |

Outputs with SHA-256 provenance are promoted to
`artifacts/tracked/cyp-pseudo-negative-rate-20261003T023544961394Z-95af3eea/`,
`artifacts/tracked/cyp-induced-defect-20261003T024225174508Z-eff829e4/`,
`artifacts/tracked/cyp-d6-official-metric-20261003T024837735720Z-008f8427/` and
`artifacts/tracked/cyp-d6-probe-preregistration-20261003T025354186349Z-621b6ad3/`.
Three superseded executions are preserved with recorded reasons and not promoted:
`20261003T024611159208Z-ef84bcde` (Part 3, crossover table contaminated by the
degenerate mixture fraction and no deployed-recipe analogue),
`20261003T025124159090Z-2aa7202b` and `20261003T025205009648Z-e55e56ba` (Part 4,
wrong prediction baseline and missing marginal-precision table).

Retention script: `scripts/retain_cyp_pseudo_negative_validation.py`.

## The circularity being resolved

All 1,497 supplied `CYP2D6_is_TDI` values belong to rows that were hits in the
single-concentration primary screen; the largest labelled D6 `log2fc_estimate` is
`-0.592`. The 750 blind rows were measured regardless of that screen. Earlier work
([test-like pool](2026-10-03-d6-testlike-pool.md)) therefore relabelled unlabelled rows
with `log2fc >= -0.25` as negatives and compared recipes on the enlarged pool. Pool-trained
recipes won there and lost on the officially labelled subset alone. That comparison cannot
settle the question, because the pseudo-negative labels are both the hypothesis and the
yardstick.

The three parts below each avoid that. Part 1 measures the assumption where labels exist.
Part 2 manufactures the D6 defect on enzymes that have those labels and judges the repair
against labels no arm saw. Part 3 rescores the frozen D6 predictions under a metric whose
composition is identified from the published aggregates rather than assumed.

## Part 1 — the assumption, measured directly

### Verification of the supplied starting facts

Every stated count reproduces exactly. `data/raw/cyp-challenge-TRAIN_Emax.csv` carries
`is_TDI` for all four enzymes; the canonical `measurements.parquet` holds 1A2 1,413 /
2C9 1,285 / 2D6 1,497 / 3A4 3,584 of them. Joining each enzyme's labels to the median
50 uM `log2fc_estimate` of that same enzyme (17,504 primary rows, 4,376 molecules per
enzyme, all development) gives, for labelled rows:

| Enzyme | Labelled (raw / canonical) | `log2fc >= -0.25`: n | positives | `[-0.5, -0.25)`: n | positives | max labelled `log2fc` |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| CYP1A2 | 1,414 / 1,413 | 130 | 0 | 122 | 1 | +2.703 |
| CYP2C9 | 1,285 / 1,285 | 161 | 0 | 249 | 3 | +0.233 |
| CYP2D6 | 1,498 / 1,497 | 0 | 0 | 0 | 0 | −0.592 |
| CYP3A4 | 3,584 / 3,584 | 347 | 0 | 178 | 1 | +0.272 |

No disagreement with the supplied table. The raw-versus-canonical difference is one 1A2
row and one 2D6 row, both without a primary measurement, consistent with the quarantined
structure; modelling uses the canonical table throughout. 1,495 of the 1,497 canonical D6
labels carry a D6 `log2fc`, so D6 really is confined to primary hits.

### The pooled rate and its exact bound

Pooling the three enzymes that have weak-primary labels:

| Stratum | n | positives | rate | 95% one-sided upper | 95% two-sided CI |
| --- | ---: | ---: | ---: | ---: | --- |
| `log2fc >= -0.25` (the relabelled region) | **638** | **0** | **.0000** | **.00468** | [.00000, .00577] |
| `[-0.5, -0.25)` (borderline, excluded by the device) | 549 | 5 | .00911 | .01905 | [.00296, .02112] |
| `log2fc < -0.5` (primary hits, the labelled region) | 3,327 | 662 | .19898 | .21071 | [.18553, .21296] |
| no primary measurement (3A4 only) | 1,768 | 158 | .08937 | .10132 | [.07648, .10364] |

Per enzyme in the relabelled region: 1A2 0/130, upper bound .0228; 2C9 0/161, upper .0184;
3A4 0/347, upper .0086. Clopper-Pearson throughout.

**The pseudo-negative assumption is correct, by direct measurement, to within half a
percent.** Among 638 labelled rows in exactly the region the D6 device relabels, across
three enzymes and two different assay readouts, there is not one positive, and the exact
upper bound on the positive rate is 0.47%. The device's own pool is 7.73% positive, so a
0.47% contamination ceiling is a factor of sixteen below the signal it is used to measure.

The adjacent band the device deliberately excludes is not clean: 5 of 549, 0.91%. Excluding
it was the right call and the margin is small. The missing-primary stratum, which the device
never relabels, is emphatically not negative: 158 of 1,768, 8.94%, with a confidence interval
nowhere near zero. Relabelling that stratum would be a different and much worse idea, and
Part 2 tests exactly that as a control.

### The cut as a point on a curve

Sweeping the cut from 0.0 down to −1.5 in steps of 0.05, pooled over 1A2, 2C9 and 3A4
(`cut_curve.csv` carries all 31 cuts per enzyme as well):

| Cut | n at or above | positives | rate | 95% upper |
| ---: | ---: | ---: | ---: | ---: |
| 0.00 | 114 | 0 | .0000 | .0259 |
| −0.10 | 212 | 0 | .0000 | .0140 |
| −0.20 | 473 | 0 | .0000 | .0063 |
| **−0.25** | **638** | **0** | **.0000** | **.0047** |
| −0.30 | 780 | 1 | .0013 | .0055 |
| −0.40 | 1,000 | 3 | .0030 | .0077 |
| −0.50 | 1,187 | 5 | .0042 | .0088 |
| −0.75 | 1,550 | 17 | .0110 | .0164 |
| −1.00 | 1,923 | 51 | .0265 | .0334 |
| −1.25 | 2,313 | 94 | .0406 | .0483 |
| −1.50 | 2,850 | 151 | .0530 | .0604 |

Differencing consecutive cuts gives the rate inside each 0.05-wide bin. Every bin from 0.00
down to −0.30 is empty of positives (114, 42, 56, 89, 172 and 165 rows in turn, zero positives
in all six); the first positive appears in `[-0.35, -0.30)`; the rate then climbs monotonically
to roughly 1% in `[-0.5, -0.3)`, 5% near −0.75, 13% near −1.0 and 10% past −1.4.

The chosen cut of −0.25 therefore sits at the far end of a 638-row zero-positive plateau, with
the gradient beginning about 0.05 to 0.10 further down. The cut is neither tuned nor
precarious: anything from 0.0 to −0.25 gives the same answer with a monotonically tightening
bound, and the first erosion costs 0.1 percentage points. This also explains why the earlier
sensitivity check that moved the cut to `log2fc >= 0` changed nothing (RL0213, Spearman .983):
both cuts land on the same plateau.

### How far the result transfers to D6, stated honestly

Three enzymes validate the rule for themselves. What carries to D6 and what does not:

**What carries.** The rule being tested is an inference about the *assay protocol*, not about
chemistry: a molecule that did not perturb the single-concentration primary readout was not
followed into dose response, and where it was followed anyway it was assigned negative. That
protocol logic is shared across the four enzymes, which were screened on the same plates at
the same 50 uM concentration from the same batches. The near-zero rate is reproduced
independently three times, at sample sizes of 130, 161 and 347, with consistent direction and
a pooled bound of 0.47%. The three enzymes have very different hit prevalences in the labelled
region (1A2 1.9%, 2C9 4.0%, 3A4 46.9%), so the result is not an artifact of one enzyme's base
rate.

**What does not carry.** The positive rate approaching the cut *from below* behaves differently
for D6 than for the other three. At `log2fc >= -0.75`, pooled 1A2/2C9/3A4 is 17 of 1,550
(1.1%), whereas D6 is 23 of 86 (26.7%); at `log2fc >= -1.5` the pooled figure is 5.3% and D6
is 24.1%. D6's positive rate does not collapse as it approaches the boundary of its observed
support, while the other enzymes' does. Two readings are possible and the data cannot separate
them: either D6 TDI genuinely extends further toward weak primary activity, in which case the
relabelled region could hold positives that the three-enzyme bound would not catch, or the
labelled D6 region is simply richer throughout because D6 labelling was conditioned on a
stricter hit criterion. The absence of any D6 labels at or above −0.59 means the extrapolation
is across a gap of about 0.33 log2 units, not a point estimate at the cut.

**What would break it.** A D6-specific mechanism that produces time-dependent inhibition without
perturbing the primary readout at 50 uM — for instance a readout interference that cancels in the
single-concentration arm — would make the relabelled rows contaminated at a rate the other
enzymes cannot bound. Nothing in the available columns detects that. The 58 D6 primary hits with
no label are the visible edge of a selection rule we have not fully reconstructed, and the two
labelled D6 rows with no primary measurement show the join itself is not complete.

The honest summary: the assumption is measured, not assumed, and it holds to 0.47% on three
enzymes; its transfer to D6 rests on a shared protocol argument plus a 0.33-unit extrapolation
that D6's own near-boundary rate makes less comfortable than the three-enzyme bound alone
suggests.

### How many hidden D6 positives would overturn the ranking

The partial-label bound machinery of run `20261002T085410573284Z-e4d60d67` returns fully
adversarial feasible intervals, which for this population are uninformative by construction
(its D6 complete-pool rows span −0.41 to +0.42). It was therefore not reused as-is. The
sharper question is parameterised by count: holding the frozen held-fold scores fixed, how
many of the 2,693 D6 pseudo-negatives would have to be truly positive for
`MCC(pool_lgbm) − MCC(official_lgbm)` to stop being positive, and where would they have to sit?

Each relabelled row moves `FP -> TP` for a model that calls it and `TN -> FN` for one that does
not, so the delta is a function of how the hidden positives distribute over the four call cells
(called by both, control only, pool only, neither). Two allocations are reported: a grid-searched
adversarial minimum, and a proportional allocation in which the hidden positives fall across the
cells in proportion to cell size.

| Split | Rate | Observed delta | Cells (both / control only / pool only / neither) | Adversarial k | as a rate | Proportional k | as a rate |
| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| umap | 6.53% | +.0614 | 27 / 130 / 49 / 2,487 | 18 | 0.67% | none at k ≤ 400 | > 14.9% |
| umap | 10% | +.0603 | 58 / 205 / 72 / 2,358 | 20 | 0.74% | 394 | 14.6% |
| scaffold | 6.53% | +.1126 | 36 / 129 / 53 / 2,475 | 30 | 1.11% | none at k ≤ 400 | > 14.9% |
| scaffold | 10% | +.1074 | 71 / 183 / 72 / 2,367 | 36 | 1.34% | none at k ≤ 400 | > 14.9% |

Read against the Part 1 bound of 0.47%:

- Under any allocation that is not chosen by an adversary, the ranking is safe by a wide
  margin. The required rate is 14.6% or more, thirty-one times the measured upper bound and
  nearly twice the pool's own positive rate. On three of four cells the ranking does not
  overturn at all within the searched range.
- Under a maximally adversarial allocation it is not safe with much margin at all. As few as
  18 hidden positives, 0.67% of the pseudo-negatives, placed entirely in the cell the control
  calls and the pool model does not, would erase the umap advantage. The pooled 0.47% bound
  excludes that, but only just, and the weakest per-enzyme bound (1A2, 2.28%) does not exclude
  it. The scaffold advantage needs 1.11% to 1.34% and is correspondingly more robust.

The adversarial case is not a physically motivated scenario — it requires the hidden positives
to be concentrated in the 129 to 205 rows the control happens to prefer — but it is the honest
worst case, and it is the reason the margin is reported both ways rather than only as the
reassuring proportional figure.

## Part 2 — induce the defect, then test the fix

### Design and two changes to the specified plan

Inside the fixed five-fold `umap` and `scaffold` assignments
(`20260909T135130268493Z-edba6fc0`), for each of 1A2, 2C9 and 3A4, labels whose primary
`log2fc >= -0.5` or is missing are hidden from fitting and from rate selection. The
relabelled rows are then added back as pseudo-negatives and both arms are judged against the
real, complete held-fold labels, which only the ceiling arm ever saw. One recipe is used for
every arm: the balanced Morgan r2/2048 plus 217 RDKit descriptor LightGBM of the D6 device,
400 trees, two seeds averaged. Calls come from every development molecule of the held fold.

Two changes were made to the plan as specified, both because the plan as written would have
answered a weaker question.

**Change 1: the dose was wrong.** The specified fix re-adds only the *labelled* weak-primary
rows, which is 130 rows against 1,161 visible for 1A2 and 347 against 1,291 for 3A4 — a dose
of 0.11 to 0.39 pseudo-negatives per visible row. D6 receives 2,693 against 1,495, a dose of
1.80. A fix tested at a fifth of the dose cannot measure the effect at full dose. The arm
`fixed_full` therefore applies the D6 device verbatim: every weak-primary development row,
labelled or not, enters as a pseudo-negative. That gives doses of 1.54 (1A2), 1.23 (2C9) and
0.39 (3A4), bracketing D6's 1.80 on two of three enzymes. The narrower `fixed` arm is kept so
that the full-dose and labelled-only versions can be separated.

**Change 2: one ceiling is not enough.** The specified ceiling trains on all real labels, which
for 3A4 includes 1,768 missing-primary rows holding 158 positives. The pseudo-negative device
cannot supply those by construction, so measuring the fix against that ceiling guarantees a
small recovery fraction for reasons unrelated to the device. `ceiling_screened`, which trains on
the real labels of every *screened* row, bounds what a weak-primary repair could possibly close.
Both are reported.

A fourth arm, `fixed_missing`, additionally relabels the missing-primary rows False. Part 1
measured that stratum at 8.94% positive, so this arm is wrong by construction; it exists to
check that the metrics respond to label quality rather than to population coverage alone. It is
defined only for 3A4.

The deployed A4 equal-seven ensemble is not refitted here. Doing so would cost most of the
runtime budget and would not sharpen the defect-fix-ceiling contrast, which needs matched arms
rather than the strongest available model. That is a deliberate reduction, recorded as such.

### Label census

| Enzyme | labelled | positives | strong (`< -0.5`) | borderline | weak labelled | weak unlabelled | missing-primary | device dose | hidden from defect | hidden positives |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CYP1A2 | 1,413 | 23 | 1,161 (22 pos) | 122 (1) | 130 (0) | 1,659 | 0 | 1.54 | 252 | 1 |
| CYP2C9 | 1,285 | 38 | 875 (35) | 249 (3) | 161 (0) | 911 | 0 | 1.23 | 410 | 3 |
| CYP3A4 | 3,584 | 764 | 1,291 (605) | 178 (1) | 347 (0) | 161 | 1,768 (158) | 0.39 | 2,293 | 159 |

The 3A4 visible set of 1,291 reproduces the `known_n = 1291` of run
`20261002T042725116238Z-fe96ab3c` exactly, confirming the two experiments mask the same rows.
Note what the induced defect removes: for 1A2 and 2C9 it hides 252 and 410 rows holding 1 and 3
positives, so it is almost purely a *negative*-coverage defect, which is exactly D6's situation.
For 3A4 it also hides 158 real positives, because 3A4's unscreened rows were labelled whereas
D6's were not — that stratum is D6-like in being unlabelled, but the information lost is larger.

### Judged on real complete labels, strong-dominated population

Nested rate selected on the other four folds using only each arm's own visible labels, MCC on
the held fold against all real labels, pooled over folds:

| Enzyme | Split | defect | fixed | fixed_full | ceiling_screened | ceiling | fixed_missing |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CYP1A2 | umap | .1129 | .1100 | .0603 | .1636 | .1636 | – |
| CYP1A2 | scaffold | .0519 | .0595 | .0553 | .0646 | .0646 | – |
| CYP2C9 | umap | .0742 | .0249 | .0394 | .1030 | .1030 | – |
| CYP2C9 | scaffold | .0522 | .0674 | .0167 | .0484 | .0484 | – |
| CYP3A4 | umap | .0890 | .1250 | .1316 | .1069 | .3125 | .3303 |
| CYP3A4 | scaffold | .0945 | .1251 | .1310 | .1238 | .3949 | .3397 |

Paired group bootstrap against the defect arm, 1,000 resamples of assignment groups, real
labels, nested operating point:

| Enzyme | Split | arm | delta | 95% interval |
| --- | --- | --- | ---: | --- |
| CYP1A2 | umap | fixed | −.0029 | [−.0592, +.0730] |
| CYP1A2 | umap | fixed_full | −.0527 | [−.1450, +.0789] |
| CYP1A2 | umap | ceiling_screened | +.0507 | [−.0013, +.1297] |
| CYP1A2 | scaffold | fixed | +.0076 | [−.0442, +.0666] |
| CYP1A2 | scaffold | fixed_full | +.0034 | [−.0492, +.0649] |
| CYP1A2 | scaffold | ceiling_screened | +.0127 | [−.0408, +.0738] |
| CYP2C9 | umap | fixed | **−.0493** | **[−.0969, −.0092]** |
| CYP2C9 | umap | fixed_full | −.0348 | [−.1122, +.0706] |
| CYP2C9 | umap | ceiling_screened | +.0288 | [−.0161, +.0940] |
| CYP2C9 | scaffold | fixed | +.0152 | [−.0620, +.0892] |
| CYP2C9 | scaffold | fixed_full | −.0355 | [−.1185, +.0562] |
| CYP2C9 | scaffold | ceiling_screened | −.0038 | [−.0927, +.0837] |
| CYP3A4 | umap | fixed | **+.0360** | **[+.0023, +.0672]** |
| CYP3A4 | umap | fixed_full | **+.0425** | **[+.0099, +.0744]** |
| CYP3A4 | umap | ceiling_screened | +.0178 | [−.0173, +.0536] |
| CYP3A4 | umap | ceiling | +.2234 | [+.1560, +.2862] |
| CYP3A4 | scaffold | fixed | +.0305 | [−.0003, +.0616] |
| CYP3A4 | scaffold | fixed_full | **+.0365** | **[+.0027, +.0685]** |
| CYP3A4 | scaffold | ceiling_screened | +.0293 | [−.0070, +.0620] |
| CYP3A4 | scaffold | ceiling | +.3004 | [+.2602, +.3422] |

On this population the result is: 3A4 positive and resolved, 1A2 and 2C9 null. On 3A4 both fix
arms exceed the screened ceiling, so the fix recovers more than 100% of the gap a weak-primary
repair could close, on both split families, with intervals excluding zero for `fixed_full`. On
1A2 and 2C9, with 23 and 38 positives respectively, every interval except one spans zero and
the one exception (2C9 umap `fixed`, −.049) is negative. Those two enzymes do not have the
power to resolve a .03 effect and should not be read as evidence either way; reported plainly
because the plan asked for them.

### The population the table above uses is not the one that matters

The broad real-label population is strong-primary dominated: 1,161 of 1,413 rows for 1A2 and
1,291 of 1,816 screened rows for 3A4. The D6 blind scored population is not. Part 3 identifies
its strong-primary-like share at about one third. Measuring a device whose whole purpose is to
place weak-primary rows correctly on a population that is 82% to 92% strong-primary
systematically understates it.

The same arms were therefore also scored with the real-label strong and weak strata reweighted
to a target strong fraction `phi`, keeping every label real. At `phi = 1` this reduces to the
strong-only metric; `phi = 1/3` matches the identified blind composition. MCC at a fixed rate,
pooled over held folds:

| Enzyme | Split | Rate | defect | fixed | fixed_full | ceiling_screened | ceiling | fixed_missing |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CYP1A2 | umap | 6.53% | .0344 | .0333 | .0381 | .0359 | .0359 | – |
| CYP1A2 | scaffold | 6.53% | .0655 | .0560 | .0372 | .0740 | .0740 | – |
| CYP2C9 | umap | 6.53% | .0110 | .0093 | .0562 | .0149 | .0149 | – |
| CYP2C9 | scaffold | 6.53% | .0422 | .0131 | .0540 | .0284 | .0284 | – |
| CYP3A4 | umap | 6.53% | .1135 | .1908 | **.2181** | .1826 | .2668 | .2399 |
| CYP3A4 | scaffold | 6.53% | .1292 | .1921 | **.2072** | .1815 | .2575 | .2314 |
| CYP3A4 | umap | 20% | .1557 | .2907 | **.3228** | .3202 | .3101 | .3282 |
| CYP3A4 | scaffold | 20% | .2333 | .3144 | **.3241** | .3156 | .3306 | .3121 |

Paired group bootstrap at `phi = 1/3`, nested operating point, real labels:

| Enzyme | Split | arm | delta | 95% interval |
| --- | --- | --- | ---: | --- |
| CYP3A4 | umap | fixed | +.2002 | [+.1607, +.2363] |
| CYP3A4 | umap | fixed_full | **+.2197** | **[+.1721, +.2617]** |
| CYP3A4 | umap | ceiling_screened | +.2015 | [+.1528, +.2471] |
| CYP3A4 | umap | fixed_missing | +.1390 | [+.0986, +.1774] |
| CYP3A4 | scaffold | fixed | +.1378 | [+.0915, +.1843] |
| CYP3A4 | scaffold | fixed_full | **+.1848** | **[+.1361, +.2323]** |
| CYP3A4 | scaffold | ceiling_screened | +.1784 | [+.1300, +.2274] |
| CYP3A4 | scaffold | fixed_missing | +.0920 | [+.0401, +.1457] |
| CYP2C9 | umap | fixed_full | +.0008 | [−.0604, +.0700] |
| CYP2C9 | scaffold | fixed_full | −.0122 | [−.0710, +.0627] |
| CYP1A2 | umap | fixed_full | −.0093 | [−.0999, +.0830] |
| CYP1A2 | scaffold | fixed_full | +.0042 | [−.0308, +.0397] |

Effective sample size at `phi = 1/3` is 732 rows for 3A4, 346 for 2C9 and 285 for 1A2, against
nominal totals of 1,638, 1,036 and 1,291. Reweighting a 130-row stratum to two thirds of the
mass costs most of 1A2's sample, which is why 1A2 stays unresolved under both metrics.

**Three things follow.**

First, on the composition that matters, the verbatim device recovers the entire recoverable gap
on 3A4 and slightly exceeds it: `fixed_full` is +.220 and +.185 against the defect where
`ceiling_screened` is +.202 and +.178, and every interval excludes zero by a wide margin. The
defect is real, large, and repaired by relabelling.

Second, the deliberately wrong arm fails, and fails in the right place. On the strong-dominated
population `fixed_missing` looked excellent (+.245 and +.241, near the complete ceiling),
because relabelling 1,768 rows that are 91% negative buys most of the population-coverage
benefit regardless of the 158 errors. Under the composition-matched metric it is clearly worse
than `fixed_full` (+.139 against +.220 on umap, +.092 against +.185 on scaffold). The metric
distinguishes correct pseudo-negatives from incorrect ones; the broad metric does not. That is
the control working as intended, and it is also a warning about reading broad MCC gains as
evidence of label quality.

Third, 2C9 flips sign between the two metrics (`fixed_full` −.035 broad, +.001 composition-matched
on umap; −.036 broad, −.012 matched on scaffold) and 1A2 stays near zero in both. The two
low-positive enzymes contribute no evidence. The claim rests on 3A4.

### Reconciliation with run `20261002T042725116238Z-fe96ab3c`

That run (RL0205) masked weak and missing-primary 3A4 classes from fitting and from threshold
selection and reported strong-region AP improving `.47839 -> .53528` while broad-known MCC
reversed `.16330 -> .00568`, paired broad delta −.158 to −.160. The reading offered for
checking was that it established the defect rather than testing the fix. That reading is
confirmed, with one component that does not replicate and one that sharpens it.

**The broad MCC collapse replicates, larger.** Ceiling to defect on 3A4 broad real labels is
`.3125 -> .0890` on umap and `.3949 -> .0945` on scaffold, a loss of .223 and .300 against
fe96ab3c's .158. Different learner, same direction, same order of magnitude. The defect is real
and not specific to a torch multi-task network.

**The strong-region AP gain does not replicate.** fe96ab3c found the masked arm *better* on
strong-region AP by +.057. Here the masked arm is slightly *worse*: 3A4 strong-region AP is
.5683 (ceiling) against .5404 (defect) on umap, and .6014 against .5996 on scaffold. Strong-region
AUROC behaves the same way (.6127 against .6042; .6407 against .6450). The AP improvement in
fe96ab3c is therefore attributable to its own masking-plus-multi-task implementation rather than
to hiding weak labels as such, and should not be carried forward as a general property.

**The reversal is in the missing-primary stratum, not the weak-primary one.** This is the
substantive reconciliation and it changes the interpretation of RL0205. Splitting the
ceiling-to-defect gap by stratum: `ceiling_screened`, which keeps the weak and borderline labels
but drops the missing-primary ones, scores .1069 and .1238, barely above the defect's .0890 and
.0945, while the complete ceiling scores .3125 and .3949. On the missing-primary population
itself, the defect arm's AUROC is .502 and .504 — exactly chance — against the complete
ceiling's .809 and .890. So essentially the whole of the broad MCC reversal comes from losing
the 1,768 missing-primary labels and their 158 positives, and almost none of it from losing the
347 weak-primary labels that hold no positives at all.

That matters because the D6 pseudo-negative device touches only the weak-primary stratum. RL0205
is therefore not evidence against the device. It is evidence about a different stratum, the one
the device deliberately leaves alone and the one Part 1 measures at 8.94% positive. The two
results are consistent: hiding labels that contain positives is damaging (RL0205), and supplying
correct negatives where there are no positives to lose is helpful (this run). RL0205's statement
that "weak/missing-primary A4 has 159 known positives, disproving primary activity as a necessary
TDI condition" is reproduced exactly — 158 missing plus 1 borderline — and it is now resolved to
its stratum: 0 of 347 in the relabelled region, 1 of 178 in the excluded band, 158 of 1,768 in
the unscreened region.

## Part 3 — a D6 evaluation that does not presuppose the answer

Nothing is refitted. The frozen held-fold predictions of run
`20261002T215345688167Z-7920b49a` cover all 6,146 development molecules for seven recipes on
both split families, so every question below is a rescoring of the same scores under a different
metric. A deployed-recipe analogue was added: `l4_rank_mean` is the equal percentile-rank mean
of `legacy_sep15_refit` and `pool_lgbm` formed inside each held fold, which is the recipe of the
October 3 payload. Its blind reconstruction reproduces the deployed 49 calls exactly (49 of 49,
checked in Part 4), so the OOF analogue is the right recipe and not an approximation of it.

### Two metrics, and a mixture fraction that is identified rather than assumed

**Metric O**, as specified: score only officially labelled rows, with the two classes reweighted
so the scored prevalence equals a target prior. Explicitly, with `n` labelled rows of which `P`
are positive and a target prior `pi`, each positive row carries weight `pi * n / P` and each
negative row `(1 - pi) * n / (n - P)`, and MCC is computed from the weighted confusion counts.
Calls still come from every development molecule of the held fold. Metric O never credits a model
for a pseudo-negative row. It also never credits correct down-ranking of weak-primary rows,
because it does not contain any, so it is a lower bound on the pool model's value rather than a
neutral referee.

**Metric M**, added here: treat the blind scored set as a two-stratum mixture. The
strong-primary-like stratum is sampled by the officially labelled rows carrying their real
labels; the weak-primary-like stratum is sampled by the pseudo-negative rows. The strong stratum
receives mass `phi`. At `phi = 1` Metric M is the unweighted official metric. Part 1 is what
licenses this: the weak stratum's labels are not the model's hypothesis, they are a measured
quantity with a 0.47% error ceiling from three other enzymes, and the rows being scored are
held-fold rows that entered no fit.

`phi` is not chosen. The public aggregates pin it. With `p_strong = 324/1497 = .21643` observed
and `p_weak` bounded by Part 1, the prevalence identity `prevalence = phi * p_strong +
(1 - phi) * p_weak` inverts to:

| Scored prevalence | `p_weak = 0` | `p_weak = .00468` (pooled bound) | `p_weak = .00577` (CI upper) |
| ---: | ---: | ---: | ---: |
| 5% | .2310 | .2140 | .2100 |
| **7.06%** | **.3262** | **.3113** | **.3078** |
| 9% | .4158 | .4029 | .3998 |

Every combination lands inside `[0, 1]`, which is itself a consistency check: the published 7.06%
prevalence, the observed official-label prevalence and the measured weak-primary rate are mutually
compatible with a single mixture, and they pin that mixture at `phi = .308` to `.326`. The
development screened pool has `(1,495 + 58) / 4,376 = .355` primary hits, so the identified blind
composition is close to, and slightly below, the development composition — which is what analog
expansion around strong 1A2/2C9/3A4 hits would produce for a D6 endpoint. Across the 5% to 9%
prior range `phi` runs .21 to .42.

### Rate optima, with the grid extended to 50%

The grid runs 1% to 50% of the development pool, 22 values, well past the 10% boundary that
truncated every earlier selection. Optima at the 7.06% prior (Metric O) and the identified
`phi = .3262` (Metric M), with the equivalent blind call count:

| Split | Model | Metric O rate / calls / MCC | Metric M rate / calls / MCC |
| --- | --- | --- | --- |
| umap | `official_lgbm` | 12.5% / 94 / .0674 | 12.5% / 94 / .0439 |
| umap | `legacy_sep15_refit` | 8% / 60 / .1132 | 8% / 60 / .0767 |
| umap | `pool_lgbm` | 15% / 112 / .0462 | 25% / 188 / .1403 |
| umap | `l4_rank_mean` | **15% / 112 / .0816** | **15% / 112 / .1199** |
| scaffold | `official_lgbm` | 3% / 22 / .0776 | 2% / 15 / .0626 |
| scaffold | `legacy_sep15_refit` | 12.5% / 94 / .1388 | 12.5% / 94 / .1000 |
| scaffold | `pool_lgbm` | 7% / 53 / .0781 | 30% / 225 / .1657 |
| scaffold | `l4_rank_mean` | **27.5% / 206 / .1158** | **27.5% / 206 / .1476** |

**The interior optimum is real and it is above 10%, but it is shallow.** For 14 of the 16
model-by-split-by-metric cells the argmax is strictly inside the grid; only the two gate variants
on umap still select the 50% edge. The earlier grid-edge selections were an artifact of stopping
at 10%. The prior (5%, 7.06%, 9%) rescales Metric O's level but does not move a single argmax, so
the rate recommendation is insensitive to the prevalence assumption.

For the deployed recipe the two metrics agree on the optimum exactly — 15% on umap and 27.5% on
scaffold — which is a useful invariance, since it means the rate recommendation does not depend
on which label set one believes. But the curve is nearly flat between those points and 49 calls.
On Metric M, `l4_rank_mean` scores .1012 at 49 calls, .1190 at 94, .1199 at 112, .1081 at 150 and
.0906 at 206 on umap; .1384 at 49, .1428 at 94, .1336 at 112, .1273 at 150 and .1476 at 206 on
scaffold. The umap family peaks near 94 to 112 calls and then declines; the scaffold family is
bimodal with 49 already within .009 of its own best. The total span from 49 to the optimum is
+.019 (umap) and +.009 (scaffold), both below the Live resolution of about .04 per scored true
positive. **The interior optimum sits at 94 to 112 calls, and the offline evidence that moving
there helps is real but smaller than the public metric can measure.**

### Does the pool-trained recipe still beat the official-only control

Paired group bootstrap against `official_lgbm` at the historical 6.53% operating point, 1,000
group resamples:

| Metric | Split | `pool_lgbm` | `l4_rank_mean` | `legacy_sep15_refit` |
| --- | --- | --- | --- | --- |
| O, prior 7.06% | umap | −.0209 [−.0622, +.0256] | +.0300 [−.0077, +.0715] | +.0622 [+.0047, +.1217] |
| O, prior 7.06% | scaffold | +.0221 [−.0212, +.0656] | +.0537 [+.0145, +.0948] | +.0469 [−.0001, +.0964] |
| M, `phi = .3262` | umap | **+.0640 [+.0215, +.1100]** | **+.0683 [+.0304, +.1118]** | +.0375 [−.0084, +.0862] |
| M, `phi = .3262` | scaffold | **+.1137 [+.0698, +.1630]** | **+.0980 [+.0616, +.1385]** | +.0432 [+.0022, +.0826] |

Sign changes along `phi` at fixed rate, scanning the non-degenerate grid (`delta_vs_phi.csv`):

| Split | Model | delta at `phi = .3262` | delta at `phi = 1` | crossings |
| --- | --- | ---: | ---: | --- |
| umap | `pool_lgbm` | +.0640 | −.0312 | .850 |
| umap | `l4_rank_mean` | +.0683 | +.0444 | none |
| umap | `legacy_sep15_refit` | +.0375 | +.0850 | none |
| scaffold | `pool_lgbm` | +.1137 | +.0345 | none |
| scaffold | `l4_rank_mean` | +.0980 | +.0760 | none |
| scaffold | `legacy_sep15_refit` | +.0432 | +.0645 | none |

**The verdict is split by metric, and the split is informative rather than embarrassing.**

Under Metric O the answer to the question as posed is no: `pool_lgbm` does not beat the
official-only control. It is −.021 [−.062, +.026] on umap and +.022 [−.021, +.066] on scaffold,
both intervals spanning zero, and at each model's own optimal rate it scores .046 and .078
against the control's .067 and .078. Under the strong-primary-only metric the pool model buys
nothing, and the legacy September 15 refit remains the best recipe by a clear margin (.113 and
.139).

Under Metric M at the identified composition the answer is yes, decisively: +.064 [+.022, +.110]
and +.114 [+.070, +.163]. The crossover is at `phi = .85` on umap and does not occur at all on
scaffold within the grid. Since `phi = .8` would imply a scored prevalence of `.8 * .2164 = 17.3%`
against the published 7.06%, the composition at which the pool model stops winning is one the
published aggregates exclude. **The pool model's advantage over the official-only control is
robust across the entire range of mixture fractions compatible with the observed prevalence, and
disappears only at compositions that prevalence rules out.**

Those two statements are not in conflict. Metric O conditions on the strong-primary region and
therefore measures a different quantity: ranking skill among primary hits, where the pool model
has no advantage and arguably a small disadvantage from diluting its positive class. Metric M
measures skill on the population the board actually scores. The earlier "the direction of every
comparison depends on which label set is assumed correct" is now resolved in a specific way: the
two label sets correspond to two populations, the population is identified by the prevalence
identity, and at the identified population the pool-trained recipe wins. The remaining
qualification is that this is an identification of composition, not of which blind rows belong to
which stratum.

One further coherence point: `l4_rank_mean`, the deployed recipe, is the only model positive
against the control under *both* metrics and both splits (+.030 to +.098, three of four intervals
excluding zero). The combination is better behaved than either component, which is consistent
with the October 3 Live outcome without being evidence for it.

## Part 4 — pre-registration of the next probe

### What is left open, and the sample-size reality

Parts 1 to 3 settle the label question offline: the pseudo-negative labels are measurably
correct (0 of 638, bound 0.47%), supplying them repairs a real defect on the one enzyme with
power (+.185 to +.220 at the identified composition, intervals excluding zero), and the
pool-trained recipe beats the official-only control at the identified blind composition. None of
that needs an upload, and no recipe change is worth one.

What is left open is the call count. Once the grid is extended past 10% the optimum moves inside
to 94 to 112 calls under both metrics, and the two split families disagree about exactly where.
That is the only remaining question an upload can answer, and it is answerable because a pure
extension of the deployed vector changes nothing else.

From [the sizing analysis](2026-10-03-d6-live-gain-analysis.md): 374 scored rows of 750, 26
scored positives, the October 3 payload fitted at 28 scored calls with 11 true positives. One
scored true positive is worth +.0399; one added scored false positive −.0111 at this state
(RL0218's +.0083 per *removed* scored non-positive is the same nonlinearity seen from the other
side). The fitted matrix's point MCC is .36165 against the published bootstrap mean .36425, so
predictions below are deltas against the point value.

### Candidate vectors and their changed calls

Each candidate is the top `k` of a frozen score column. Changes are against the deployed October
3 payload (`20261002T223229553466Z-a08301d2`). Nothing was created, prepared or uploaded.

| Candidate | calls | added | removed | changed | monotone | expected added scored |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| `l4_rank_mean` top 49 | 49 | 0 | 0 | 0 | yes | 0.0 |
| `l4_rank_mean` top 60 | 60 | 11 | 0 | 11 | yes | 5.5 |
| `l4_rank_mean` top 75 | 75 | 26 | 0 | 26 | yes | 13.0 |
| **`l4_rank_mean` top 94** | **94** | **45** | **0** | **45** | **yes** | **22.4** |
| `l4_rank_mean` top 112 | 112 | 63 | 0 | 63 | yes | 31.4 |
| `l4_rank_mean` top 150 | 150 | 101 | 0 | 101 | yes | 50.4 |
| `pool_lgbm` top 49 | 49 | 14 | 14 | 28 | no | 7.0 |

The `pool_lgbm` top 49 figure of 14 added and 14 removed reproduces RL0218's "14 calls changed,
overlap 35" exactly. The reconstructed `l4_rank_mean` top 49 equals the deployed 49 calls
exactly, which validates the reconstruction used throughout Part 3.

### Public-metric consequence as a function of marginal precision

Predicted delta against the deployed .36165 point MCC, as a function of the marginal precision of
the added scored calls:

| Candidate | 0.00 | 0.07 | 0.10 | 0.17 | 0.25 | 0.30 | 0.393 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| top 60 | −.042 | −.028 | −.022 | −.008 | +.008 | +.018 | +.037 |
| top 75 | −.087 | −.057 | −.044 | −.013 | +.022 | +.044 | +.084 |
| **top 94** | **−.131** | **−.083** | **−.062** | **−.014** | **+.042** | **+.076** | **+.140** |
| top 112 | −.164 | −.101 | −.074 | −.011 | +.062 | +.107 | +.191 |
| top 150 | −.218 | −.127 | −.088 | +.003 | +.107 | +.172 | +.293 |

Break-even marginal precision is **.188** for top 94 and between .17 and .19 for every monotone
extension. Deployed precision is 11/28 = .393 over ranks 1 to 49; the base prevalence is .0706.

### Offline estimate of the marginal precision, with a calibration step

Held-fold precision of each rank band of `l4_rank_mean`, on officially labelled rows and under the
identified mixture (`marginal_precision.csv`):

| Split | Band | official rows | positives | official precision | pseudo rows | mixture precision |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| umap | 1–49 | 128 | 46 | .359 | 116 | **.176** |
| umap | 50–60 | 38 | 10 | .263 | 25 | .150 |
| umap | 61–75 | 29 | 11 | .379 | 39 | .149 |
| umap | 76–94 | 49 | 13 | .265 | 52 | .120 |
| umap | 95–112 | 39 | 12 | .308 | 59 | .112 |
| umap | 113–150 | 85 | 16 | .188 | 115 | .074 |
| scaffold | 1–49 | 129 | 55 | .426 | 108 | **.217** |
| scaffold | 50–60 | 24 | 7 | .292 | 27 | .127 |
| scaffold | 61–75 | 37 | 11 | .297 | 48 | .119 |
| scaffold | 76–94 | 47 | 17 | .362 | 64 | .141 |
| scaffold | 95–112 | 36 | 7 | .194 | 57 | .069 |
| scaffold | 113–150 | 82 | 20 | .244 | 115 | .093 |

Pooling bands 50–94 gives a mixture precision of .136 (umap) and .131 (scaffold), mean **.133**.
The same estimator applied to the deployed band 1–49 gives .176 and .217, mean **.197**, against
a realised Live precision of .393. **The offline mixture estimator understates the realised
precision of the deployed calls by a factor of 1.99.** That single calibration point is the only
one available, and it is one observation, but it is directly relevant: it is the same estimator
applied to the same recipe on the same board.

Two routes therefore bracket the prediction for ranks 50–94:

- Uncalibrated, marginal precision .133, predicted delta **−.039**.
- Calibration-scaled by 1.99, marginal precision .265, predicted delta **+.052**.

### The pre-registered prediction

Fixed before any upload exists, for `l4_rank_mean` top 94 against the current D6 Live MCC of
.364248:

- **Primary prediction.** Realised D6 Live MCC delta of **+.005**, i.e. a realised MCC of about
  **.369**, with a pre-registered 80% interval of **[−.080, +.080]** (MCC .284 to .444),
  corresponding to a marginal precision of ranks 50–94 between .10 and .31.
- **Implied marginal precision.** The realised delta identifies it: delta `-.131` means .00,
  `-.062` means .10, `-.014` means .17, `+.042` means .25 and `+.140` means .393.
- **Decision rules.** A delta above **+.080** establishes that the marginal precision of ranks
  50–94 is at or above .31 and that the interior optimum found in Part 3 transfers; raise the
  deployed count. A delta below **−.080** establishes a marginal precision at or below .07, i.e.
  ranks 50–94 are at the blind base rate, and refutes the extension hypothesis; restore 49 calls
  and stop asking. A delta inside `[-.080, +.080]` means the marginal precision is between .07 and
  .31 and one upload cannot localise it further; do not spend a second opportunity on the same
  question.
- **Also tested by the same observation.** The 1.99 calibration factor between the offline mixture
  estimator and the realised Live precision. If the realised delta lands near the calibrated
  prediction of +.052, the factor transfers out of sample; if it lands near the uncalibrated −.039,
  the factor was an artifact of ranks 1–49.
- **What this probe does not test.** It changes no recipe, so it says nothing further about the
  pseudo-negative device; Parts 1 to 3 have already settled that offline as far as offline evidence
  can. It also cannot distinguish the umap preference (112 calls) from the scaffold one (206); 94 is
  chosen because both families are within .009 of their own best there and the downside is smaller.

**Expected value is approximately zero and the information gain is the point.** The honest
statement is that this probe is predicted to move the board by less than the noise, while
identifying a quantity — the marginal precision of ranks 50 to 112 — that no amount of offline
work can pin down, because every offline estimate of it depends on the label set under test. The
conservative alternative, top 60, has a span of only .079 across scenarios against top 94's .271,
so it costs the same opportunity for a third of the information.

Ranked against the alternatives in RL0218: `pool_lgbm` alone at 49 calls remains the better probe
if the open question were still *which recipe*, but Part 3 has now answered that offline at the
identified composition, so its information per opportunity has dropped. Rank-mean weight variation
and pseudo-negative cut changes remain below resolution and are now also unnecessary: Part 1's cut
curve shows the whole plateau from 0.0 to −0.25 is equivalent.

No candidate, prepared submission or upload was created by any part of this work, and the
challenge API was not contacted.

## Design changes made, and why

1. **The Part 2 fix was tested at the wrong dose, so a full-dose arm was added.** The specified
   fix re-adds only labelled weak-primary rows, a dose of 0.11 to 0.39 pseudo-negatives per
   visible row against D6's 1.80. `fixed_full` applies the device verbatim and brackets D6's dose
   on two of three enzymes. The narrower arm is retained for separation.
2. **A screened-only ceiling was added.** The specified ceiling includes the missing-primary
   labels, which the device cannot supply, so recovery against it is bounded below 1 for
   structural reasons. `ceiling_screened` bounds what a weak-primary repair could close, and the
   difference between the two ceilings turned out to be the whole reconciliation with fe96ab3c.
3. **A deliberately wrong arm was added as a metric control.** `fixed_missing` relabels the 8.94%
   positive missing-primary stratum. It looks excellent on the broad metric and clearly worse on
   the composition-matched one, which is what demonstrates the composition-matched metric responds
   to label correctness rather than population coverage.
4. **The composition-matched metric replaced the specified prevalence reweighting as the primary
   Part 3 estimator, and the mixture fraction was identified rather than assumed.** Reweighting
   classes inside the official stratum changes the positive mass but leaves the negative population
   strong-primary-only, which is not the blind population. Mixing the two strata and pinning the
   mixture with the published prevalence identity gives `phi = .308` to `.326` from three
   independent quantities, turns the label-set conflict into a statement about population, and
   makes the sensitivity analysis (priors 5%, 7%, 9% give `phi` .21 to .42) a statement about
   composition rather than about an arbitrary weight. Metric O is still reported exactly as
   specified, and is explicitly labelled a lower bound rather than a neutral referee.
5. **The Part 3 comparison was extended to the deployed recipe.** `l4_rank_mean`, the OOF analogue
   of the October 3 payload, was absent from the frozen bank. It is reconstructible from two frozen
   columns and its blind form reproduces the deployed 49 calls exactly, so Part 4 could be
   pre-registered against the actual deployed vector instead of a proxy.
6. **The partial-label bound machinery of `e4d60d67` was not reused.** Its fully adversarial
   intervals for this population are uninformative by construction (−.41 to +.42 on the D6
   complete pool). The overturn question was re-derived parameterised by count and by allocation,
   which gives a usable number (0.67% adversarial, 14.6% proportional) instead of a vacuous range.
7. **Part 3 and Part 4 fit nothing.** Both rescore frozen predictions, which kept total compute
   under ten minutes and left the seed and fold design of Part 2 untouched rather than reduced.
8. **The A4 equal-seven ensemble was not refitted in Part 2.** Matched arms matter more than the
   strongest model for a defect-fix-ceiling contrast, and refitting it would have consumed the
   budget. Recorded as a reduction, not an oversight.

## What this does and does not establish

Established:

- The pseudo-negative assumption is correct where it is measurable: 0 positives among 638 labelled
  rows at `log2fc >= -0.25` across 1A2, 2C9 and 3A4, exact 95% one-sided upper bound 0.47%,
  two-sided CI [0, 0.58%]. Per enzyme 0/130, 0/161 and 0/347.
- The cut of −0.25 sits at the end of a 638-row zero-positive plateau; the first positive bin is
  `[-0.35, -0.30)` and the rate reaches 0.42% at −0.5, 2.65% at −1.0 and 5.30% at −1.5.
- The stratum the device excludes is 0.91% positive (5 of 549); the stratum it never touches is
  8.94% positive (158 of 1,768). Both exclusions are correct and the first is a narrow margin.
- Under a proportional allocation, 14.6% of the D6 pseudo-negatives would have to be positive to
  overturn the pool-versus-control ranking, thirty-one times the measured bound. Under a maximally
  adversarial allocation, 0.67% (18 rows) suffices on umap and 1.11% on scaffold, which the pooled
  bound excludes but the weakest per-enzyme bound does not.
- Inducing the D6 coverage defect on 3A4 costs .223 to .300 broad MCC against the complete ceiling,
  reproducing RL0205's reversal with a different learner and a larger magnitude.
- Re-adding the weak-primary rows as pseudo-negatives recovers the entire recoverable gap on 3A4:
  +.036 to +.043 broad and +.185 to +.220 at the identified composition, against a screened ceiling
  of +.018 to +.029 broad and +.178 to +.202 matched, with intervals excluding zero on both splits.
  **The pseudo-negative device is validated by evidence that does not use it as the yardstick.**
- That validation rests on 3A4 alone. 1A2 and 2C9, with 23 and 38 positives, produce intervals
  spanning zero in eleven of twelve comparisons and one significantly negative result (2C9 umap
  `fixed`, −.049 [−.097, −.009]).
- RL0205's A4 reversal is located in the missing-primary stratum, not the weak-primary one:
  `ceiling_screened` is only +.018 and +.029 above the defect while the complete ceiling is +.223
  and +.300, and the defect arm's AUROC on the missing stratum is .502 and .504, exactly chance.
  RL0205 therefore is not evidence against the device.
- fe96ab3c's strong-region AP improvement does not replicate under a LightGBM learner (.568 against
  .540 on umap, .601 against .600 on scaffold, both favouring the ceiling), so it should not be
  carried forward as a general consequence of masking weak labels.
- The published 7.06% prevalence, the observed 21.64% official-label prevalence and the measured
  weak-primary rate are mutually compatible with one mixture and identify its strong-primary share
  at `phi = .308` to `.326` (.21 to .42 across 5% to 9% priors).
- Under the official-labels-only reweighted metric, `pool_lgbm` does not beat the official-only
  control (−.021 [−.062, +.026] umap, +.022 [−.021, +.066] scaffold) and `legacy_sep15_refit`
  remains the best recipe. Under the composition-matched metric at the identified `phi` it beats
  it by +.064 [+.022, +.110] and +.114 [+.070, +.163], and the crossover lies at `phi = .85` or
  beyond, a composition the published prevalence excludes.
- With the grid extended to 50% the rate optimum is interior in 14 of 16 cells and lies at 94 to
  112 blind calls for the deployed recipe, identically under both metrics. The prior does not move
  a single argmax. The curve from 49 calls to the optimum spans only +.009 to +.019.

Not established:

- That the measured rate transfers to D6. Three enzymes validate it for themselves; D6's positive
  rate does not collapse toward the cut the way theirs does (26.7% at `>= -0.75` against a pooled
  1.1%), and the extrapolation crosses a 0.33-unit gap in which D6 has no labels at all.
- That the device helps D6. Part 2's evidence is from 3A4, with a dose of 0.39 against D6's 1.80,
  and from 2C9 and 1A2 at doses of 1.23 and 1.54 where the measurement has no power.
- Which blind rows belong to which stratum. Metric M identifies the *composition* of the scored
  set, not the membership. The blind rows were never primary screened.
- Any Live effect size for a recipe change. Part 3's advantage is an offline quantity on a
  reweighted development population; the Live resolution remains about .04 per scored true positive
  against 26 scored positives.
- That the 1.99 calibration factor between the offline mixture precision and the realised Live
  precision is real. It is one observation on one recipe and one band.

## Limitations of the procedure itself

- Part 1's bound is a binomial bound on an observed stratum. It does not account for the labelling
  selection rule itself: 58 D6 primary hits and 2,207 3A4 primary hits were never labelled, so the
  labelled weak-primary rows are a non-random subset of weak-primary molecules.
- Part 2 uses one learner for every arm, so arm differences are learner-specific to the extent that
  LightGBM's response to added negatives is. The reconciliation with fe96ab3c shows that one
  component of that earlier result was in fact learner-specific.
- Part 2's composition-matched metric reweights a 130-row stratum to two thirds of the mass for 1A2,
  giving an effective sample size of 285 against a nominal 1,291. 2C9 is 346 and 3A4 is 732.
- `fixed_full` trains on unlabelled weak-primary rows whose labels are bounded by Part 1 rather than
  observed, so its training signal inherits Part 1's assumption. Its scoring does not: all judged
  rows carry real labels.
- The umap family has 50 assignment groups, so its group bootstrap is coarse; scaffold has 5,367.
  Where the two families disagree, umap's interval is the less trustworthy.
- Metric M assumes each blind stratum is sampled by its development counterpart. Analog expansion
  around strong 1A2/2C9/3A4 hits makes that an approximation of unknown quality.
- Part 4's prediction treats the added ranks as scored at the overall scored fraction of .499,
  whereas the fitted 49-call payloads are scored at .57 to .65 and the 175-call payload at .509.
  If ranks 50–94 are enriched in scored rows the predicted deltas scale up in both directions.
- The Live state underlying Part 4 is a reconstruction from published bootstrap aggregates, not
  disclosed labels.
