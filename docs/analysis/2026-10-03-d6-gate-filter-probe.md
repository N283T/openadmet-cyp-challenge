# Diagnosing the D6 activity gate and testing a gate-filtered vector

Three recorded classification runs. None fitted a model, created a candidate, prepared a
submission or contacted the challenge API.

| Part | Run | Seconds | Code |
| --- | --- | ---: | --- |
| 1 Gate diagnosis | `20261003T041642823714Z-1ab63580` | 0.2 | `scripts/evaluate_cyp_d6_gate_diagnosis.py`, `scripts/launch_cyp_d6_gate_diagnosis.py`, `configs/classification/d6-gate-diagnosis-v1.json` |
| 2 Filtered vector | `20261003T042011147407Z-15430428` | 3.9 | `scripts/evaluate_cyp_d6_gate_filter_probe.py`, `scripts/launch_cyp_d6_gate_filter_probe.py`, `configs/classification/d6-gate-filter-probe-v1.json` |
| 3 Option comparison | `20261003T042954902078Z-fa7d7e77` | 0.5 | `scripts/evaluate_cyp_d6_option_comparison.py`, `scripts/launch_cyp_d6_option_comparison.py`, `configs/classification/d6-option-comparison-v1.json` |

Outputs with SHA-256 provenance are promoted to
`artifacts/tracked/cyp-d6-gate-diagnosis-20261003T041642823714Z-1ab63580/`,
`artifacts/tracked/cyp-d6-gate-filter-probe-20261003T042011147407Z-15430428/` and
`artifacts/tracked/cyp-d6-option-comparison-20261003T042954902078Z-fa7d7e77/`.
Retention script: `scripts/retain_cyp_d6_gate_filter_probe.py`.
Two executions are not promoted. `20261003T041633081427Z-648209ce` failed immediately on a
misplaced array reduction in the circularity table, produced no scientific output, and was
removed with a recorded discard note. `20261003T042325736922Z-c0341b5d` is the first execution
of Part 3, run before its source was reformatted; its outputs are byte-identical to the retained
run and it is preserved as a run record with a recorded reason.

All quantities below are columns in organiser-supplied training files, frozen score tables,
or reconstructions from published leaderboard aggregates. No biological or clinical claim is
made or needed.

## The stated contradiction, and why it is not one

Every supplied starting fact reproduces exactly. The deployed D6 payload
`20261002T223229553466Z-a08301d2` calls 49 of 750 blind rows positive. Measured against
`gate_x_tdi.p_active` of run `20261002T215345688167Z-7920b49a`, 40 of those 49 sit below
0.5, 21 below 0.3 and 7 below 0.2; among the 701 uncalled rows 128 are above 0.5 and 251
above 0.3. The reconstructed equal percentile-rank mean of the deployed September 15 blind
score and `pool_lgbm` reproduces the 49 deployed calls exactly, and the September 15 score
alone reproduces the 49 September 15 calls exactly.

The claimed tension is that RL0219 found 0 positives among 638 labelled rows at primary
`log2fc >= -0.25`, so the 40 low-gate calls ought to contribute almost nothing, while
RL0218 attributes 11 scored true positives to the payload out of about 28 scored calls.

**The two statements are about different quantities and are both true.** `p_active` is a
calibrated estimate of `P(primary log2fc < -1)`. The cross-enzyme rule concerns
`log2fc >= -0.25`. The band between them, `[-1, -0.25)`, is neither: 297 of the 1,497
officially labelled D6 rows have `log2fc >= -1` and 22.2% of them are positive, against an
official base rate of 21.6%. A low `p_active` therefore does not place a row in the
relabelled region; it places it somewhere in a distribution that is still mostly above the
gate threshold but mostly below the relabelling threshold. Part 1 measures each link
separately and finds that the premise joining them fails.

## Part 1 — diagnosing the gate

### The blind p_active distribution is shifted up, and it looks like the unscreened stratum

Held-fold `gate_x_tdi.p_active` on all 6,146 development molecules against the 750 blind
molecules (`p_active_distributions.csv`, `p_active_shift.csv`):

| Population | n | mean | q10 | q50 | q90 | >= 0.2 | >= 0.3 | >= 0.5 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **blind** | **750** | **.3094** | .0916 | .2140 | .7205 | .533 | .372 | .183 |
| umap development, all | 6,146 | .2812 | .0820 | .2103 | .5960 | .526 | .353 | .158 |
| umap screened hit (`log2fc < -0.5`) | 1,553 | .4172 | .1308 | .3729 | .7875 | .779 | .609 | .348 |
| umap screened non-hit (`log2fc >= -0.5`) | 2,823 | .1958 | .0704 | .1483 | .3937 | .343 | .182 | .049 |
| umap unscreened (no `log2fc`) | 1,770 | .2982 | .0952 | .2427 | .5828 | .595 | .401 | .164 |
| scaffold development, all | 6,146 | .2778 | .0702 | .1976 | .6431 | .496 | .335 | .165 |
| scaffold screened hit | 1,553 | .4401 | .1236 | .3895 | .8511 | .772 | .619 | .375 |
| scaffold screened non-hit | 2,823 | .1826 | .0605 | .1314 | .3756 | .304 | .158 | .050 |
| scaffold unscreened | 1,770 | .2873 | .0827 | .2265 | .6069 | .560 | .368 | .166 |

The direction is up. The blind mean exceeds the whole-development mean by +.028 (umap) and
+.032 (scaffold), KS .056 and .095. The shift is small against the stratum differences: the
screened-hit and screened-non-hit means differ by .22 to .26, roughly eight times the
blind-versus-development gap.

The blind distribution is closest to the unscreened development stratum, which is the right
comparator because those rows also carry no primary measurement: mean difference +.011
(umap) and +.022 (scaffold), KS .076 and .059, and the AUROC for ranking a blind row above
an unscreened development row is .481 and .510, that is, indistinguishable. Against the
screened hits the blind set is clearly lower (AUROC .350 and .343) and against the screened
non-hits clearly higher (.646 and .682). **The blind rows are not out of the gate's observed
range; they occupy the part of it that the unscreened development rows occupy.**

### On development held folds the gate is well calibrated

Binned predicted `p_active` against the observed rate of `log2fc < -1` on the 4,376 screened
development molecules, held-fold scores, Wilson intervals (`gate_reliability.csv`):

| Bin | umap n | predicted | observed | reliability | scaffold n | predicted | observed | reliability |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.0–0.1 | 866 | .0770 | .0520 | **−.0250** | 1,081 | .0714 | .0444 | **−.0270** |
| 0.1–0.2 | 1,333 | .1447 | .1283 | −.0165 | 1,239 | .1435 | .1340 | −.0095 |
| 0.2–0.3 | 718 | .2448 | .2660 | +.0212 | 649 | .2435 | .2666 | +.0231 |
| 0.3–0.4 | 468 | .3465 | .3590 | +.0125 | 416 | .3477 | .3822 | +.0345 |
| 0.4–0.5 | 312 | .4485 | .4712 | +.0226 | 269 | .4496 | .5279 | **+.0783** |
| 0.5–0.6 | 236 | .5490 | .5805 | +.0315 | 197 | .5484 | .5381 | −.0103 |
| 0.6–0.7 | 166 | .6477 | .6807 | +.0330 | 148 | .6526 | .6351 | −.0174 |
| 0.7–0.8 | 127 | .7474 | .7480 | +.0006 | 138 | .7501 | .7681 | +.0180 |
| 0.8–0.9 | 100 | .8511 | .8600 | +.0089 | 130 | .8484 | .8154 | −.0330 |
| 0.9–1.0 | 50 | .9235 | .9000 | −.0235 | 109 | .9342 | .8991 | −.0351 |

| Split | screened | base rate | mean p_active | ECE | max deviation | Brier | Brier of base rate | AUROC | bins containing the prediction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| umap | 4,376 | .2738 | .2744 | **.0198** | .0330 | .1487 | .1988 | .8045 | 9 of 10 |
| scaffold | 4,376 | .2738 | .2740 | **.0243** | .0783 | .1418 | .1988 | .8252 | 8 of 10 |

The mean prediction matches the base rate to four decimals on both families, the expected
calibration error is 2.0 and 2.4 percentage points, and the only systematic deviation is a
mild over-prediction in the lowest bin (predicted .071–.077 against observed .044–.052).
The gate is a competent, well-calibrated predictor of what it was trained to predict.

### The unscreened development rows are not systematically different

| Split | unscreened | mean | screened | mean | difference | AUROC | Mann-Whitney p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| umap | 1,770 | .2982 | 4,376 | .2744 | +.0238 | .5594 | 1.5e-07 |
| scaffold | 1,770 | .2873 | 4,376 | .2740 | +.0133 | .5510 | 4.6e-05 |

The unscreened rows sit .013 to .024 above the screened ones. The difference is significant
at this sample size and negligible in magnitude: an AUROC of .55 against .50. The 1,770
figure is the 1,768 `unscreened` rows plus the two labelled rows with no primary measurement
already noted in RL0219.

### The decisive measurement: p_active carries no positive information about the endpoint

Officially labelled development rows binned by held-fold `p_active`, with the observed
`CYP2D6_is_TDI` positive rate and Wilson intervals (`tdi_rate_by_p_active.csv`):

| Bin | umap n | positives | rate | 95% CI | scaffold n | positives | rate | 95% CI |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |
| 0.0–0.2 | 325 | 75 | **.2308** | [.188, .280] | 335 | 77 | **.2299** | [.188, .278] |
| 0.2–0.3 | 255 | 70 | .2745 | [.223, .332] | 227 | 60 | .2643 | [.211, .325] |
| 0.3–0.5 | 388 | 82 | .2113 | [.174, .255] | 366 | 84 | .2295 | [.189, .275] |
| 0.5–1.0 | 529 | 97 | **.1834** | [.153, .219] | 569 | 103 | **.1810** | [.152, .215] |

| Split | official rows | positives | base rate | AUROC of p_active for TDI | Spearman | p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| umap | 1,497 | 324 | .2164 | **.4514** | **−.0694** | .0072 |
| scaffold | 1,497 | 324 | .2164 | **.4608** | **−.0559** | .0304 |

**Among the rows that carry a TDI label, `p_active` is significantly negatively associated
with being positive.** The AUROC is below one half on both families and the Spearman
correlation is negative with p = .007 and p = .030. The lowest gate bin is 23.0–23.1%
positive and the highest is 18.1–18.3%. This is the measurement that decides the proposal:
a filter that removes low-gate calls removes rows drawn from a 23% stratum and admits rows
drawn from an 18% stratum, before any loss of score rank is counted.

### Composing the gate with the measured region rates

Among screened development rows, the share falling in each primary region by `p_active`
band, and the positive rate that composition implies. Region rates are `.2150` for
`log2fc < -1` and `.2222` for `[-1, -0.5)` from the D6 labels themselves, `.0091` for
`[-0.5, -0.25)` from the pooled 1A2/2C9/3A4 measurement (5 of 549), and three scenarios for
the relabelled region `>= -0.25`, where D6 has no labels at all
(`region_decomposition.csv`, `implied_blind_rate.csv`):

| Split | Band | n | share `< -1` | share `[-1,-0.5)` | share `[-0.5,-0.25)` | share `>= -0.25` | implied rate, clean | at the .0047 bound | at a .0706 relabelled rate |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| umap | 0.0–0.2 | 2,199 | .098 | .058 | .023 | .820 | .0343 | .0381 | .0922 |
| umap | 0.2–0.3 | 718 | .266 | .102 | .038 | .595 | .0801 | .0829 | .1221 |
| umap | 0.3–0.5 | 780 | .404 | .114 | .040 | .442 | .1125 | .1146 | .1438 |
| umap | 0.5–1.0 | 679 | .701 | .096 | .031 | .172 | .1723 | .1731 | .1844 |
| scaffold | 0.0–0.2 | 2,320 | .092 | .060 | .023 | .825 | .0334 | .0373 | .0917 |
| scaffold | 0.2–0.3 | 649 | .267 | .100 | .037 | .596 | .0799 | .0827 | .1220 |
| scaffold | 0.3–0.5 | 685 | .439 | .114 | .048 | .399 | .1202 | .1221 | .1484 |
| scaffold | 0.5–1.0 | 722 | .706 | .100 | .028 | .166 | .1743 | .1751 | .1860 |

Applying those shares to the blind `p_active` distribution:

| Split | Scenario | implied blind prevalence | ratio to the published .0706 | implied rate, blind `p_active < 0.2` | 0.2–0.5 | `>= 0.5` |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| umap | relabelled region clean | .0817 | 1.157 | .0343 | .0976 | .1723 |
| umap | relabelled at the .0047 bound | .0845 | 1.197 | .0381 | .1000 | .1731 |
| umap | relabelled at the blind base rate | .1236 | 1.751 | .0922 | .1338 | .1844 |
| scaffold | relabelled region clean | .0831 | 1.177 | .0334 | .1017 | .1743 |
| scaffold | relabelled at the .0047 bound | .0859 | 1.216 | .0373 | .1040 | .1751 |
| scaffold | relabelled at the blind base rate | .1245 | 1.764 | .0917 | .1362 | .1860 |

Two things follow.

First, **a blind row below `p_active` 0.2 has an implied positive rate of .033 to .038, about
half the blind base rate, not nothing.** Only 82% of such rows fall in the relabelled region;
the remaining 18% are primary hits, where the positive rate is 21–22%. Seven of the deployed
49 calls sit below 0.2 and 21 below 0.3; at an implied rate near .03 to .08 those rows are
expected to contribute roughly one scored positive between them, which is entirely consistent
with 11 scored true positives arising mostly from the 28 calls above 0.3.

Second, **the "D6 is genuinely different" horn is the one the prevalence identity disfavours.**
With the relabelled region clean, the decomposition implies a blind prevalence of .082 to .083
against the published .0706, a 16–18% overshoot that the difference between full-set and
scored-subset prevalence can easily absorb. Pushing the relabelled region to the blind base
rate implies .124, a 75% overshoot that would have to be cancelled by the scored mask being
strongly depleted in positives, which the identity returning .0706 from all eight states makes
unattractive.

### A pseudo-negative-labelled metric cannot adjudicate a gate filter

| Split | official | pseudo-negative | AUROC of p_active separating them | mean official | mean pseudo-negative | official share above 0.3 | above 0.5 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| umap | 1,497 | 2,693 | **.8047** | .4194 | .1902 | .665 | .819 |
| scaffold | 1,497 | 2,693 | **.8233** | .4431 | .1768 | .704 | .826 |

`p_active` separates the official stratum from the pseudo-negative stratum with AUROC .80 to
.82, and 82% of the development rows above `p_active` 0.5 are official rows. Pseudo-negative
labels are defined by primary `log2fc` and `p_active` is a calibrated predictor of primary
`log2fc`, so raising a gate threshold mechanically concentrates calls on the stratum whose
labels are real and away from the stratum that is labelled negative by construction. **Any
metric that scores pseudo-negative rows rewards a gate filter partly by construction.** Metric
M is reported in Part 2 exactly as specified, and read with that discount.

### Part 1 verdict

**Neither horn as stated is supported, because the premise that joins them is false.** The
gate is well calibrated on development held folds (ECE .020 and .024, mean prediction equal to
the base rate to four decimals, prediction inside the bin interval in 9 and 8 of 10 bins). Its
blind distribution is shifted up by .028 to .032 and is statistically indistinguishable from
the unscreened development stratum, so blind extrapolation is interpolation within the gate's
observed range rather than a leap. And a low `p_active` does not imply membership of the clean
relabelled region: 18% of low-gate rows are primary hits, giving an implied positive rate of
.033 to .038 rather than zero. The direct measurement closes the question: among labelled
rows, `p_active` is *negatively* associated with TDI positivity (AUROC .451 and .461, Spearman
−.069 and −.056, p = .007 and .030). There is nothing to reconcile, and no need to invoke
either a miscalibrated gate or a D6-specific weak-region excess. Of the two offered horns, the
second is the one the prevalence identity actively disfavours (implied blind prevalence .124
against .0706).

## Part 2 — the filtered vector evaluated offline

### Design

The control is `l4_control`, the out-of-fold analogue of the deployed recipe: the equal
per-fold percentile-rank mean of `legacy_sep15_refit` and `pool_lgbm`, rebuilt from the frozen
held-fold scores of run `7920b49a`. Four filtered variants restrict eligibility to
`p_active >= tau` for tau in {0.2, 0.3, 0.4, 0.5} and fill any shortfall from the excluded rows
by control rank, implemented as `control_rank + 10 * (p_active >= tau)`, which is the same
lexicographic ordering used for the L5 variants of run `2d8d20f5`. The two multiplier forms
already tested against the legacy control are included as `l4_x_p_active` and
`l4_x_sqrt_p_active`. Nothing is refitted.

Metric O scores officially labelled rows only with the classes reweighted to a target scored
prevalence; Metric M is the composition-matched two-stratum mixture. Both are the estimators of
run `008f8427`, with the same rate grid extended to 50% and the same paired group bootstrap
(1,000 resamples of assignment groups). Priors .05/.0706/.09 and mixture fractions
.23/.3262/.42 are all evaluated; the primary cells are prior .0706 and phi .3262.

The verdict rule was written into `configs/classification/d6-gate-filter-probe-v1.json` before
any Part 2 number was read: candidate-worthy requires a positive Metric M delta with a 95%
interval excluding zero on both split families, no significantly negative Metric O delta on
either family, and no contradiction from the historical transitions.

### Paired group bootstrap against the unfiltered control at 49/750

| Variant | Metric M phi .3262, umap | Metric M phi .3262, scaffold | Metric O prior .0706, umap | Metric O prior .0706, scaffold |
| --- | --- | --- | --- | --- |
| `l4_filter_tau0.2` | +.0201 [−.0054, +.0478] | +.0084 [−.0129, +.0317] | −.0090 [−.0310, +.0149] | −.0157 [−.0367, +.0059] |
| `l4_filter_tau0.3` | −.0034 [−.0378, +.0293] | +.0269 [−.0101, +.0633] | **−.0365 [−.0728, −.0031]** | −.0220 [−.0566, +.0105] |
| `l4_filter_tau0.4` | −.0101 [−.0500, +.0347] | +.0228 [−.0242, +.0692] | **−.0525 [−.0984, −.0098]** | −.0374 [−.0816, +.0029] |
| `l4_filter_tau0.5` | −.0093 [−.0571, +.0366] | −.0118 [−.0652, +.0415] | **−.0712 [−.1298, −.0249]** | **−.0849 [−.1323, −.0415]** |
| `l4_x_p_active` | −.0080 [−.0584, +.0501] | +.0066 [−.0479, +.0644] | **−.0776 [−.1308, −.0299]** | **−.0818 [−.1299, −.0364]** |
| `l4_x_sqrt_p_active` | −.0223 [−.0636, +.0232] | +.0078 [−.0416, +.0548] | **−.0787 [−.1242, −.0379]** | **−.0694 [−.1148, −.0275]** |

Metric M sensitivity to the mixture fraction, at 49/750 (`bootstrap_deltas.csv`):

| Variant | phi .23, umap | phi .42, umap | phi .23, scaffold | phi .42, scaffold |
| --- | --- | --- | --- | --- |
| `l4_filter_tau0.2` | +.0231 [−.0002, +.0484] | +.0163 [−.0093, +.0460] | +.0118 [−.0083, +.0336] | +.0046 [−.0180, +.0291] |
| `l4_filter_tau0.3` | +.0034 [−.0285, +.0321] | −.0104 [−.0491, +.0254] | +.0337 [+.0002, +.0676] | +.0194 [−.0204, +.0577] |
| `l4_filter_tau0.4` | −.0001 [−.0364, +.0399] | −.0200 [−.0633, +.0291] | +.0331 [−.0109, +.0760] | +.0122 [−.0368, +.0611] |
| `l4_filter_tau0.5` | +.0061 [−.0359, +.0491] | −.0240 [−.0786, +.0269] | +.0070 [−.0429, +.0553] | −.0291 [−.0846, +.0280] |

Lowering phi towards .23 moves every delta up and raising it towards .42 moves every delta
down, which is exactly what the Part 1 circularity measurement predicts: the more mass the
pseudo-negative stratum carries, the more a gate filter is rewarded for separating strata
rather than for ranking the endpoint. Even at phi .23 no variant has a positive interval
excluding zero on both families.

### Rate optima and the value at a fixed 49 calls

| Split | Metric | Variant | best rate | best calls of 750 | best MCC | MCC at 49 calls |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| umap | O .0706 | `l4_control` | 15% | 112 | .0816 | **.0678** |
| umap | O .0706 | `l4_filter_tau0.2` | 2% | 15 | .0783 | .0588 |
| umap | O .0706 | `l4_filter_tau0.3` | 2% | 15 | .0426 | .0313 |
| umap | O .0706 | `l4_filter_tau0.4` | 1% | 8 | .0222 | .0153 |
| umap | O .0706 | `l4_filter_tau0.5` | 50% | 375 | .0381 | −.0034 |
| scaffold | O .0706 | `l4_control` | 27.5% | 206 | .1158 | **.1067** |
| scaffold | O .0706 | `l4_filter_tau0.2` | 2% | 15 | .1081 | .0910 |
| scaffold | O .0706 | `l4_filter_tau0.3` | 4% | 30 | .0915 | .0846 |
| scaffold | O .0706 | `l4_filter_tau0.4` | 1% | 8 | .0967 | .0693 |
| scaffold | O .0706 | `l4_filter_tau0.5` | 1% | 8 | .0869 | .0218 |
| umap | M .3262 | `l4_control` | 15% | 112 | .1199 | .1008 |
| umap | M .3262 | `l4_filter_tau0.2` | 50% | 375 | .1678 | .1210 |
| umap | M .3262 | `l4_filter_tau0.3` | 40% | 300 | .1552 | .0974 |
| scaffold | M .3262 | `l4_control` | 27.5% | 206 | .1476 | .1384 |
| scaffold | M .3262 | `l4_filter_tau0.2` | 35% | 262 | .1916 | .1468 |
| scaffold | M .3262 | `l4_filter_tau0.3` | 25% | 188 | .1807 | .1653 |

Under Metric O the unfiltered control has the highest MCC at 49 calls on both families and
every variant is below it, the shortfall growing monotonically with the threshold from −.009
to −.071 on umap and from −.016 to −.085 on scaffold. The filters' own Metric O optima
collapse to 8–30 calls, which is the curve telling us that once ineligible rows are removed
the remaining ordering is worth very little past the first few dozen rows. Under Metric M the
filters push the selected rate to 25–50% and gain there, which is the circularity operating
at full strength and not a deployable statement about 49 calls.

### Nested rate selection

Rates chosen on the other four folds under the same metric (`nested_selection.csv`). Under
Metric M the filters select 22.5% to 50% of the pool and beat the control by +.05 to +.08 with
intervals excluding zero; under Metric O the picture reverses, with `l4_filter_tau0.4` at
−.0655 [−.1043, −.0266] on umap. Those nested figures move the operating point by a factor of
four to eight and therefore do not speak to the fixed-count proposal; they are reported because
the design asked for them.

### Historical transition concordance

The frozen design of run `2d8d20f5` applied to the eight scored transitions, with each
variant's blind ordering converted to a percentile rank over the 750 rows
(`transitions.csv`, `concordance_agreement.csv`):

| Column | sign agreements | transitions with permutation p < .05 | agrees on sep15 to oct03 | contains deployed components |
| --- | ---: | ---: | --- | --- |
| `l4_control` | **8 of 8** | 6 | yes | yes |
| `l4_filter_tau0.4` | 8 of 8 | 6 | yes | yes |
| `l4_filter_tau0.5` | 8 of 8 | 5 | yes | yes |
| `l4_filter_tau0.2` | 7 of 8 | 5 | yes | yes |
| `l4_filter_tau0.3` | 7 of 8 | 4 | yes | yes |
| `l4_x_p_active` | 7 of 8 | 4 | yes | yes |
| `l4_x_sqrt_p_active` | 6 of 8 | 5 | yes | yes |

**This test cannot distinguish the proposal from the control, because the control itself scores
8 of 8.** Every column here contains the deployed September 15 blind score and `pool_lgbm`,
both components of the score that selected six of the eight payloads and the eighth in
particular, so the agreement is near-tautological in exactly the sense RL0218 flagged for
`pool_lgbm`. The three columns RL0218 recorded at 7 of 8 are `gate_x_tdi.p_active`,
`gate_x_tdi_chemeleon.p_active` and `pool_lgbm.score`; the filtered variants here reach 6 to 8
of 8 but none of them beats the unfiltered control, and the two that tie it at 8 of 8
(tau = 0.4 and 0.5) are the two with the worst Metric O deltas. The honest reading is that the
concordance criterion is saturated for anything built from the deployed score and supplies no
evidence either way.

### Blind swap inventory at 49 calls

Against the deployed payload, with the similarity-weighted top-10 labelled positive rate reused
from run `20261002T223515148996Z-ccb05cf3` (`blind_swaps.csv`):

| Variant | overlap | added | removed | mean p_active added | mean p_active removed | neighbour rate retained | added | removed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `l4_control` | 49 | 0 | 0 | – | – | .4407 | – | – |
| `l4_filter_tau0.2` | 42 | 7 | 7 | .4823 | .1304 | .4360 | **.3536** | .4693 |
| `l4_filter_tau0.3` | 28 | 21 | 21 | .5265 | .2122 | .4398 | **.2785** | .4420 |
| `l4_filter_tau0.4` | 20 | 29 | 29 | .5651 | .2516 | .4588 | **.2627** | .4282 |
| `l4_filter_tau0.5` | 9 | 40 | 40 | .6520 | .3046 | .4855 | **.1962** | .4306 |
| `l4_x_p_active` | 19 | 30 | 30 | .6874 | .2565 | .4695 | **.2153** | .4225 |
| `l4_x_sqrt_p_active` | 24 | 24 | 24 | .6128 | .2272 | .4522 | **.2608** | .4287 |

The mechanism works as intended on the gate and against the endpoint. Every variant raises the
mean `p_active` of its calls, and every variant swaps in rows whose similarity-weighted labelled
neighbourhood is much poorer in positives than the rows it drops: .354 against .469 at
tau = 0.2, and .196 against .431 at tau = 0.5, where the retained calls sit at .44 to .49. This
is an independent confirmation of the Part 1 finding, computed on labelled-neighbour support
rather than on the gate itself. The deployed 49 calls have mean `p_active` .3682 with a minimum
of .1070, so the deployed vector is already enriched in gate probability relative to the blind
mean of .3094, which is why tau = 0.2 touches only 7 calls.

### Verdict against the declared rule

| Variant | Metric M positive with interval excluding zero on both splits | min Metric M delta | Metric O not significantly negative on either split | min Metric O delta | candidate-worthy |
| --- | --- | ---: | --- | ---: | --- |
| `l4_filter_tau0.2` | no | +.0084 | yes | −.0157 | **no** |
| `l4_filter_tau0.3` | no | −.0034 | no | −.0365 | **no** |
| `l4_filter_tau0.4` | no | −.0101 | no | −.0525 | **no** |
| `l4_filter_tau0.5` | no | −.0118 | no | −.0849 | **no** |
| `l4_x_p_active` | no | −.0080 | no | −.0818 | **no** |
| `l4_x_sqrt_p_active` | no | −.0223 | no | −.0787 | **no** |

**Every variant fails, and most fail both conditions.** The best threshold, tau = 0.2, is the
only one that avoids a significantly negative Metric O delta, and its Metric M deltas are
+.0201 [−.0054, +.0478] and +.0084 [−.0129, +.0317], both intervals spanning zero. Thresholds
0.3 and above are significantly negative under Metric O on at least one family; both multiplier
forms are significantly negative under Metric O on both. No candidate is built.

## Part 3 — the three options in public-metric units

The sizing basis is the one established in RL0218 and used in RL0219: 374 scored rows of 750,
26 scored positives, the deployed payload fitted at 28 scored calls with 11 true positives, 17
false positives and 15 false negatives, point MCC .36165 against a reported Live .364248. Each
option is simulated exactly rather than linearised: the true positives among the scored removed
calls are drawn from the deployed realised precision 11/28 = .3929, those among the scored added
calls from each marginal-precision scenario, and MCC is recomputed from the resulting integer
confusion matrix, 200,000 draws per cell (`option_comparison.csv`). The quoted interval is the
10th to 90th percentile of outcomes at a fixed assumed precision; it is a spread over outcomes,
not a posterior over the precision.

| Option | calls | changed | scored changed | marginal-precision scenario | expected delta | sd | 80% interval |
| --- | ---: | ---: | ---: | --- | ---: | ---: | --- |
| gate filter tau 0.2 | 49 | 7 | 4.0 | replacements as precise as the calls they displace (.393) | +.0000 | .055 | [−.080, +.080] |
| gate filter tau 0.2 | 49 | 7 | 4.0 | offline band ratio (.265) | **−.0204** | .053 | [−.080, +.040] |
| gate filter tau 0.2 | 49 | 7 | 4.0 | labelled-neighbour ratio (.296) | −.0154 | .053 | [−.080, +.040] |
| gate filter tau 0.2 | 49 | 7 | 4.0 | replacements at the base rate (.071) | −.0515 | .044 | [−.120, +.000] |
| gate filter tau 0.3 | 49 | 21 | 12.0 | offline band ratio (.265) | −.0612 | .091 | [−.160, +.040] |
| gate filter tau 0.3 | 49 | 21 | 12.0 | labelled-neighbour ratio (.248) | −.0699 | .090 | [−.200, +.040] |
| gate filter tau 0.5 | 49 | 40 | 22.9 | labelled-neighbour ratio (.179) | −.1859 | .107 | [−.320, −.040] |
| extension to 94 | 94 | 45 | 22.5 | uncalibrated mixture (.133) | −.0389 | .049 | [−.098, +.025] |
| extension to 94 | 94 | 45 | 22.5 | break-even (.188) | −.0015 | .057 | [−.067, +.087] |
| extension to 94 | 94 | 45 | 22.5 | calibrated by 1.99 (.265) | +.0509 | .064 | [−.037, +.149] |
| **hold** | 49 | 0 | 0 | – | **.0000** | .000 | [.000, .000] |

Sensitivity to the scored fraction of the changed ranks, which the reconstruction cannot pin
down (49-call payloads are fitted at .57 to .65 scored against an overall .499):

| Option | scored fraction .499 | .571 | .650 |
| --- | ---: | ---: | ---: |
| gate filter tau 0.2, offline band ratio | −.0152 | −.0203 | −.0255 |
| extension to 94, uncalibrated mixture | −.0388 | −.0413 | −.0430 |

How often the pre-registered extension produces an outcome its own decision rule calls
decisive, that is, a delta outside [−.080, +.080] (`extension_decisiveness.csv`):

| Scenario | marginal precision | P(delta > +.080) | P(delta < −.080) | P(decisive) | P(inside, no follow-up permitted) |
| --- | ---: | ---: | ---: | ---: | ---: |
| uncalibrated mixture | .133 | .020 | .190 | **.210** | .790 |
| break-even | .188 | .102 | .062 | **.164** | .836 |
| calibrated by 1.99 | .265 | .359 | .010 | **.370** | .630 |

### Recommendation: hold the current vector

**Hold.** Of the three options it is the only one with a non-negative expected change, and it
is the only one that does not consume a 12-hour opportunity.

The gate filter is ruled out by Part 2 on its own declared rule and by Part 3 on expected value:
at the one threshold that is not significantly negative under the board's own label set it
changes 7 calls for an expected −.015 to −.020 with an 80% interval of [−.080, +.040], and the
only scenario in which it breaks even is the one Part 1 directly measured to be false, namely
that the replacements are as good as the calls they displace. Every larger threshold is worse
in proportion to how many calls it moves.

### On the 94-row extension, I agree with the stated scepticism, and for a further reason

The three grounds together:

1. **The favourable arm rests on one observation extrapolated in the wrong direction.** The
   1.99 factor is the ratio of the realised Live precision of ranks 1–49 to the offline mixture
   estimate of the same band. Applying it to ranks 50–94 assumes the bias is multiplicative and
   rank-invariant, in a region where the underlying offline estimate is itself falling
   monotonically (.176, .150, .149, .120, .112 across the umap bands). Nothing in the available
   data supports rank-invariance, and a bias that is partly additive would leave the lower bands
   near the uncalibrated −.039.
2. **The action the probe licenses is worth less than the measurement noise.** RL0219's own rate
   curve puts the Metric M gain from 49 calls to the offline optimum at +.019 (umap) and +.009
   (scaffold), against a per-scored-true-positive quantum of .0399 and a Null B spread of .0435.
   A probe whose purpose is to decide whether to move to an operating point worth less than one
   fifth of the metric's resolution cannot pay for itself.
3. **The probe is unlikely to be decisive, by its own pre-registered rule.** The simulation above
   puts the probability of landing outside [−.080, +.080] at .16 to .37 depending on which
   bracket is right. In 63% to 84% of outcomes the result is "the marginal precision is between
   .07 and .31", which is approximately the prior, and the pre-registration then forbids a second
   upload on the same question. Expected information per opportunity is low, and the identifying
   assumption is fragile: a scored fraction of .571 instead of .499 rescales every predicted
   delta by about 1.14, which alone moves outcomes across the decision boundary.

A fourth point is worth recording because it cuts the other way and still does not rescue the
probe. The board displays the latest submission rather than the best, and the project has
already restored a superseded best once, so a loss costs one opportunity plus a re-send rather
than the record. That makes the probe cheap, not valuable. With roughly sixty opportunities left
before 2026-11-03, the binding constraint is which questions are worth asking.

### What the next offline step should be instead

1. **Close the gate as a call filter.** Part 1 settles it: `p_active` is well calibrated for
   primary activity and negatively associated with the endpoint among rows that carry a label.
   It should not reappear as a filter, a multiplier or a tie-break, and the L5 result of RL0213
   (the one variant with a non-negative official-label interval, on scaffold only) should be read
   as the split-dependent noise that it is rather than as a lead.
2. **Earn the right to act on an offline estimate by testing the estimators where the Live
   measurement has power.** Every D6 methodology question is currently judged against 26 scored
   positives and a .04 quantum. A4 is fitted at about 368 scored rows with 249 calls and a Live
   MCC of .5299, and MA aggregates both endpoints. Replaying the Metric O and Metric M apparatus
   against the A4 and MA historical transitions is purely retrospective, needs no upload, and is
   the only available way to find out whether the composition-matched metric predicts Live
   direction at all. If it does not, no amount of D6 offline work is actionable; if it does, the
   D6 rate question can be decided offline rather than probed.
3. **Measure the seed stability of the deployed vector.** The +.0556 gain is one scored molecule
   and sits at the 76th percentile of the exchangeable null. If the 17 additions that produced it
   are fragile under a reseed or refit of `pool_lgbm`, the current .364248 is partly luck, which
   changes the risk calculus of any future probe and of any re-send. This is a refit of one
   frozen recipe at a handful of seeds, not a new model family.

## Design changes made, and why

1. **Part 1 was extended with the measurement that actually decides the question.** The
   specified diagnosis asked for the blind-versus-development distribution, the gate's
   reliability against `log2fc < -1`, and the unscreened stratum. All three are reported, and all
   three turn out to be favourable to the gate: it is well calibrated and the blind rows are
   inside its observed range. None of them settles the proposal, because the proposal depends on
   whether `p_active` predicts the *endpoint*, which is a different question. The positive rate
   by `p_active` bin among officially labelled rows was therefore added, and it is negative. The
   region decomposition was added for the same reason: to convert `p_active` into an implied
   positive rate rather than leaving it as a probability about a different variable.
2. **The two offered horns were not adopted as the only options.** The brief asked which of a
   miscalibrated gate or a genuinely higher D6 weak-region rate is better supported. Neither is
   needed, because the premise that a low `p_active` places a row in the cross-enzyme clean
   region is false; the decomposition puts 18% of low-gate rows in the primary-hit region. The
   second horn is additionally disfavoured by the prevalence identity. Reporting a forced choice
   between two unsupported explanations would have been worse than reporting that the question
   dissolves.
3. **Metric M is reported as specified and explicitly discounted.** Part 1 measures that
   `p_active` separates the official from the pseudo-negative stratum at AUROC .80 to .82, so a
   gate filter raises Metric M precision partly by construction. The phi sensitivity is reported
   in both directions to show the effect: every delta rises as phi falls towards .23 and falls as
   phi rises towards .42. Metric M remains the primary condition of the declared verdict rule, as
   specified, and the rule fails anyway.
4. **The concordance check is reported and declared uninformative.** The brief asked for the
   filtered variant's tally and for any dependence on components of the deployed vector. The
   dependence is total: all seven columns contain the deployed September 15 score, the unfiltered
   control itself scores 8 of 8, and the two variants that tie it are the two with the worst
   Metric O deltas. The tally is given; no weight is placed on it.
5. **A third run was added for the option comparison.** The brief asked for expected change and
   spread for three options in public-metric units. Doing that by linearising at .0399 per scored
   true positive, as RL0219 did, understates the spread, so each option is simulated from integer
   confusion matrices instead, and the probability that the pre-registered extension is decisive
   under its own rule is computed rather than asserted.
6. **The filter is implemented as a lexicographic ordering rather than as a masked re-rank.**
   `control_rank + 10 * (p_active >= tau)` reproduces "restrict eligibility, take the top k, fill
   any shortfall from the excluded rows by control rank" exactly, matches the L5 implementation of
   run `2d8d20f5` so the two are comparable, and avoids a separate fill path. For the concordance
   test the ordering is converted to a percentile rank over the 750 blind rows so that its mean
   is on the same scale as the other columns.

## What this does and does not establish

Established:

- The gate is well calibrated on development held folds for its own target: expected calibration
  error .0198 (umap) and .0243 (scaffold) over ten bins, mean prediction equal to the .2738 base
  rate to four decimals, Brier .149 and .142 against .199 for the base rate, AUROC .805 and .825,
  and the bin prediction inside the Wilson interval in 9 and 8 of 10 bins. The single systematic
  deviation is over-prediction in the lowest bin by .025 to .027.
- The blind `p_active` distribution is shifted up against development by +.028 and +.032 in the
  mean (KS .056 and .095) and is statistically indistinguishable from the unscreened development
  stratum (mean difference +.011 and +.022, AUROC .481 and .510).
- The 1,770 unscreened development rows have `p_active` .013 to .024 above the screened ones,
  AUROC .551 and .559: significant at this sample size, negligible in size.
- Among the 1,497 officially labelled rows, `p_active` is negatively associated with TDI
  positivity: AUROC .4514 and .4608, Spearman −.0694 (p = .0072) and −.0559 (p = .0304), positive
  rate .231/.230 below `p_active` 0.2 against .183/.181 above 0.5.
- Composing the gate's primary-region shares with the measured region-wise positive rates gives a
  blind row below `p_active` 0.2 an implied positive rate of .033 to .038, and reproduces the
  published .0706 scored prevalence to within a factor 1.16 to 1.22 when the relabelled region is
  taken as clean. Setting the relabelled region to the blind base rate implies .124, a factor
  1.75 overshoot.
- `p_active` separates the official from the pseudo-negative development stratum at AUROC .805
  and .823, with 82% of rows above `p_active` 0.5 being official, so any metric that scores
  pseudo-negative labels rewards a gate filter partly by construction.
- No gate-filtered or gate-multiplied variant of the deployed recipe satisfies the declared
  verdict rule. At 49/750, the best threshold tau = 0.2 is +.0201 [−.0054, +.0478] and
  +.0084 [−.0129, +.0317] under Metric M and −.0090 [−.0310, +.0149] and −.0157 [−.0367, +.0059]
  under Metric O. Thresholds 0.3 and above and both multiplier forms are significantly negative
  under Metric O on at least one family, up to −.0849 [−.1323, −.0415].
- Under Metric O the unfiltered control has the best MCC at 49 calls on both families, and the
  shortfall of the filtered variants grows monotonically with the threshold.
- Every variant swaps in blind rows with a markedly poorer similarity-weighted labelled
  neighbourhood than the rows it drops: .354 against .469 at tau = 0.2, .196 against .431 at
  tau = 0.5, with retained calls at .44 to .49.
- The eight-transition concordance criterion is saturated for any column built from the deployed
  score: the unfiltered control scores 8 of 8 and no filtered variant beats it.
- On the recorded sizing basis the gate filter at tau = 0.2 has an expected delta of −.015 to
  −.020 with an 80% interval of [−.080, +.040]; the pre-registered 94-call extension spans −.039
  to +.051 across its scenarios; holding is zero by construction. The extension produces an
  outcome its own rule calls decisive in only .16 to .37 of draws.

Not established:

- Anything about the gate's calibration on the blind rows themselves. They were never primary
  screened, so no reliability curve exists for them; the distributional match to the unscreened
  development stratum is suggestive and is not a calibration check.
- That the region decomposition transfers. It applies development region shares conditional on
  `p_active` to blind rows, which assumes the gate's conditional distribution of `log2fc` given
  `p_active` is the same in both populations.
- Any Live effect size. Every number in Part 3 is built on confusion matrices reconstructed from
  published bootstrap aggregates, and the scored fraction of the changed ranks is assumed.
- That the negative association between `p_active` and TDI among labelled rows extends to
  unlabelled rows. It cannot, by construction: the labelled rows are all primary hits, and for an
  unconditional population the association must be positive, as the decomposition shows. The two
  facts together are why the deployed score, which already carries mean `p_active` .368 against a
  blind mean of .309, has extracted most of the available gate signal and why pushing further
  trades rank for nothing.
- That no gate-derived quantity could ever help. What is established is that this gate, used as a
  fixed-count eligibility filter or as a multiplier on this score, does not.

## Limitations of the procedure itself

- The officially labelled rows are all primary hits, so the positive-rate-by-`p_active` table
  conditions on a selected population. It is the population whose labels the board scores, which
  is why it is the relevant one, but it is not a random sample of blind-like molecules.
- Region rates for `log2fc < -1` and `[-1, -0.5)` come from D6 labels and those for
  `[-0.5, -0.25)` and above from the pooled other three enzymes, so the decomposition mixes two
  label sources.
- The implied blind prevalence is over all 750 rows while the published .0706 is the scored
  subset's prevalence, and RL0218 found the top-ranked calls modestly enriched in scored rows, so
  the two are not required to agree exactly.
- The umap family has 50 assignment groups against scaffold's 5,367, so its group bootstrap is
  coarse; where the two disagree, umap's interval is the less trustworthy.
- Nested rate selection moves the operating point by a factor of four to eight for the filtered
  variants, so its deltas are not comparable with the fixed-count proposal.
- The 200,000-draw option simulation holds the marginal precision fixed in each cell and varies
  only which molecules the swap touches, so its interval is a sampling spread and not an
  uncertainty over the precision.
- The three scenarios for the replacement precision are a bracket chosen by the analyst from
  recorded quantities, not a fitted distribution.
