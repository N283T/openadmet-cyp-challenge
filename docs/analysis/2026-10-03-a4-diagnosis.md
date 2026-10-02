# CYP3A4 TDI diagnosis: Live lead, population composition and pool/gate variants

Three recorded classification runs, none of which created a candidate, prepared a
submission or uploaded anything:

| Part | Run | Seconds | Code |
| --- | --- | ---: | --- |
| 1 Live versus Interim | `20261002T225609263717Z-1daab817` | 0.6 | `scripts/evaluate_cyp_a4_live_interim.py`, `scripts/launch_cyp_a4_live_interim.py`, `configs/classification/a4-live-interim-v1.json` |
| 2 OOF strata and reweighting | `20261002T231843369550Z-3db76e6a` | 21.2 | `scripts/evaluate_cyp_a4_oof_strata.py`, `scripts/launch_cyp_a4_oof_strata.py`, `configs/classification/a4-oof-strata-v1.json` |
| 3 Pool and gate variants | `20261002T231635959631Z-518d8252` | 38.3 | `scripts/evaluate_cyp_a4_pool_gate.py`, `scripts/launch_cyp_a4_pool_gate.py`, `configs/classification/a4-pool-gate-v1.json` |

Compact outputs with SHA-256 provenance are promoted to
`artifacts/tracked/cyp-a4-live-interim-20261002T225609263717Z-1daab817/`,
`artifacts/tracked/cyp-a4-oof-strata-20261002T231843369550Z-3db76e6a/` and
`artifacts/tracked/cyp-a4-pool-gate-20261002T231635959631Z-518d8252/`.
Earlier attempts are preserved with recorded notes. Two failed outright and produced
no scientific output, so they are recorded as discarded:
`20261002T225951320341Z-84225916` (asserted that the equal-seven OOF covers all
compounds when it covers the 6,146 development molecules only) and
`20261002T230028064290Z-5d56e1ce` (inverted log2fc bin interval). Four more are
recorded as superseded observations: `20261002T225510891345Z-86429ecc` for Part 1,
`20261002T230128781419Z-79a29be6` and `20261002T231501082013Z-0402eb9d` for Part 2,
and `20261002T230753297900Z-33cf990e` with `20261002T231528746692Z-e54b37f5` for
Part 3. Each was superseded by a source change that does not alter the result — an
added split-assignment assertion, a launcher line-length fix, deterministic
concordance row order and a cosmetic reformat — and every decision output of the
superseded Part 2 and Part 3 runs is byte-identical to the retained one. The retained
runs are the ones whose recorded source snapshot matches the checked-in code.

The vector under study is unchanged since September 15: the equal percentile-rank
mean of seven A4 members, top 249 of 750 called positive. Its blind score vector is
`artifacts/runs/classification/20260915T125939829893Z-1d921ffc/outputs/test-scores.parquet`
(column `a4_prediction`). Both Part 2 and Part 3 assert that its top 249 equals the
A4 positives of submission `20261002T092348406577Z-4ea606c1` exactly; the 249th to
250th score gap is `.00019`, so the boundary involves no tie.

## Part 1 — Is the Live A4 lead real?

The frozen September 26 snapshot has 127 Live and 114 Interim classification
entrants; 114 names are common to both boards.

| Endpoint | MCC Spearman | Rank Spearman | Top-10 overlap | Mean Live − Interim | SD of that difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| CYP3A4_is_TDI | .756 | .756 | 3 | −.0074 | .0393 |
| CYP2D6_is_TDI | .915 | .915 | 6 | +.0364 | .0568 |
| MA | .898 | .898 | 6 | +.0145 | .0371 |

A4 is by far the least stable of the three boards between the half-test and
full-test populations. Only three of the Live A4 top 10 appear in the Interim top 10.

### Separation in bootstrap standard deviation units

| Board | N283T | Comparison | Their MCC (std) | Delta | In N283T std | In their std | In pooled std |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| Live A4 | rank 1, .5299 (.0477) | rank 2 harinu123 | .5154 (.0509) | +.0145 | +0.30 | +0.28 | **+0.21** |
| Interim A4 | rank 4, .5023 (.0342) | rank 1 EleanorRigby | .5207 (.0286) | −.0185 | −0.54 | −0.64 | **−0.41** |
| Interim A4 | rank 4, .5023 (.0342) | rank 3 stir_bar | .5035 (.0318) | −.0013 | −0.04 | −0.04 | −0.03 |

On the Live board 16 entrants (N283T included) sit within one pooled bootstrap std
of N283T and 68 within two. On the Interim board the figures are 20 and 60. The
lead over rank 2 is one fifth of a pooled standard deviation; the deficit to the
Interim leader is four tenths of one.

### Inferred operating points

Per-entrant scored positive counts come from the recall bootstrap spread,
`P ≈ r(1 − r)/std_r²`, and call fractions from
`prevalence ≈ (1 − accuracy)/(1 − 2r + r/p)`, `call fraction ≈ prevalence · r/p`.
Dividing the inferred positive count by the inferred prevalence gives the implied
size of the scored subset, which is a useful consistency check.

| Board | Rank | Entrant | MCC | Precision | Recall | Inferred scored positives | Inferred prevalence | Inferred call fraction | Implied scored molecules |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Live | 1 | **N283T** | .5299 | .6284 | .7292 | 107.4 | .2916 | .338 | 368 |
| Live | 2 | harinu123 | .5154 | .6585 | .6542 | 100.4 | .2916 | .290 | 344 |
| Live | 3 | TeamPrescience | .5018 | .6493 | .6438 | 101.7 | .2916 | .289 | 349 |
| Live | 4 | EleanorRigby | .4987 | .5267 | .8701 | 97.9 | .2915 | .482 | 336 |
| Live | 5 | nova | .4938 | .6776 | .5873 | 108.3 | .2916 | .253 | 371 |
| Interim | 1 | EleanorRigby | .5207 | .5633 | .8623 | 234.4 | .3097 | .474 | 757 |
| Interim | 2 | CYPblend | .5134 | .6673 | .6591 | 235.3 | .3096 | .306 | 760 |
| Interim | 3 | stir_bar | .5035 | .5846 | .7883 | 232.1 | .3096 | .417 | 750 |
| Interim | 4 | **N283T** | .5023 | .6393 | .6848 | 230.4 | .3097 | .332 | 744 |
| Interim | 5 | srijitseal | .4981 | .5559 | .8376 | 226.2 | .3096 | .467 | 730 |

The full tables for the top 15 of each board are in `operating_points.csv`. Three
facts follow. First, every entrant's inferred prevalence is the same to three
decimals — .2916 on Live and .3096 on Interim — which is what should happen if the
identity is sound, because prevalence is a property of the scored population, not of
the entrant. Second, the implied scored-subset size is about 368 on Live and about
750 on Interim, so the Interim A4 board scores essentially every blinded molecule
while Live scores roughly half. Third, N283T's inferred call fraction is .338 (Live)
and .332 (Interim) against an exact full-pool fraction of 249/750 = .332, so the
inference reproduces the known quantity.

The blind A4 prevalence implied by the Interim board is therefore **.3097**, which
confirms rather than merely supports the ~29% working estimate and places it at the
upper end: the training labelled prevalence is .2132.

### A4 Live history

Each row is the earliest saved Live A4 snapshot for that submitted timestamp.
Snapshots were identified by their CSV content; endpoint boards carry a bare `MCC`
column while the macro board carries `MA-MCC`, and the N283T value was checked
against the known A4 set .4307 / .3750 / .4350 / .4496 / .5299 / .4843.

| Submitted | Entrants | Rank | MCC | Delta | Delta in std | Inferred call fraction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2026-09-10 13:32 | 78 | 18 | .4307 | — | — | .343 |
| 2026-09-12 12:27 | 81 | 45 | .3750 | −.0556 | −1.24 | .529 |
| 2026-09-13 11:38 | 81 | 15 | .4350 | +.0599 | +1.22 | .341 |
| 2026-09-14 02:08 | 82 | 13 | .4496 | +.0147 | +0.30 | .341 |
| 2026-09-15 13:00 | 90 | **1** | .5299 | +.0803 | **+1.69** | .338 |
| 2026-09-18 09:45 | 107 | 1 | .5299 | 0 | 0 | .338 |
| 2026-09-23 18:03 | 125 | 1 | .5299 | 0 | 0 | .338 |
| 2026-09-24 06:24 | 126 | 1 | .5299 | 0 | 0 | .338 |
| 2026-09-26 15:10 | 127 | 9 | .4843 | −.0456 | −0.95 | .373 |
| 2026-10-01 14:13 | 143 | 1 | .5299 | +.0456 | +0.96 | .338 |
| 2026-10-02 12:08 | 143 | 1 | .5299 | 0 | 0 | .338 |

Two things moved A4. The call rate moved it twice: the September 12 vector called
378 of 750 (inferred fraction .529) and lost .0556, and restoring 249 calls on
September 13 recovered .0599; the September 26 probe at 275 calls (inferred fraction
.373) lost .0456 and was reverted. The model moved it once decisively: the
September 14 → 15 swap of 18 calls for 18 at a constant 249 gained **+.0803**, the
only A4 transition larger than one reported bootstrap std. Everything after
September 15 is either an exact-vector control or the rejected rate probe.

### Verdict

**Live rank 1 is not distinguishable from the Live top tier.** The margin over rank
2 is +0.21 pooled bootstrap std and 15 other entrants lie within one pooled std of
N283T. The full-test picture, which is the one the Interim board measures, is rank 4
of 114 at .5023, −0.41 pooled std behind the leader and −0.03 behind rank 3 — again
indistinguishable, in the other direction. The honest statement is that the current
A4 vector is a strong top-tier entry on both populations, that its Live rank 1 is an
artifact of which half of the test set Live scores, and that it is not the best A4
model on the full test set as far as the published aggregates can resolve.

## Part 2 — Where does the submitted ensemble get its score?

The frozen complete-pool UMAP OOF of the submitted ensemble covers all 6,146
development molecules; 3,584 carry an A4 label and 764 of those are positive
(.2132). 2,335 labelled molecules have a direct pIC50 and 1,249 are assigned
negatives with none. 1,816 labelled molecules carry an A4 primary log2fc, so 1,768
do not. 2,207 unlabelled development molecules are primary hits (log2fc < −0.5) and
only 353 unlabelled screened molecules are non-hits, of which 161 sit at
log2fc ≥ −0.25. The highest log2fc of any A4 positive is −0.483, but unlike D6 the
labelled set extends well into the inactive region: 525 labelled molecules have
log2fc ≥ −0.5 and exactly one of them is positive.

### Strata, with the 249/750 contribution

Calls are the top 33.2% of each held fold's complete pool; the table reports only
the labelled rows. Nearest-blind similarity is the maximum Morgan r2/2048 Tanimoto
of each training molecule to any of the 750 blind molecules; quartile edges from the
labelled distribution are .2625 / .2958 / .3390.

| Stratification | Stratum | n | Positive rate | AUROC | AP | Calls | TP | Precision | Mean nearest-blind |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| primary log2fc | ≤ −3 | 399 | .599 | .625 | .728 | 252 | 167 | .663 | .335 |
| primary log2fc | (−3, −2] | 325 | .532 | .650 | .656 | 177 | 114 | .644 | .316 |
| primary log2fc | (−2, −1] | 329 | .505 | .595 | .586 | 150 | 85 | .567 | .308 |
| primary log2fc | (−1, −0.5) | 238 | .113 | .647 | .201 | 74 | 11 | .149 | .307 |
| primary log2fc | ≥ −0.5 | 525 | .002 | .947 | .034 | 49 | 1 | .020 | .305 |
| primary log2fc | none | 1768 | .089 | .853 | .415 | 283 | 98 | .346 | .311 |
| provenance | direct pIC50 | 2335 | .327 | .741 | .564 | 901 | 476 | .528 | .319 |
| provenance | assigned negative | 1249 | .000 | – | – | 84 | 0 | .000 | .300 |
| nearest-blind | Q1 | 904 | .179 | .772 | .426 | 231 | 92 | .398 | .236 |
| nearest-blind | Q2 | 902 | .192 | .832 | .536 | 236 | 111 | .470 | .280 |
| nearest-blind | Q3 | 888 | .221 | .801 | .524 | 262 | 128 | .489 | .316 |
| nearest-blind | Q4 | 890 | .262 | .817 | .601 | 256 | 145 | .566 | .420 |
| predicted log2fc | ≤ −3 | 814 | .154 | .857 | .501 | 164 | 78 | .476 | .319 |
| predicted log2fc | (−3, −2] | 1015 | .294 | .816 | .616 | 431 | 234 | .543 | .318 |
| predicted log2fc | (−2, −1] | 1068 | .258 | .722 | .475 | 362 | 155 | .428 | .308 |
| predicted log2fc | (−1, −0.5) | 411 | .134 | .664 | .269 | 27 | 8 | .296 | .303 |
| predicted log2fc | ≥ −0.5 | 276 | .036 | .591 | .187 | 1 | 1 | 1.000 | .305 |
| gate τ = −1 | p_active ≥ .5 | 2558 | .249 | .808 | .557 | 879 | 437 | .497 | .316 |
| gate τ = −1 | p_active < .5 | 1026 | .123 | .727 | .306 | 106 | 39 | .368 | .304 |
| gate τ = −2 | p_active ≥ .5 | 1747 | .226 | .842 | .583 | 553 | 287 | .519 | .319 |
| gate τ = −2 | p_active < .5 | 1837 | .201 | .765 | .452 | 432 | 189 | .438 | .306 |

The whole-population AUROC is .806 and AP .524. The decomposition shows that almost
none of that comes from resolving which strong primary actives are TDI positives:
inside the three strongly active bins the AUROC is .595 to .650, barely above chance
and far below the pooled figure. What the ensemble does well is coarse: it ranks the
525 weak-primary molecules correctly (AUROC .947, one positive, 49 calls), it ranks
the 1,768 molecules with no primary measurement well (AUROC .853, AP .415), and it
places only 84 of its 985 labelled calls on the 1,249 assigned negatives. The pooled
AUROC is largely a between-stratum effect.

Performance does not degrade on blind-like chemistry. Q4, the most blind-like
quartile, has the highest positive rate (.262), the highest AUROC among the
quartiles (.817) and the best precision at the operating point (.566); Q1 is the
worst on every count. All four quartiles have a nearly identical share of the calls.

### Reweighting to the blind population

Weights are `p_blind(stratum)/p_labelled(stratum)`. Blind-side mixes come from two
independent estimators: `nn1`/`nn5` transfer each blind molecule the stratum of its
nearest (or similarity-weighted five nearest) development molecule, and `direct`
applies the fold-excluded log2fc gate to the blind molecules themselves and bins both
sides by the *predicted* log2fc so no stratum is empty. The gate regressor's inner
Spearman against observed log2fc is .693 to .717 across folds, much better than the
.48 to .54 the D6 version reached.

| Weighting | Effective n | Projected prevalence | AUROC | AP | 200 | 225 | 249 | 275 | 300 | 325 | 350 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| none (labelled OOF) | 3584 | .213 | .806 | .524 | .375 | .386 | .406 | **.427** | .423 | .420 | .401 |
| primary log2fc : nn1 | 1152 | .471 | .708 | .674 | .274 | .276 | .307 | .327 | .328 | **.329** | .321 |
| primary log2fc : nn5 | 1988 | .386 | .740 | .626 | .311 | .318 | .347 | .370 | .371 | **.374** | .360 |
| provenance : nn1 | 2367 | .325 | .742 | .563 | .307 | .317 | .341 | **.368** | .364 | .365 | .347 |
| provenance : nn5 | 2841 | .293 | .763 | .554 | .330 | .341 | .365 | **.390** | .386 | .387 | .369 |
| nearest-blind : nn1 | 890 | .262 | .817 | .601 | .423 | .439 | .440 | **.457** | .455 | .443 | .419 |
| nearest-blind : nn5 | 978 | .260 | .816 | .597 | .420 | .436 | .439 | **.455** | .454 | .442 | .419 |
| predicted log2fc : direct | 3529 | .216 | .810 | .531 | .383 | .394 | .414 | **.432** | .428 | .426 | .406 |
| predicted log2fc : nn1 | 2836 | .231 | .819 | .553 | .400 | .413 | .433 | **.445** | .441 | .439 | .419 |
| gate τ = −1 : direct | 3542 | .219 | .807 | .530 | .377 | .389 | .409 | **.428** | .425 | .423 | .404 |
| gate τ = −1 : nn1 | 2946 | .240 | .809 | .550 | .380 | .395 | .414 | **.432** | .431 | .429 | .410 |
| gate τ = −2 : direct | 3539 | .215 | .811 | .531 | .382 | .394 | .414 | **.432** | .429 | .427 | .407 |
| gate τ = −2 : nn1 | 2898 | .219 | .824 | .554 | .407 | .421 | .441 | **.451** | .449 | .446 | .425 |

The `nn5` rows for the gate weightings are in `reweighted_metrics.csv` and sit
between their `direct` and `nn1` counterparts.

Two conclusions, both negative for the current operating point.

**The 29% prevalence estimate is bracketed but not reproduced.** The
leaderboard-derived blind prevalence is .2916 (Live) and .3097 (Interim). The
activity-gate and predicted-activity reweightings project only .215 to .240, i.e.
they say the blind set is barely more positive than the training pool. The
observed-log2fc nearest-neighbour transfer projects .386 to .471, far above the
leaderboard figure. Only the label-provenance transfer lands on it (.293 and .325).
The spread of .215 to .471 across reasonable estimators means no reweighting here
pins the blind prevalence; the .29–.31 figure derived from the public aggregates is
the better estimate, and it is consistent with the plausible range.

**The 249-call choice is not the reweighted OOF optimum under any weighting
tested.** Every weighting except the two observed-log2fc transfers peaks at 275, and
those two peak at 325. Projected MCC at 249 is .406 unweighted against .427 at 275,
and the gap persists after every reweighting. This reproduces the September 26 sweep
rather than correcting it, and the September 26 Live probe at 275 lost .0456. The
reweighted OOF therefore does not explain why 249 wins on the leaderboard; the
disagreement between the complete-pool OOF and the Live board on call rate survives
population reweighting and remains unexplained.

**The advantage is concentrated in weakly active and unmeasured strata, not in
strongly active ones.** Within-stratum AUROC falls from .947 (log2fc ≥ −0.5) and
.853 (no primary) to .595–.650 in the three strong bins. Reweighting towards the
observed-activity mix of the blind set drops pooled AUROC from .806 to .708–.740.
Reweighting by chemical proximity instead leaves it at .816–.824. The blind set is
built from analogs of strong 3A4/1A2/2C9 hits, so the activity composition shifts
towards exactly the region where the ensemble is weakest, while the chemistry stays
in a region where it is fine.

On the blind side the submitted 249 are already strongly enriched for predicted
primary activity: their mean p_active(−1) is .827 against .677 for the 501
non-calls, and mean predicted log2fc is −2.51 against −2.07. Their mean nearest
development Tanimoto is .596 against .599, so the submitted calls are not simply
the molecules closest to the training set.

## Part 3 — Can the ensemble be improved?

Fitted inside each outer fold on the fixed five-fold `umap` and `scaffold`
assignments, with nothing from the held fold entering any fit; features are fresh
Morgan r2/2048 plus 217 RDKit descriptors (2,265 columns) and frozen CheMeleon
2,048-d; LightGBM is fixed at 400 trees, 31 leaves, lr 0.03, colsample 0.5,
min_child_samples 20, balanced class weights, two seeds averaged; the CheMeleon
logistic penalty C = 0.01 was selected by inner AP in all eleven fits.

| Model | Recipe |
| --- | --- |
| A `official_lgbm` | balanced Morgan+RDKit LightGBM on official A4 labels |
| B `pool_lgbm` | A plus the 161 unlabelled screened non-hits (log2fc ≥ −0.25) as pseudo-negatives |
| C `gate_x_tdi` | p_active = P(log2fc < −1) from a Platt-calibrated LightGBM log2fc regressor, times A |
| D `chemeleon_lr` | balanced CheMeleon logistic regression on official labels |
| E `submitted_equal7` | comparator: frozen complete-pool equal-seven OOF, **umap only** |

The comparator's alignment with this run was verified over all 6,146 rows: 0 missing
scores, 0 fold mismatches, 0 label mismatches. E has no scaffold counterpart, so the
scaffold family uses A as its reference and the deltas against E are UMAP-only.

Pool composition on A4 differs sharply from D6. The official labels already include
the primary-inactive population (525 labelled molecules at log2fc ≥ −0.5, one
positive), only 161 pseudo-negatives exist, and 66.5% of screened development
molecules satisfy log2fc < −1, against 27.4% for D6. Pool positive rate .2040 against
official .2132. Evaluation is on official labels only.

### Ranking and decision metrics, official labels, pooled held folds

| Split | Model | AUROC | AP | 25% | 29% | 33.2% | 37% | 40% | Nested MCC |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| umap | E submitted_equal7 | **.806** | **.524** | .370 | .387 | .406 | .429 | .423 | **.429** |
| umap | equal7_filter_0.1 | .807 | .522 | .372 | .387 | .406 | .429 | .423 | .429 |
| umap | equal7_filter_0.2 | .805 | .521 | .370 | .385 | .405 | .426 | .420 | .426 |
| umap | equal7_filter_0.3 | .794 | .516 | .368 | .383 | .410 | .418 | .418 | .404 |
| umap | rank_mean(E, A) | .802 | .508 | .380 | .390 | .389 | .404 | .403 | .397 |
| umap | equal7 × sqrt(p_active) | .798 | .504 | .381 | .389 | .385 | .393 | .396 | .372 |
| umap | rank_mean(E, D) | .782 | .494 | .365 | .367 | .376 | .365 | .362 | .376 |
| umap | equal7 × p_active | .781 | .483 | .368 | .359 | .361 | .367 | .360 | .354 |
| umap | C gate_x_tdi | .769 | .457 | .317 | .324 | .330 | .332 | .342 | .342 |
| umap | B pool_lgbm | .769 | .454 | .337 | .338 | .343 | .344 | .336 | .322 |
| umap | A official_lgbm | .767 | .449 | .332 | .329 | .330 | .329 | .330 | .312 |
| umap | D chemeleon_lr | .698 | .370 | .221 | .228 | .237 | .231 | .236 | .233 |
| scaffold | A official_lgbm | .808 | .509 | .383 | .383 | .387 | .391 | .395 | .389 |
| scaffold | B pool_lgbm | .806 | .507 | .372 | .384 | .389 | .391 | .390 | .383 |
| scaffold | C gate_x_tdi | .809 | .513 | .377 | .377 | .389 | .390 | .387 | .360 |
| scaffold | D chemeleon_lr | .765 | .430 | .285 | .306 | .320 | .330 | .338 | .338 |

The submitted ensemble leads every new model on UMAP by .037 to .108 AUROC and .067
to .154 AP. Its own rate curve rises to 37% (.429) and only then falls, consistent
with Part 2 and with the September 26 sweep. The gate filter at τ = 0.1 is effectively a
no-op: no blind molecule and very few development molecules fall below p_active 0.1,
so the variant reproduces E to three decimals and its calls are identical at the
reference rate.

### Paired group bootstrap against E, official labels, 1,000 group resamples

| Variant | 33.2% delta [95%] | Nested delta [95%] | Blind overlap with the 249 |
| --- | --- | --- | ---: |
| equal7_filter_0.1 | +.0000 [.0000, .0000] | +.0000 [−.0012, +.0012] | 249 |
| equal7_filter_0.2 | −.0006 [−.0041, +.0025] | −.0034 [−.0079, +.0005] | 247 |
| equal7_filter_0.3 | **+.0040 [−.0050, +.0120]** | −.0257 [−.0427, −.0097] | 243 |
| rank_mean(E, A) | −.0167 [−.0392, +.0041] | −.0322 [−.0551, −.0103] | 217 |
| equal7 × sqrt(p_active) | −.0210 [−.0486, +.0096] | −.0569 [−.0833, −.0318] | 216 |
| rank_mean(E, D) | −.0304 [−.0654, +.0053] | −.0537 [−.0888, −.0189] | 212 |
| equal7 × p_active | −.0452 [−.0786, −.0082] | −.0752 [−.1052, −.0429] | 207 |
| B pool_lgbm | −.0631 [−.0996, −.0288] | −.1078 [−.1439, −.0750] | 198 |
| A official_lgbm | −.0761 [−.1167, −.0418] | −.1169 [−.1511, −.0809] | 194 |
| C gate_x_tdi | −.0758 [−.1187, −.0426] | −.0869 [−.1202, −.0555] | 186 |
| D chemeleon_lr | −.1688 [−.2164, −.1266] | −.1966 [−.2362, −.1592] | 178 |

**No variant has a positive delta whose interval excludes zero at either operating
point.** The largest positive point estimate anywhere in the table is +.0040 for the
τ = 0.3 filter at the reference rate, with an interval spanning zero, and that same
variant is significantly negative under nested rate selection
(−.0257 [−.0427, −.0097]). The τ = 0.1 filter is identically zero because it changes
nothing. Every other variant is negative, and all four standalone new models are
significantly worse than the comparator.

Model B versus model A, the pseudo-negative question, behaves as the pool
composition predicts: +.0130 [−.0043, +.0309] at the reference rate and
+.0090 [−.0128, +.0315] nested on UMAP, +.0016 [−.0135, +.0186] and
−.0061 [−.0237, +.0120] on scaffold. Adding 161 pseudo-negatives to 3,584 official
labels changes nothing measurable, which is the expected result and the opposite of
the large D6 effect, because the A4 official labels already cover the primary-inactive
population. The gate as a multiplier helps A slightly on UMAP under nested selection
(+.0300 [+.0025, +.0548]) and hurts it on scaffold (−.0295 [−.0529, −.0083]), so
that gain does not transfer across split families.

### Blind activity gate

Refitted on all 6,146 development molecules and applied to the 750 blind molecules,
the gate gives mean p_active(−1) .727, median .861, quartiles .525 and .969, with
572 of 750 above 0.5 and only 13 below 0.1, 68 below 0.2 and 106 below 0.3. The
development base rate of log2fc < −1 among screened molecules is 66.5%. The gate
therefore has very little to say about the blind A4 set: it judges most of it
primary-active, which is consistent with the blind set's construction and leaves the
filters with almost nothing to remove. This is the structural reason the A4 gate
cannot reproduce the D6 gate's usefulness.

### Historical A4 transition concordance

Mean blind score over the calls removed from the earlier payload and over the calls
added by the later one, with a two-sided permutation p over 10,000 relabelings of
the union (fixed seed). Payloads: `20260912T122715261670Z-363bf7b8` (378 calls),
`20260913T113807374506Z-77f4825a`, `20260914T020759952207Z-37dd7bd8` and
`20260915T130004513006Z-75b25070` (249 calls each); transition sizes 131/2 and
18/18 were reproduced from the payload files.

| Transition | LB delta | Column | removed → added | Difference | perm p | Sign agrees |
| --- | ---: | --- | --- | ---: | ---: | --- |
| sep12 → sep13 | +.0599 | E submitted_equal7 | .541 → .725 | +.184 | .0041 | yes* |
| sep12 → sep13 | +.0599 | gate p_active | .660 → .869 | +.209 | .354 | yes |
| sep12 → sep13 | +.0599 | A official_lgbm | .277 → .306 | +.029 | .849 | yes |
| sep12 → sep13 | +.0599 | B pool_lgbm | .280 → .280 | +.000 | .999 | yes |
| sep12 → sep13 | +.0599 | C gate_x_tdi | .204 → .262 | +.058 | .697 | yes |
| sep12 → sep13 | +.0599 | D chemeleon_lr | .425 → .647 | +.222 | .243 | yes |
| sep12 → sep13 | +.0599 | equal7 × p_active | .362 → .656 | +.294 | .0070 | yes |
| sep12 → sep13 | +.0599 | equal7 × sqrt(p_active) | .433 → .704 | +.271 | .0008 | yes |
| sep12 → sep13 | +.0599 | rank_mean(E, A) | .544 → .680 | +.136 | .125 | yes |
| sep12 → sep13 | +.0599 | rank_mean(E, D) | .542 → .743 | +.201 | .034 | yes |
| sep14 → sep15 | +.0803 | E submitted_equal7 | .608 → .655 | +.047 | .0000 | yes* |
| sep14 → sep15 | +.0803 | gate p_active | .760 → .749 | −.011 | .895 | **no** |
| sep14 → sep15 | +.0803 | A official_lgbm | .442 → .295 | −.147 | .011 | **no** |
| sep14 → sep15 | +.0803 | B pool_lgbm | .448 → .284 | −.164 | .0076 | **no** |
| sep14 → sep15 | +.0803 | C gate_x_tdi | .335 → .244 | −.091 | .127 | **no** |
| sep14 → sep15 | +.0803 | D chemeleon_lr | .496 → .441 | −.055 | .523 | **no** |
| sep14 → sep15 | +.0803 | equal7 × p_active | .474 → .520 | +.047 | .414 | yes |
| sep14 → sep15 | +.0803 | equal7 × sqrt(p_active) | .538 → .588 | +.050 | .173 | yes |
| sep14 → sep15 | +.0803 | rank_mean(E, A) | .665 → .624 | −.041 | .154 | **no** |
| sep14 → sep15 | +.0803 | rank_mean(E, D) | .610 → .618 | +.008 | .835 | yes |

The filter columns are in `transition_concordance.csv`; τ = 0.1 is identical to E and
τ = 0.2 / 0.3 give +.011 at p ≈ 1 on the second transition.

\* **The comparator's own agreement is tautological, not evidence.** The September 15
call set *is* the top 249 of the `a4_prediction` column, so the 18 added calls lie
inside that top 249 by definition and the 18 removed calls lie outside it, giving
perfect separation and a permutation p of zero. The same applies to the τ = 0.1
filter, which reproduces E exactly. Any variant that is a mild monotone perturbation
of E inherits part of this. The first transition is also weak evidence: removing 131
of 378 calls drops the lowest-ranked part of an over-called vector, so almost any
score agrees in sign, and indeed all ten columns do, with A and B showing no
separation at all (+.029 at p = .85 and +.000 at p = .999).

The informative row is the September 14 → 15 swap, and there the result is uniformly
against the new models. The gate, A, B, C and D all rank the September 15 additions
*below* the September 14 removals, two of them significantly, while the move gained
+.0803 on the Live board — the largest A4 gain on record. Mixing A into E by rank
mean inverts E's direction (−.041). Only the two product forms and rank_mean(E, D)
keep the right sign, and neither separation is significant.

## Plain-language answers

**Why is A4 winning on Live, and is the lead real?** The lead is not real in any
statistical sense. On the September 26 Live board N283T is +0.21 pooled bootstrap
std ahead of rank 2 with 15 other entrants inside one pooled std; on the Interim
board, which scores the whole test set, the same predictions place rank 4 at −0.41
pooled std behind the leader. A4's Live-to-Interim MCC Spearman is .756 with only
3 of 10 top names shared, against .915 for D6 and .898 for the macro board, so the
A4 half-test ranking is the noisiest of the three. Rank 1 on Live is a favourable
draw of which 368 of the 750 molecules Live happens to score, and the honest
statement is "strong top-tier A4 entry on both populations", not "best A4 model".

**What the score is actually made of.** The ensemble's pooled .806 AUROC is mostly
coarse triage. It is near-perfect at pushing down the 525 labelled molecules with no
meaningful primary 3A4 inhibition (AUROC .947), good on the 1,768 labelled molecules
with no primary measurement at all (.853), and it wastes only 84 of 985 labelled
calls on the 1,249 assigned negatives. Inside the strongly active population, where
the positive rate is 50–60%, it is close to uninformative (AUROC .595 to .650). The
blind set was built as analog expansions around strong 3A4/1A2/2C9 hits and its
scored prevalence is .3097 on the Interim board against .2132 in training, so the
test population is shifted precisely into the region the ensemble resolves worst.
Reweighting the labelled OOF to the blind activity mix drops AUROC to .708–.740;
reweighting by chemical proximity does not drop it at all (.816–.824). Both can be
true: the chemistry transfers, the difficulty mix does not.

**Is 249 defensible?** Not from the OOF. Every reweighting tested, including the two
built from the activity gate, puts the projected optimum at 275 or 325, not 249, and
the September 26 Live probe at 275 lost .0456 and was reverted. The OOF-versus-Live
disagreement on A4 call rate is reproduced, not resolved, by this work. 249 is
retained because it is the leaderboard-validated count, not because the OOF endorses
it.

**Is any variant a credible improvement?** No. The decision rule was a positive
paired delta against E on UMAP whose interval excludes zero, not contradicted by the
historical transitions. Nothing clears the first half: the best point estimate is
+.0040 [−.0050, +.0120] for the τ = 0.3 gate filter at the reference rate, and that
variant is −.0257 [−.0427, −.0097] under nested rate selection. The second half fails
too: on the one informative historical transition the gate and all four standalone
models rank the September 15 winners below the September 14 losers. **Not
candidate-worthy** — none of the eleven variants, and none of A, B, C or D on their
own.

## Limitations

- The published metrics are separately averaged bootstrap estimates. Inferred
  positive counts, prevalences and call fractions are approximations, not other
  teams' exact vectors, even though the internal consistency checks pass.
- The reported `MCC std` is the bootstrap spread of one fixed prediction vector on
  one scored population. Expressing between-entrant gaps in those units is a
  reasonable resolution scale, not a hypothesis test, and the pooled form assumes an
  independence that two entrants scored on the same molecules do not have.
- Live and Interim are one observation each. The Interim board was captured once, on
  September 26, for the September 24 resubmission of the September 15 vector.
- The equal-seven OOF is consumed frozen, so the small TabICL/TabPFN reload
  differences documented for run `e74748cd` (Spearman .99999989 against the saved
  labelled-only scores, maximum difference .000968) are inherited here.
- Blind stratum mixes are inferred from training chemistry and a fitted log2fc
  model, never measured. The spread across estimators (.215 to .471 projected
  prevalence) is itself the main uncertainty, and importance weights cannot correct
  for a stratum the labelled pool does not cover.
- The gate reweightings bin both sides by *predicted* log2fc, so they measure
  predicted-activity transfer rather than true-activity transfer.
- The comparator exists on UMAP folds only. Every delta against E is single-family,
  and the one cross-family check available (C against A) changes sign between
  families, which is a direct warning about single-family deltas.
- The fixed-rate intervals are conditional on these OOF predictions and are not
  adjusted for having inspected a five-point rate grid.
- Two A4 historical transitions is not a concordance benchmark. One is tautological
  for the comparator and the other is near-automatic in sign for every score, so the
  check can refute a variant but cannot endorse one.
- No hyperparameter search was run; LightGBM settings are inherited from the D6
  probe, and only the CheMeleon logistic penalty was selected inside each fold.
- Nothing here addresses the blind scored-label mask, i.e. which of the 750
  molecules are scored on either board.
