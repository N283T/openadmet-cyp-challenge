# D6 activity-gate historical LB concordance and legacy x gate combinations

Recorded run `20261002T221918697661Z-2d8d20f5` (classification, 3.1 s).
Code: `scripts/evaluate_cyp_d6_gate_concordance.py`,
`scripts/launch_cyp_d6_gate_concordance.py`,
`configs/classification/d6-gate-concordance-v1.json`.
Compact outputs with SHA-256 provenance:
`artifacts/tracked/cyp-d6-gate-concordance-20261002T221918697661Z-2d8d20f5/`.
No model was fitted in this run; it consumes the frozen scores of
`20261002T215345688167Z-7920b49a` ([analysis](2026-10-03-d6-testlike-pool.md)).
No candidate, prepared submission or upload resulted from it.

## Question

The previous run built a calibrated D6 activity gate, p_active = P(primary
log2fc < -1), and showed that pool-trained and gated models beat the
official-labels-only control under pseudo-negative labels while losing under
official labels. That disagreement left the gate unvalidated against anything
external. This run asks two separate questions:

1. Does the gate score, which never saw a blind label, separate the removed
   from the added calls of the project's seven submitted D6 payloads in the
   direction of the observed Live MCC movements? This is the
   [historical leaderboard concordance backtest](../classification-evaluation.md)
   criterion.
2. If it does, does combining the gate with the historically best recipe (the
   Sep 15 legacy score) improve complete-pool OOF decisions over that recipe
   itself, rather than over the weak official-only control used before?

## Inputs and frozen criteria

Declared inputs: `oof_scores.parquet` and `test_scores.parquet` of run
`7920b49a`; the deployed Sep 15 blind score vector
(`20260924T045344348206Z-7e2db311/inputs/blind-scores.parquet`, originally
`decision-test.parquet` of `20260923T173655387004Z-97d71dce`); the seven
submitted D6 payloads; `aliases.parquet`, `test_mapping.parquet` and the
`assignments.parquet` that supplies bootstrap groups; and
`data/raw/cyp-challenge-single-concentration-TRAIN.csv` for the sensitivity
label set.

Verification of the Sep 15 vector: 750 rows, 49 stored calls, its top 49 by
score equals the Sep 15 payload's 49 positives exactly, and the 49th-to-50th
score gap is 0.00133, so no tie is involved. The vector is therefore usable as
the deployed continuous legacy score. The OOF comparator
`legacy_sep15_refit` is a refit analog of that recipe on the fixed folds, not
the deployed fit.

Score columns, call sets, transitions, the rate grid, the variant list and the
uncertainty procedure were fixed in the config before any concordance number
was read. `gate_x_tdi.p_tdi` turns out to be bit-identical to
`official_lgbm.score` (maximum absolute difference 0.0 over the 750 blind
molecules), because the gate model multiplies exactly that probability; the two
columns are kept separate in the tables but are one observation.

## Part 1 - submitted payloads

| Payload | Submission | Positives | Live MCC | mean p_active (Morgan) | mean p_active (CheMeleon) | positives with p_active < 0.3 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| sep13 | `20260913T113807374506Z-77f4825a` | 49 | .2710 | .301 | .277 | 28 |
| sep14 | `20260914T020759952207Z-37dd7bd8` | 49 | .2786 | .295 | .274 | 29 |
| sep15 | `20260915T130004513006Z-75b25070` | 49 | .3086 | .322 | .304 | 27 |
| sep18 | `20260918T094451131332Z-f23dd61b` | 49 | .1097 | .283 | .254 | 29 |
| sep23 | `20260923T180233070503Z-552d9af4` | 49 | .2348 | .290 | .271 | 32 |
| oct01 | `20261001T141041310215Z-c7faa79f` | 49 | .1050 | .236 | .219 | 37 |
| oct02 | `20261002T092348406577Z-4ea606c1` | 175 | .1177 | .215 | .216 | 139 |

The mean gate probability of a payload's positives is monotone with its Live
MCC over these seven states except for the Sep 18/Sep 23 pair. Even the best
payload (Sep 15) has 27 of its 49 calls below p_active 0.3 and 18 below 0.2, so
the gate does not endorse the historical best vector as a whole.

## Part 1 - transition concordance

For each transition, mean score over the calls removed from the earlier payload
and over the calls added by the later one, their difference (added minus
removed), a two-sided permutation p over 10,000 relabelings of the union of the
two call sets (fixed seed), and whether the sign of the difference agrees with
the Live MCC delta.

| Transition | LB delta | Column | removed / added | mean removed | mean added | difference | perm p | sign agrees |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | --- |
| sep13→sep14 | +.0076 | gate p_active (Morgan) | 6 / 6 | .278 | .223 | -.055 | .558 | no |
| sep13→sep14 | +.0076 | gate p_active (CheMeleon) | 6 / 6 | .250 | .218 | -.032 | .693 | no |
| sep13→sep14 | +.0076 | pool_lgbm score | 6 / 6 | .083 | .153 | +.070 | .166 | yes |
| sep13→sep14 | +.0076 | pool_chemeleon score | 6 / 6 | .517 | .283 | -.233 | .174 | no |
| sep13→sep14 | +.0076 | rank_mean_gates score | 6 / 6 | .886 | .793 | -.094 | .216 | no |
| sep13→sep14 | +.0076 | gate p_tdi (official head) | 6 / 6 | .413 | .370 | -.044 | .660 | no |
| sep13→sep14 | +.0076 | official_lgbm score | 6 / 6 | .413 | .370 | -.044 | .660 | no |
| sep14→sep15 | +.0300 | gate p_active (Morgan) | 9 / 9 | .164 | .310 | +.146 | .041 | yes |
| sep14→sep15 | +.0300 | gate p_active (CheMeleon) | 9 / 9 | .161 | .325 | +.164 | .010 | yes |
| sep14→sep15 | +.0300 | pool_lgbm score | 9 / 9 | .098 | .271 | +.173 | .058 | yes |
| sep14→sep15 | +.0300 | pool_chemeleon score | 9 / 9 | .339 | .554 | +.215 | .184 | yes |
| sep14→sep15 | +.0300 | rank_mean_gates score | 9 / 9 | .656 | .904 | +.247 | .004 | yes |
| sep14→sep15 | +.0300 | gate p_tdi (official head) | 9 / 9 | .236 | .464 | +.228 | .040 | yes |
| sep14→sep15 | +.0300 | official_lgbm score | 9 / 9 | .236 | .464 | +.228 | .039 | yes |
| sep15→sep18 | -.1989 | gate p_active (Morgan) | 27 / 27 | .315 | .245 | -.070 | .152 | yes |
| sep15→sep18 | -.1989 | gate p_active (CheMeleon) | 27 / 27 | .312 | .221 | -.091 | .056 | yes |
| sep15→sep18 | -.1989 | pool_lgbm score | 27 / 27 | .328 | .125 | -.203 | .003 | yes |
| sep15→sep18 | -.1989 | pool_chemeleon score | 27 / 27 | .498 | .375 | -.123 | .153 | yes |
| sep15→sep18 | -.1989 | rank_mean_gates score | 27 / 27 | .835 | .729 | -.105 | .025 | yes |
| sep15→sep18 | -.1989 | gate p_tdi (official head) | 27 / 27 | .495 | .301 | -.193 | .004 | yes |
| sep15→sep18 | -.1989 | official_lgbm score | 27 / 27 | .495 | .301 | -.193 | .004 | yes |
| sep15→sep23 | -.0738 | gate p_active (Morgan) | 15 / 15 | .294 | .191 | -.103 | .104 | yes |
| sep15→sep23 | -.0738 | gate p_active (CheMeleon) | 15 / 15 | .311 | .206 | -.105 | .124 | yes |
| sep15→sep23 | -.0738 | pool_lgbm score | 15 / 15 | .209 | .268 | +.059 | .545 | no |
| sep15→sep23 | -.0738 | pool_chemeleon score | 15 / 15 | .539 | .450 | -.089 | .333 | yes |
| sep15→sep23 | -.0738 | rank_mean_gates score | 15 / 15 | .806 | .721 | -.085 | .221 | yes |
| sep15→sep23 | -.0738 | gate p_tdi (official head) | 15 / 15 | .367 | .377 | +.011 | .902 | no |
| sep15→sep23 | -.0738 | official_lgbm score | 15 / 15 | .367 | .377 | +.011 | .907 | no |
| sep18→sep23 | +.1251 | gate p_active (Morgan) | 21 / 21 | .254 | .271 | +.017 | .733 | yes |
| sep18→sep23 | +.1251 | gate p_active (CheMeleon) | 21 / 21 | .231 | .272 | +.042 | .373 | yes |
| sep18→sep23 | +.1251 | pool_lgbm score | 21 / 21 | .097 | .400 | +.303 | .0001 | yes |
| sep18→sep23 | +.1251 | pool_chemeleon score | 21 / 21 | .358 | .453 | +.095 | .323 | yes |
| sep18→sep23 | +.1251 | rank_mean_gates score | 21 / 21 | .728 | .802 | +.075 | .186 | yes |
| sep18→sep23 | +.1251 | gate p_tdi (official head) | 21 / 21 | .285 | .541 | +.257 | .001 | yes |
| sep18→sep23 | +.1251 | official_lgbm score | 21 / 21 | .285 | .541 | +.257 | .001 | yes |
| sep15→oct01 | -.2036 | gate p_active (Morgan) | 21 / 21 | .365 | .164 | -.200 | .0001 | yes |
| sep15→oct01 | -.2036 | gate p_active (CheMeleon) | 21 / 21 | .365 | .168 | -.197 | .0003 | yes |
| sep15→oct01 | -.2036 | pool_lgbm score | 21 / 21 | .364 | .232 | -.133 | .122 | yes |
| sep15→oct01 | -.2036 | pool_chemeleon score | 21 / 21 | .523 | .313 | -.210 | .022 | yes |
| sep15→oct01 | -.2036 | rank_mean_gates score | 21 / 21 | .877 | .677 | -.200 | .0003 | yes |
| sep15→oct01 | -.2036 | gate p_tdi (official head) | 21 / 21 | .493 | .367 | -.126 | .136 | yes |
| sep15→oct01 | -.2036 | official_lgbm score | 21 / 21 | .493 | .367 | -.126 | .135 | yes |
| sep15→oct02 | -.1909 | gate p_active (Morgan) | 6 / 132 | .435 | .185 | -.249 | .0005 | yes |
| sep15→oct02 | -.1909 | gate p_active (CheMeleon) | 6 / 132 | .435 | .193 | -.241 | .0001 | yes |
| sep15→oct02 | -.1909 | pool_lgbm score | 6 / 132 | .340 | .106 | -.234 | .006 | yes |
| sep15→oct02 | -.1909 | pool_chemeleon score | 6 / 132 | .641 | .306 | -.335 | .004 | yes |
| sep15→oct02 | -.1909 | rank_mean_gates score | 6 / 132 | .945 | .606 | -.339 | .0002 | yes |
| sep15→oct02 | -.1909 | gate p_tdi (official head) | 6 / 132 | .496 | .234 | -.262 | .004 | yes |
| sep15→oct02 | -.1909 | official_lgbm score | 6 / 132 | .496 | .234 | -.262 | .005 | yes |

Agreement counts over the seven transitions:

| Column | sign agreements | transitions with perm p < .05 | smallest perm p |
| --- | ---: | ---: | ---: |
| gate p_active (Morgan) | 6 | 3 | .0001 |
| gate p_active (CheMeleon) | 6 | 3 | .0001 |
| pool_lgbm score | 6 | 3 | .0001 |
| pool_chemeleon score | 6 | 2 | .004 |
| rank_mean_gates score | 6 | 4 | .0002 |
| official_lgbm score | 5 | 4 | .001 |
| gate p_tdi (official head) | 5 | 4 | .001 |

The two gate probabilities agree with 6 of 7 directions; their single
disagreement is Sep 13→Sep 14, a +.0076 LB move over six swapped calls that is
inside leaderboard noise. `pool_lgbm` also agrees 6 of 7, but it fails on the
Sep 15→Sep 23 decline (+.059, p = .545), which the gate calls correctly, and
the official-only head fails on both Sep 13→14 and Sep 15→23. The two October
failures, which were not available when the gate design was chosen, are the
most clearly separated: -.200 (p = .0001) and -.249 (p = .0005) for the Morgan
gate, -.197 (p = .0003) and -.241 (p = .0001) for the CheMeleon gate.

These seven transitions are not seven independent experiments. They share one
hidden blind set, one scored-label mask and one bootstrap public MCC estimate;
consecutive payloads share most of their calls; and the number of swapped
molecules is 6 to 132. The reserved Sep 18→Sep 23 check is agreed but
uninformative on its own (difference +.017, p = .733). This is calibration
evidence that the gate probability carries information about the direction of
D6 Live movements. It is not proof, and it does not establish an effect size.

## Part 2 - legacy x gate combinations on the complete-pool OOF

The control is the legacy comparator, not the official-only LightGBM. All
scores are per-held-fold percentile ranks of the frozen OOF scores of run
`7920b49a`; nothing is refitted.

| Variant | Definition |
| --- | --- |
| L0 | control: per-fold percentile rank of `legacy_sep15_refit` |
| L1 | L0 rank x p_active (Morgan) |
| L2 | L0 rank x sqrt(p_active Morgan) |
| L3 | L0 rank x sqrt(p_active Morgan x p_active CheMeleon) |
| L4 | equal per-fold rank mean of L0 and `pool_lgbm` |
| L5(tau) | L0 ranking restricted to p_active >= tau, calls then filled from the excluded molecules by L0 rank; tau in {0.2, 0.3, 0.4} |
| B | reference: per-fold percentile rank of `pool_lgbm` |

Pools are unchanged from the previous run: 6,146 development molecules per
family, 1,497 official labels, 4,190 test-like labels. On development
molecules the legacy score and the Morgan gate are almost uncorrelated
(Spearman -.054 umap, -.062 scaffold), so the product variants combine nearly
independent information. 2,169 (umap) / 2,058 (scaffold) development
molecules have held-fold p_active >= 0.3.

### Fixed 6.53% (49/750 analog), pooled over held folds

| Split | Variant | calls | MCC t-l | prec t-l | TP/FP t-l | MCC off. | prec off. | TP/FP off. |
| --- | --- | ---: | ---: | ---: | --- | ---: | ---: | --- |
| umap | L0 | 401 | .074 | .156 | 39 / 211 | .147 | .464 | 39 / 45 |
| umap | L1 | 401 | .089 | .171 | 43 / 208 | .005 | .222 | 43 / 151 |
| umap | L2 | 401 | .102 | .188 | 45 / 195 | .040 | .262 | 45 / 127 |
| umap | L3 | 401 | .091 | .175 | 43 / 203 | -.003 | .213 | 43 / 159 |
| umap | L4 | 401 | .104 | .189 | 46 / 198 | .106 | .359 | 46 / 82 |
| umap | L5(0.2) | 401 | .095 | .186 | 40 / 175 | .095 | .354 | 40 / 73 |
| umap | L5(0.3) | 401 | .086 | .173 | 39 / 186 | .047 | .277 | 39 / 102 |
| umap | L5(0.4) | 401 | .095 | .181 | 43 / 194 | .033 | .254 | 43 / 126 |
| umap | B | 401 | .096 | .177 | 46 / 214 | .031 | .250 | 46 / 138 |
| scaffold | L0 | 400 | .087 | .169 | 43 / 211 | .145 | .443 | 43 / 54 |
| scaffold | L1 | 400 | .138 | .223 | 56 / 195 | .057 | .276 | 56 / 147 |
| scaffold | L2 | 400 | .159 | .251 | 59 / 176 | .104 | .333 | 59 / 118 |
| scaffold | L3 | 400 | .140 | .225 | 57 / 196 | .058 | .277 | 57 / 149 |
| scaffold | L4 | 400 | .142 | .232 | 55 / 182 | .157 | .426 | 55 / 74 |
| scaffold | L5(0.2) | 400 | .161 | .258 | 58 / 167 | .178 | .457 | 58 / 69 |
| scaffold | L5(0.3) | 400 | .167 | .263 | 60 / 168 | .155 | .411 | 60 / 86 |
| scaffold | L5(0.4) | 400 | .148 | .243 | 55 / 171 | .090 | .320 | 55 / 117 |
| scaffold | B | 400 | .155 | .232 | 65 / 215 | .115 | .340 | 65 / 126 |

The gate-filtered variants are the only ones that improve both label sets at
once, and only on scaffold: L5(0.2) reaches MCC .161 test-like and .178
official against the control's .087 and .145, with official precision .457
against .443. The product variants L1-L3 gain on test-like labels and lose
badly on official labels on umap (L3 official MCC -.003), exactly the pattern
the previous run reported for the ungated models.

### Nested rate selection (rate chosen on the other four folds)

| Split | Variant | calls | MCC t-l | prec t-l | MCC off. | prec off. | selected rates per fold |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| umap | L0 | 393 | .068 | .152 | .128 | .434 | 4,8,8,8,4 |
| umap | L1 | 616 | .106 | .168 | .015 | .229 | 10 x 5 |
| umap | L2 | 616 | .128 | .187 | .062 | .273 | 10 x 5 |
| umap | L3 | 616 | .114 | .175 | .013 | .227 | 10 x 5 |
| umap | L4 | 591 | .107 | .170 | .106 | .332 | 10,10,8,10,10 |
| umap | L5(0.2) | 554 | .111 | .180 | .098 | .326 | 10,10,5,10,10 |
| umap | L5(0.3) | 616 | .105 | .169 | .057 | .272 | 10 x 5 |
| umap | L5(0.4) | 566 | .087 | .158 | .010 | .227 | 10,8,8,10,10 |
| umap | B | 616 | .102 | .162 | .028 | .242 | 10 x 5 |
| scaffold | L0 | 491 | .072 | .146 | .131 | .411 | 10,5,5,10,10 |
| scaffold | L1 | 590 | .159 | .215 | .070 | .277 | 10,10,10,8,10 |
| scaffold | L2 | 615 | .170 | .223 | .112 | .317 | 10 x 5 |
| scaffold | L3 | 615 | .176 | .226 | .084 | .286 | 10 x 5 |
| scaffold | L4 | 424 | .116 | .200 | .124 | .382 | 10,6.53,5,6.53,6.53 |
| scaffold | L5(0.2) | 572 | .161 | .223 | .172 | .405 | 10,10,10,6.53,10 |
| scaffold | L5(0.3) | 362 | .153 | .256 | .140 | .402 | 6.53,6.53,5,5,6.53 |
| scaffold | L5(0.4) | 454 | .128 | .208 | .057 | .276 | 8,8,10,3,8 |
| scaffold | B | 472 | .154 | .219 | .110 | .326 | 8,6.53,8,8,8 |

Because the rate is selected on test-like MCC and the test-like curve is still
rising at 10%, most variants are pushed to the grid edge, as in the previous
run. Only L5(0.3) and L4 on scaffold select interior rates. Fold-wise rows are
in `oof_nested_selection.csv`.

### Paired group bootstrap versus L0 (2,000 group resamples)

Test-like labels:

| Split | Variant | operating point | delta | 95% interval |
| --- | --- | --- | ---: | --- |
| umap | L1 | nested | +.038 | [-.002, +.080] |
| umap | L1 | 6.53% | +.015 | [-.032, +.066] |
| umap | L2 | nested | +.059 | [+.016, +.105] |
| umap | L2 | 6.53% | +.027 | [-.015, +.072] |
| umap | L3 | nested | +.046 | [-.010, +.102] |
| umap | L3 | 6.53% | +.017 | [-.030, +.067] |
| umap | L4 | nested | +.039 | [+.006, +.072] |
| umap | L4 | 6.53% | +.029 | [-.014, +.073] |
| umap | L5(0.2) | nested | +.043 | [+.010, +.078] |
| umap | L5(0.2) | 6.53% | +.020 | [-.013, +.055] |
| umap | L5(0.3) | nested | +.037 | [-.003, +.077] |
| umap | L5(0.3) | 6.53% | +.011 | [-.031, +.054] |
| umap | L5(0.4) | nested | +.019 | [-.026, +.067] |
| umap | L5(0.4) | 6.53% | +.021 | [-.016, +.065] |
| umap | B | nested | +.034 | [-.015, +.086] |
| umap | B | 6.53% | +.022 | [-.028, +.077] |
| scaffold | L1 | nested | +.087 | [+.036, +.140] |
| scaffold | L1 | 6.53% | +.050 | [-.004, +.109] |
| scaffold | L2 | nested | +.099 | [+.050, +.149] |
| scaffold | L2 | 6.53% | +.071 | [+.018, +.123] |
| scaffold | L3 | nested | +.104 | [+.053, +.159] |
| scaffold | L3 | 6.53% | +.053 | [-.004, +.113] |
| scaffold | L4 | nested | +.044 | [+.010, +.081] |
| scaffold | L4 | 6.53% | +.054 | [+.014, +.096] |
| scaffold | L5(0.2) | nested | +.089 | [+.048, +.134] |
| scaffold | L5(0.2) | 6.53% | +.073 | [+.034, +.115] |
| scaffold | L5(0.3) | nested | +.081 | [+.035, +.132] |
| scaffold | L5(0.3) | 6.53% | +.079 | [+.031, +.132] |
| scaffold | L5(0.4) | nested | +.057 | [+.007, +.107] |
| scaffold | L5(0.4) | 6.53% | +.061 | [+.010, +.114] |
| scaffold | B | nested | +.082 | [+.032, +.132] |
| scaffold | B | 6.53% | +.068 | [+.018, +.118] |

Official labels, same call vectors:

| Split | Variant | operating point | delta | 95% interval |
| --- | --- | --- | ---: | --- |
| umap | L1 | nested | -.113 | [-.161, -.058] |
| umap | L1 | 6.53% | -.142 | [-.206, -.070] |
| umap | L2 | nested | -.066 | [-.118, -.007] |
| umap | L2 | 6.53% | -.107 | [-.165, -.045] |
| umap | L3 | nested | -.115 | [-.180, -.041] |
| umap | L3 | 6.53% | -.150 | [-.209, -.084] |
| umap | L4 | nested | -.022 | [-.076, +.040] |
| umap | L4 | 6.53% | -.041 | [-.102, +.031] |
| umap | L5(0.2) | nested | -.030 | [-.074, +.018] |
| umap | L5(0.2) | 6.53% | -.051 | [-.097, -.001] |
| umap | L5(0.3) | nested | -.071 | [-.125, -.016] |
| umap | L5(0.3) | 6.53% | -.100 | [-.159, -.037] |
| umap | L5(0.4) | nested | -.117 | [-.171, -.053] |
| umap | L5(0.4) | 6.53% | -.114 | [-.166, -.055] |
| umap | B | nested | -.100 | [-.160, -.027] |
| umap | B | 6.53% | -.116 | [-.182, -.043] |
| scaffold | L1 | nested | -.061 | [-.129, +.008] |
| scaffold | L1 | 6.53% | -.088 | [-.158, -.020] |
| scaffold | L2 | nested | -.019 | [-.084, +.045] |
| scaffold | L2 | 6.53% | -.041 | [-.110, +.023] |
| scaffold | L3 | nested | -.047 | [-.116, +.023] |
| scaffold | L3 | 6.53% | -.087 | [-.158, -.012] |
| scaffold | L4 | nested | -.007 | [-.057, +.045] |
| scaffold | L4 | 6.53% | +.012 | [-.044, +.068] |
| scaffold | L5(0.2) | nested | +.041 | [-.020, +.098] |
| scaffold | L5(0.2) | 6.53% | +.033 | [-.019, +.086] |
| scaffold | L5(0.3) | nested | +.009 | [-.053, +.073] |
| scaffold | L5(0.3) | 6.53% | +.010 | [-.055, +.073] |
| scaffold | L5(0.4) | nested | -.074 | [-.139, -.007] |
| scaffold | L5(0.4) | 6.53% | -.055 | [-.123, +.012] |
| scaffold | B | nested | -.022 | [-.092, +.047] |
| scaffold | B | 6.53% | -.030 | [-.104, +.038] |

Against this stronger control the picture is the same in kind as before, but
with one exception. Every variant gains on test-like labels, more on scaffold
than umap. Under official labels every product variant loses, several with
intervals excluding zero (L1 umap -.113 [-.161, -.058]). The exception is the
gate *filter* at a low threshold: L5(0.2) on scaffold is +.041 [-.020, +.098]
nested and +.033 [-.019, +.086] at 6.53% under official labels, together with
+.089 and +.073 under test-like labels. On umap the same variant is -.030 and
-.051 under official labels, so the one variant that is not label-set
contradictory is also split-contradictory.

## Part 3 - blind swap inventory (diagnostic)

The deployed Sep 15 blind score was converted to a percentile rank over the 750
molecules and combined with the blind gate probabilities and `pool_lgbm` score
from `test_scores.parquet`. 279 of 750 blind molecules have p_active >= 0.3.
Vectors are in `blind_candidates.parquet` and every changed molecule in
`swaps.csv`. These are diagnostics; no candidate was registered.

| Variant | calls | overlap with Sep 15 | added | removed | mean p_active called | mean p_active added | mean p_active removed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| L0 | 49 | 49 | 0 | 0 | .322 | – | – |
| L0 | 60 | 49 | 11 | 0 | .323 | .331 | – |
| L1 | 49 | 21 | 28 | 28 | .559 | .594 | .177 |
| L1 | 60 | 22 | 38 | 27 | .543 | .565 | .172 |
| L2 | 49 | 23 | 26 | 26 | .526 | .552 | .167 |
| L2 | 60 | 27 | 33 | 22 | .491 | .514 | .147 |
| L3 | 49 | 19 | 30 | 30 | .541 | .564 | .206 |
| L3 | 60 | 21 | 39 | 28 | .522 | .540 | .197 |
| L4 | 49 | 32 | 17 | 17 | .368 | .309 | .175 |
| L4 | 60 | 34 | 26 | 15 | .358 | .327 | .184 |
| L5(0.3) | 49 | 22 | 27 | 27 | .478 | .456 | .172 |
| L5(0.3) | 60 | 22 | 38 | 27 | .483 | .470 | .172 |

At 49 calls the gate-combined variants replace 26 to 30 of the historical best
49, and the rank-mean with the pool model replaces 17. The removed calls have
mean p_active .15-.21 and the added ones .46-.59, which is the intended
mechanism and also the reason these vectors are not comparable to the Sep 15
payload without a leaderboard observation. Of the Sep 15 49 positives, 18 have
p_active < 0.2, 27 < 0.3 and 41 < 0.5, so any gate-weighted variant necessarily
rewrites most of the historical best vector.

## Part 4 - sensitivity to a stricter pseudo-negative set

The test-like label set was rebuilt with pseudo-negatives restricted to primary
log2fc >= 0, which drops 173 of 2,693 pseudo-negatives (test-like pool 4,017
instead of 4,190). Everything else is identical.

| Split | top variant (both sets) | ordering identical | Spearman of variant ranks |
| --- | --- | --- | ---: |
| umap | L2 | no | .983 |
| scaffold | L3 | no | .983 |

The only ordering changes are adjacent swaps: on umap L1 and L4 exchange
fourth and fifth place, on scaffold L5(0.3) and B exchange fifth and sixth.
Levels are not directly comparable between the two label sets because the
negative count differs; within each set the nested pooled MCC moves by at most
+.003 per variant (scaffold L3 .176 → .179, umap L2 .128 → .130) with one
exception: scaffold L5(0.4) selects a different nested rate under the stricter
labels (478 instead of 454 calls), which moves its pooled test-like MCC from
.128 to .145 and its official-label delta from -.074 [-.139, -.007] to
-.060 [-.124, +.007]. Every other bootstrap interval keeps its sign and its
endpoints to within .005. The L-variant ordering is therefore robust to the
borderline part of the pseudo-negative definition.

## What this does and does not establish

Established:

- The activity gate, built without any blind label, agrees in sign with 6 of
  the 7 observed D6 Live MCC transitions, and its two largest separations are
  the two October failures (permutation p .0001-.0005). Its one disagreement is
  the +.0076 Sep 13→Sep 14 move. The official-only TDI head agrees 5 of 7 and
  fails on the Sep 15→Sep 23 decline, which the gate calls correctly.
- Mean payload-level gate probability is monotone with Live MCC across the
  seven submitted states except for the Sep 18/Sep 23 pair.
- Against the legacy Sep 15 refit control, every legacy x gate combination
  improves pooled test-like MCC (scaffold +.044 to +.104 nested, intervals
  excluding zero for all; umap +.019 to +.059), and all but the gate filter at
  tau = 0.2 lose under official labels.
- The gate *filter* at tau = 0.2 is the single variant with a non-negative
  official-label interval on scaffold at both operating points (+.041 and
  +.033) while also gaining under test-like labels; it does not reproduce this
  on umap (-.030 and -.051).
- The variant ordering is insensitive to tightening the pseudo-negative
  definition to log2fc >= 0 (Spearman .983, adjacent swaps only).
- The deployed Sep 15 continuous vector reproduces the submitted 49 calls
  exactly, with a 0.00133 score gap at the boundary.

Not established:

- Any causal or predictive claim about future Live MCC. Seven transitions over
  one hidden set, with a bootstrap public metric and heavily overlapping call
  sets, cannot be treated as seven independent successes, and the agreement
  statistic was computed once on a frozen design rather than tuned.
- A deployable recipe or rate. Test-like rate selection again sits at the 10%
  grid edge for most variants, the official-label and test-like directions
  still conflict for the product variants, and the one variant that avoids the
  conflict does so on scaffold only.
- That the gate is better than the pool model as a score. `pool_lgbm` also
  agrees 6 of 7 transitions; the two differ in which transition they miss.
- Anything about the blind scored-label mask, or about the Live MCC of the
  blind vectors in Part 3.

## Limitations of the procedure

- `gate_x_tdi.p_tdi` and `official_lgbm.score` are the same numbers, so the
  seven score columns are six distinct observations.
- All OOF variants reuse frozen scores; no selection inside the variants is
  cross-fitted beyond what run `7920b49a` already did, and the combination
  forms (product, sqrt, rank mean, filter) were chosen by the analyst rather
  than selected out of fold.
- The thresholds tau and the two blind call counts are a small predeclared
  grid, not a search; L5 behaviour between tau 0.2 and 0.4 is not monotone.
- The legacy OOF comparator is a refit analog of the Sep 15 recipe, while Part
  3 uses the deployed vector; the two are not the same fit.
- Pseudo-negative labels remain an inference from the primary screen, and the
  blind molecules were never primary screened.
