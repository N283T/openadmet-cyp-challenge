# D6 test-like pool with primary-screen pseudo-negatives and an activity gate

Recorded run `20261002T215345688167Z-7920b49a` (classification, 467 s).
Code: `scripts/evaluate_cyp_d6_testlike_pool.py`,
`scripts/launch_cyp_d6_testlike_pool.py`,
`configs/classification/d6-testlike-pool-v1.json`.
Compact outputs with SHA-256 provenance:
`artifacts/tracked/cyp-d6-testlike-pool-20261002T215345688167Z-7920b49a/`.
No candidate, prepared submission or upload resulted from this run.

## Question

All 1,497 supplied CYP2D6_is_TDI labels belong to primary-screen hits (TDI arm,
50 uM, log2fc < -0.5). A label-positive needs TDI-arm pIC50 > 4.301, which
implies roughly half inhibition at 50 uM (log2fc about -1 or lower), so screened
non-hits would almost certainly have been assigned negative under the official
rule. The 750 blind molecules were measured in dose response regardless of
primary screening. This run asks whether evaluating and training on
official + screened-non-hit labels changes the model ranking and the selected
call rate, and whether an explicit activity gate helps.

## Pool and counts

D6 log2fc is absent from `measurements.parquet`, so
`data/raw/cyp-challenge-single-concentration-TRAIN.csv` is declared as an input
and joined through `aliases.parquet` (17,504 rows, 4,376 molecules per enzyme,
all at 50 uM).

| Role (development, n = 6,146) | Rule | Molecules | Use |
| --- | --- | ---: | --- |
| `official` | supplied D6 label | 1,497 (324 positive) | fitted and scored |
| `pseudo_negative` | unlabelled, log2fc >= -0.25 | 2,693 | label False, fitted and scored |
| `borderline` | unlabelled, -0.5 <= log2fc < -0.25 | 130 | excluded; calls counted |
| `unlabelled_hit` | unlabelled, log2fc < -0.5 | 58 | excluded; calls counted |
| `unscreened` | unlabelled, no log2fc | 1,768 | never scored; callable |

Test-like pool = 4,190 molecules, positive rate 7.73%. Verification: 1,495 of
1,497 labelled molecules carry a D6 log2fc and the largest labelled log2fc is
-0.592, so labels really are confined to primary hits. The converse fails for
58 hits that have no label. 66 of the 324 positives sit at log2fc >= -1, so the
gate threshold of -1 discounts about a fifth of known positives; 1,198 of the
4,376 screened development molecules satisfy log2fc < -1.

## Design

Outer splits are the fixed five-fold `umap` and `scaffold` assignments
(`20260909T135130268493Z-edba6fc0`); every development molecule of a held fold
receives a score. Inner selection uses GroupKFold(3) on the same family's
groups inside the training folds only. Features are fresh Morgan r2/2048 plus
217 RDKit 2D descriptors (2,265 columns) and frozen CheMeleon 2,048-d
(`20260910T012630927160Z-e9576045`). LightGBM is fixed at 400 trees, 31 leaves,
learning rate 0.03, colsample 0.5, min_child_samples 20, balanced class weights,
two seeds averaged; no hyperparameter search was run.

| Model | Recipe |
| --- | --- |
| A `official_lgbm` | control: Morgan+RDKit LightGBM, official labels only |
| B `pool_lgbm` | same model on official + pseudo-negative labels |
| C `gate_x_tdi` | p_active x p_tdi; p_tdi from A |
| D `pool_chemeleon` | balanced logistic regression on CheMeleon, pool labels |
| E `gate_x_tdi_chemeleon` | as C with CheMeleon for gate and TDI head |
| F `rank_mean_gates` | equal percentile-rank mean of C and E |
| `legacy_sep15_refit` | comparator, see below |

The gate is a LightGBM regressor for D6 log2fc fitted on all screened training
molecules, then Platt-calibrated to P(log2fc < -1) on inner group-fold OOF
predictions of that regressor. Inner Spearman of predicted against observed
log2fc is .48-.54 across folds. The CheMeleon logistic penalty C = 0.01 was
selected by inner AP in all eleven fits.

Calls: in each held fold the top round(r x n_fold) molecules of the complete
development pool are called positive, for r in {3, 4, 5, 6, 7, 8, 10}% plus the
historical analog 49/750 = 6.53%. The nested rate is the grid value with the
best pooled test-like MCC on the other four folds. Uncertainty is a paired
group bootstrap (1,000 resamples of assignment groups) of pooled
MCC(model) - MCC(A).

Comparator: the complete-pool nested OOF of the refitted Sep 15 recipe
(`20261002T120413515280Z-f2330789`), averaged over its two seeds. Its coverage,
fold assignment and labels match this run exactly (0 missing, 0 fold mismatches,
0 label mismatches over 12,292 rows), so it is scored on the same pools under
the same rate procedure. It contributes no blind scores here.

## Ranking metrics (pooled held-fold scores)

| Split | Model | AUROC off. | AP off. | AUROC test-like | AP test-like |
| --- | --- | ---: | ---: | ---: | ---: |
| umap | A official_lgbm | .534 | .251 | .509 | .085 |
| umap | B pool_lgbm | .528 | .239 | .639 | .133 |
| umap | C gate_x_tdi | .516 | .228 | .601 | .116 |
| umap | D pool_chemeleon | .527 | .230 | .611 | .113 |
| umap | E gate_x_tdi_chemeleon | .532 | .239 | .619 | .120 |
| umap | F rank_mean_gates | .534 | .239 | .632 | .126 |
| umap | legacy_sep15_refit | .586 | .301 | .542 | .101 |
| scaffold | A official_lgbm | .584 | .280 | .536 | .095 |
| scaffold | B pool_lgbm | .579 | .282 | .685 | .167 |
| scaffold | C gate_x_tdi | .566 | .262 | .644 | .144 |
| scaffold | D pool_chemeleon | .569 | .259 | .648 | .132 |
| scaffold | E gate_x_tdi_chemeleon | .572 | .259 | .653 | .139 |
| scaffold | F rank_mean_gates | .583 | .277 | .672 | .159 |
| scaffold | legacy_sep15_refit | .638 | .348 | .587 | .123 |

The two label sets rank the models differently. On official labels alone the
legacy refit leads both families (AP .301 / .348) and no new model beats the
control by more than .002 AP. On test-like labels the ordering inverts: every
pool-trained or gated model gains .02-.07 AP over the control, the legacy refit
falls behind all of them, and `pool_lgbm` leads (AP .133 umap, .167 scaffold).

## Decision metrics

Nested rate selection (rate chosen on the other four folds, MCC on the held
fold, pooled over folds) and the fixed 6.53% analog:

| Split | Model | nested MCC test-like | nested MCC off. | nested calls | 6.53% MCC test-like | 6.53% MCC off. |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| umap | A official_lgbm | .042 | .093 | 616 | .034 | .062 |
| umap | B pool_lgbm | .102 | .028 | 616 | .096 | .031 |
| umap | C gate_x_tdi | .078 | .012 | 579 | .081 | .026 |
| umap | D pool_chemeleon | .076 | .025 | 491 | .050 | -.004 |
| umap | E gate_x_tdi_chemeleon | .068 | .014 | 554 | .076 | .031 |
| umap | F rank_mean_gates | .085 | .017 | 517 | .089 | .033 |
| umap | legacy_sep15_refit | .068 | .128 | 393 | .074 | .147 |
| scaffold | A official_lgbm | .053 | .099 | 257 | .043 | .081 |
| scaffold | B pool_lgbm | .155 | .110 | 454 | .155 | .115 |
| scaffold | C gate_x_tdi | .112 | .066 | 578 | .119 | .071 |
| scaffold | D pool_chemeleon | .092 | .053 | 468 | .094 | .051 |
| scaffold | E gate_x_tdi_chemeleon | .123 | .075 | 615 | .113 | .062 |
| scaffold | F rank_mean_gates | .123 | .073 | 590 | .121 | .074 |
| scaffold | legacy_sep15_refit | .075 | .136 | 442 | .087 | .145 |

Nested-selected rates are mostly at the top of the grid: 10% in all five umap
folds for A, B and C, and 7-10% for the gated models on scaffold. A on scaffold
is the exception (3-5%). The pooled test-like rate curve is flat to slightly
increasing out to 10% for every model, so the grid boundary, not an interior
optimum, sets the rate. `pool_lgbm` on scaffold peaks at 7-8% (.162) and is
.155 at 6.53%.

At 6.53% the call composition differs between models. The control places 88
(umap) / 89 (scaffold) of its ~400 calls on officially labelled molecules,
whereas the pool and gated models place 158-191 there and correspondingly fewer
on pseudo-negatives. Test-like precision at that rate is .114 for the control
versus .177 (`pool_lgbm`, umap) and .232 (`pool_lgbm`, scaffold). Calls landing
on excluded molecules stay small: 3-11 borderline and 99-145 unscreened per
model and split.

## Paired group bootstrap versus the control

Test-like labels, 1,000 group resamples, pooled over held folds:

| Split | Model | operating point | delta | 95% interval |
| --- | --- | --- | ---: | --- |
| umap | B pool_lgbm | nested | +.060 | [+.024, +.099] |
| umap | B pool_lgbm | 6.53% | +.061 | [+.018, +.109] |
| umap | C gate_x_tdi | nested | +.036 | [+.010, +.066] |
| umap | C gate_x_tdi | 6.53% | +.047 | [+.014, +.083] |
| umap | D pool_chemeleon | nested | +.034 | [+.001, +.069] |
| umap | D pool_chemeleon | 6.53% | +.015 | [-.022, +.060] |
| umap | E gate_x_tdi_chemeleon | nested | +.026 | [-.027, +.083] |
| umap | E gate_x_tdi_chemeleon | 6.53% | +.042 | [-.019, +.111] |
| umap | F rank_mean_gates | nested | +.043 | [+.011, +.079] |
| umap | F rank_mean_gates | 6.53% | +.054 | [+.010, +.108] |
| umap | legacy_sep15_refit | nested | +.026 | [-.018, +.074] |
| scaffold | B pool_lgbm | nested | +.102 | [+.054, +.156] |
| scaffold | B pool_lgbm | 6.53% | +.113 | [+.068, +.163] |
| scaffold | C gate_x_tdi | nested | +.059 | [+.018, +.101] |
| scaffold | C gate_x_tdi | 6.53% | +.076 | [+.039, +.121] |
| scaffold | D pool_chemeleon | nested | +.039 | [-.009, +.091] |
| scaffold | D pool_chemeleon | 6.53% | +.051 | [+.001, +.104] |
| scaffold | E gate_x_tdi_chemeleon | nested | +.070 | [+.020, +.119] |
| scaffold | E gate_x_tdi_chemeleon | 6.53% | +.070 | [+.020, +.119] |
| scaffold | F rank_mean_gates | nested | +.070 | [+.026, +.117] |
| scaffold | F rank_mean_gates | 6.53% | +.078 | [+.035, +.127] |
| scaffold | legacy_sep15_refit | nested | +.022 | [-.031, +.104] |

The same bootstrap on official labels only reverses most of these signs:
`pool_lgbm` is -.065 [-.112, -.010] on umap and +.011 [-.056, +.079] on
scaffold, `gate_x_tdi` is -.081 [-.132, -.035] on umap, and the legacy refit is
the only model with a positive official-label interval at 6.53%
(+.085 [+.005, +.168] umap, +.064 [+.000, +.130] scaffold). The direction of
every comparison therefore depends on which label set is assumed correct, and
the pseudo-negative assumption is exactly what is at stake.

## Blind test diagnostics

Each recipe was refit on all 6,146 development molecules (gate and calibration
included, inner groups taken from the umap family) and applied to the 750 blind
molecules. Keyed scores are in `test_scores.parquet`
(molecule_name, compound_id, model, score, p_active, p_tdi). The reference is
the Sep 15 payload `20260924T050705396203Z-10e1d3ad` with 49 positives
(Live MCC .3086). The rate column below is each model's pooled-OOF optimum for
that family; the fold-wise nested rates are in `blind_summary.csv` and
`nested_selection.csv`.

| Model | p_active > 0.5 | top-49 overlap | umap rate / calls / overlap | scaffold rate / calls / overlap |
| --- | ---: | ---: | --- | --- |
| A official_lgbm | n/a | 31 | 10% / 75 / 36 | 5% / 38 / 25 |
| B pool_lgbm | n/a | 22 | 10% / 75 / 25 | 8% / 60 / 23 |
| C gate_x_tdi | 137 | 27 | 10% / 75 / 33 | 10% / 75 / 33 |
| D pool_chemeleon | n/a | 17 | 8% / 60 / 20 | 6% / 45 / 16 |
| E gate_x_tdi_chemeleon | 144 | 24 | 10% / 75 / 30 | 10% / 75 / 30 |
| F rank_mean_gates | n/a | 27 | 7% / 53 / 30 | 10% / 75 / 34 |

The gate predicts primary activity for 137 (Morgan) and 144 (CheMeleon) of the
750 blind molecules, mean p_active .31 and .30; the development base rate of
log2fc < -1 among screened molecules is 1,198/4,376 = 27.4%, so the gate does
not single out the blind set as unusually active.

Top-49 overlap with the Sep 15 positives is 17-31 of 49. The control overlaps
most (31) and the two pool-trained models least (22 and 17), which is the same
direction as the official-label metrics and the opposite of the test-like ones.
Transferring each family's held-fold score threshold instead of its rank gives
quite different blind counts (for instance 52 and 49 calls for `gate_x_tdi`
against 119 and 106 for `gate_x_tdi_chemeleon`), so score-level calibration does
not transfer consistently across feature views.

## What this does and does not establish

Established by this run:

- The label-coverage structure is exactly as assumed: all 1,497 labels lie
  inside the primary hits (largest labelled log2fc -0.592), and the
  official + pseudo-negative pool has a 7.73% positive rate, close to the
  7.06% estimate for the blind set.
- Under test-like labels, training on the pool or multiplying by a calibrated
  activity gate beats the official-labels-only control on both split families,
  by +.06 to +.11 pooled MCC with intervals excluding zero at both operating
  points for `pool_lgbm`, and by +.04 to +.08 for the gated and rank-mean
  variants. The plain pool-trained LightGBM, not the gate, is the strongest
  model, so the gate is not the ingredient that produces the gain.
- Model ordering is label-set dependent. The legacy Sep 15 refit remains the
  best model on official labels (AP .301 / .348, MCC .147 / .145 at 6.53%) and
  the worst or near-worst among the new models on test-like labels. The sign of
  most control comparisons flips between the two label sets.
- The test-like criterion pushes the selected call rate to the top of the
  evaluated grid (mostly 7-10%, i.e. 53-75 blind calls), for every model and
  both families.

Not established:

- That the pseudo-negatives are correct. They are an inference from the primary
  screen, not measurements. 66 known positives have log2fc >= -1 and 58 primary
  hits were never labelled, so neither the gate threshold nor the -0.25 cut is
  a clean separation. Any gain measured against pseudo-negatives partly
  rewards predicting primary inactivity, which the official blind labels do not
  reward directly.
- That this transfers to the leaderboard. The blind molecules were not primary
  screened, so the "test-like" pool is an analogy about label prevalence, not
  about acquisition. The test-like gain and the official-label loss are
  consistent with the known failure mode in which complete-pool diagnostics
  missed the Live direction.
- Anything about the blind scored-label mask or about which of the 750
  molecules would be scored.
- A deployable rate. Selection was performed on umap and scaffold only, the
  curves are flat, and the chosen rate sits at the grid boundary; a wider grid
  would likely move it.

Consequently this run is evidence that the historical labelled-only evaluation
and the pseudo-negative evaluation disagree systematically, and that the
disagreement is large enough to reverse model choice. It is not evidence that
any model here would improve the D6 Live MCC. No upload is recommended or
implied.

## Limitations of the procedure itself

- Pseudo-negative labels enter both training and scoring, so a model that
  learns primary-screen activity is rewarded twice.
- `borderline` and `unlabelled_hit` molecules are excluded from scoring but
  remain callable; their call counts (3-11 and small) are reported per rate in
  `rate_curves.csv`.
- The legacy comparator is an OOF produced by another run under the same fold
  assignments and seeds, not a refit inside this script; it has no blind
  predictions here.
- Two LightGBM seeds are averaged; a single CheMeleon logistic fit per split
  uses no seed averaging because the solver is deterministic.
