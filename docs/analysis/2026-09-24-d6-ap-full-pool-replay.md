# D6 AP-mass90 full-pool replay and remaining LB mismatch

## Question and frozen comparison

The September 23 AP-mass90 submission lost D6 leaderboard MCC against the
September 15 control (`.23483` versus `.30864`) despite favorable development
OOF. This replay asks whether the OOF decision changes when each fold ranks all
held-out molecules, including those without a D6 TDI label, and when AP weights
and the positive-call rate are chosen without the evaluated outer fold.

The comparison here uses the AP ensemble's leading ExtraTrees atom as its
control. It is **not** the September 15 submitted mean-block control. The
historical LB contrast remains the independent external observation.

## Reproduced predictions

Runs `c42e6f27`, `48b841b3`, `eaac6dd8`, `62dde800`, `95f78765` and
`ec9287a1` produced the seven AP member predictions on all 6,146 UMAP-fold
training molecules. The first six members reproduce the original labelled OOF
predictions exactly, except the TabPFN 3.5 member, whose full-batch inference
changes labelled probabilities by at most `.00162` without changing labelled
fixed-rate calls. The legacy TabPFN 2.6 refit has per-fold Spearman correlation
at least `.9986` to its old labelled OOF, with mean absolute probability
difference `.00130`. Its saved fit states could not be restored in the current
TabPFN version because `date_transformer_` is missing; the failed restoration
run `3dc862a3` is retained, and `ec9287a1` refits the frozen fold recipe.

Run `b719ba35` ranks each member within the complete held-fold pool and
combines all seven with either the submitted full-fit weights (descriptive) or
the saved outer-excluded AP weights. The source ensemble, exact member names,
weights and fold assignments are retained in its immutable recorded inputs.
The 1,497 labelled D6 molecules contribute to MCC; the other 4,649 molecules
consume call budget but contribute no class label.

## Operating-point result

At the historical 49/750 call fraction and estimated 7.06% scored-population
prior, the full-fit seven-member blend appears `+.01035` MCC better than
ExtraTrees if calls are chosen among labelled molecules only. Ranking the
complete held-fold pool cuts that estimate to `+.00759` (paired 95% bootstrap
interval `[-.03532,+.05052]`). With outer-excluded weights the labelled-only
gain is `+.00514`, but the full-pool delta is exactly `0.00000` (interval
`[-.04159,+.03802]`). None supports a replacement.

Selecting each model's call rate on other outer folds gives the full-fit blend
`−.00493` MCC relative to ExtraTrees on the full pool, versus `+.02677` in the
labelled-only diagnostic. The outer-excluded seven-member blend gives
`−.02695` (interval `[-.06241,+.00446]`) on the full pool. The selected rates
are unstable and reach the grid boundaries (15% for ExtraTrees, 20% for the
outer-excluded blend). They do not justify applying those rates to the blind
test. A six-member replay omitting the 2.7%-weight legacy TabPFN is more
negative (`−.06442`, interval `[-.11936,−.00317]`) under outer-excluded UMAP
rate selection; that stronger estimate must not be substituted for the exact
seven-member result.

These results identify two concrete sources of OOF optimism in comparisons
against AP's leading atom: scoring a labelled-only call budget instead of the
complete molecular pool, and treating an AP-optimized full-fit ensemble and
rate as though they were chosen before the evaluated fold. They do not prove
that either source caused the public LB loss.

## Direct September 15 control check

The September 15 recipe was also reconstructed on all UMAP held-fold molecules.
Frozen LightGBM and graph fits reproduce their original labelled OOF scores
exactly; the two TabPFN 2.6 variants were refit with their saved preprocessing,
and the fixed Morgan tail fits were restored. Run `914bb06c` composes both a
raw-score analog of the submitted blind recipe and a rank-before-core analog of
the historical labelled OOF. These analogs have Spearman correlations `.9874`
and `.9781`, respectively, to the saved deployed labelled OOF. They are not an
exact re-run of its original cross-fitted recipe selection.

Directly comparing the seven-member AP replay with the raw-score September 15
analog at 49/750 and 7.06% prior gives **AP `+.00711` MCC** with outer-excluded
weights on the complete pool (95% group-bootstrap interval
`[-.03670,+.04621]`; run `638c65b7`). The full-fit weights give `+.01470`
(`[-.02026,+.04692]`; run `f3e60b06`). The rank-before-core control sensitivity
still gives AP `+.01278` (`[-.02796,+.05086]`; run `ae37ff29`). Thus the direct
UMAP replay **does not predict the actual LB direction** (`.23483` AP versus
`.30864` September 15). This supersedes any interpretation that full-pool
ranking alone resolved the OOF-to-LB mismatch. The small local deltas and broad
intervals support no AP promotion, but the source of the public reversal remains
unidentified. A common alternative outer split, with both recipes genuinely
refit and reselected, is still required. The Sep 18 → Sep 23 transition remains
reserved as a later check of a frozen protocol.

## Alternative split and blind-domain check

Both recipes were refit on the saved scaffold folds. The September 15 atoms
were refit as LightGBM (`45b73abd`), Morgan tails (`0f606ad8`), PCA/raw
TabPFN 2.6 (`b92b2044`, `3f9d3bb6`) and TDI-only graph with inner epoch
selection (`e3a605c0`). Run `e17778e3` composes the fixed submitted control
recipe. All seven AP atoms are available on the same scaffold folds;
`8db19b4b` freezes the seven submitted members and chooses their AP weights
using only other scaffold folds. The member set itself originated from UMAP
labels, so this is a robustness check rather than a fully independent recipe
selection.

At 49/750 and 7.06% prior, scaffold OOF favors AP by `+.04548` MCC over the
raw-score control analog (95% interval `[+.00137,+.09008]`, run `99cb9161`).
When the call rate is selected on other folds, the difference is `−.00130`
(`[-.04455,+.04149]`), and AP's chosen rates range from 2% to 20%. A
rank-before-core control sensitivity still favors AP at 49/750 (`+.05874`,
run `99139d5f`). Thus scaffold fixed-rate OOF also misses the public reversal;
rate selection reduces the apparent gain but is too unstable to set a blind
rate.

Run `4135649a` uses a five-fold cross-fitted molecular domain classifier to
estimate the train-to-blind density ratio from Morgan fingerprints and eight
descriptors, without hidden labels. Its domain AUROC is `.877`; clipping weights
to `[.2,5]` lowers the labelled effective sample from 1,497 to about 782.
Density weighting increases the AP-versus-control fixed-rate delta from
`+.00711` to `+.04511` on UMAP and from `+.04548` to `+.08766` on scaffold.
This prespecified covariate weighting worsens concordance with LB and is
rejected as a transfer correction.

The recorded conservative evidence gate `49bb67b7` requires the paired lower
interval to exceed zero for both the historical 49/750 rate and outer-fold
rate selection on UMAP and scaffold. AP clears only the scaffold fixed-rate
check (one of four); its status is `INSUFFICIENT_EVIDENCE`. This gate is a
minimum standard for considering another blind probe, not a claim that a
passing candidate will improve LB. It makes the observed uncertainty and
split dependence part of the candidate decision rather than treating a single
OOF point estimate as sufficient.

## Leading atom and label-route checks

ExtraTrees alone was compared directly with the September 15 control in run
`25c958cf`. At the 49/750 rate it is only `+.00711` MCC on UMAP (95% interval
`[-.04486,+.06242]`) and `+.03159` on scaffold
(`[-.01938,+.08371]`). Outer-fold rate selection favors ExtraTrees on UMAP
(`+.06802`, interval `[+.00520,+.12865]`) but favors the control on scaffold
(`−.02095`, interval `[−.06910,+.02651]`). Its selected rate is 15% in every
UMAP fold and ranges from 6% to 20% on scaffold. De-ensembling AP therefore
does not give a stable replacement at the submitted boundary.

Run `c55d672a` joins fixed-rate complete-pool decisions to training assay
routes. AP exchanges 211 calls each way on UMAP, with 17 supplied positives
added and 15 removed; scaffold exchanges 178 each way, with 20 added and 13
removed. On UMAP, AP's net two positives comprise +3 in the `overlap_zero`
shift-bound stratum and −1 outside it. On scaffold, five of the net seven are
in that stratum. These large exchanges and small labelled gains make the OOF
ranking fragile. `overlap_zero` is a training interval descriptor, **not** a
verified replica of the hidden scoring mask; the blind compounds' label routes
are unknown. Run `622cb62b` additionally evaluates a hypothetical mask that
keeps only training labels whose shift interval excludes the twofold threshold.
It leaves 485 scored rows and 105 positives. AP still leads at 49/750: `+.02725`
on UMAP (95% interval `[−.07179,+.12143]`) and `+.09523` on scaffold
(`[−.00318,+.20267]`). That sensitivity therefore **does not recover** the
public LB direction. It also cannot establish the actual hidden scoring mask.

## Conservative control-anchored blend

A prespecified 75% September 15 plus 25% AP rank blend was composed on both
complete-pool splits in run `7352d993` and evaluated in `2fe73ea9`. At 49/750
and 7.06% prior it changes projected MCC by `−.00136` on UMAP (95% interval
`[−.02425,+.01741]`) and `+.02918` on scaffold
(`[+.00383,+.05655]`). With outer-fold rate selection the deltas are `+.01382`
(`[−.03011,+.04936]`) and `+.01259` (`[−.01159,+.03600]`), respectively.
The selected blend rates vary by fold. Its point estimates do not establish a
stable improvement at the historical boundary. For an explicitly exploratory
leaderboard probe, run `85d10ccf` nevertheless prepared a valid 750-row
prediction package with 49 D6 calls, four additions and four removals relative
to September 15, and A4 unchanged. Candidate
`20260924T043703823086Z-51106137` is frozen but **not uploaded**. Run
`d41baf66` records the eight changed compound IDs. The additions' mean
similarity-weighted positive-neighbor rate is `.315` versus `.310` for
removals, giving no independent support for an expected improvement. This
candidate fails the dual-split evidence gate; any upload would be a deliberate
small probe with limited downside in call changes, not a validated upgrade.

## Applied decision

Keep the September 15 submitted D6 vector as the leaderboard-validated
control. Do not promote AP-mass90 or a new D6 package on this replay. The
evaluation API now rejects an incomplete molecular pool by default; a
labelled-only calculation requires an explicit diagnostic opt-in. Future
candidate comparisons must use a complete-pool fold-excluded bank, a stated
prior and rate grid, and outer-excluded recipe selection before deployment.

Compact predictions and comparison tables are retained in
`artifacts/tracked/d6-ap-full-pool-replay-20260924/`. The full run IDs are
`20260924T035023873825Z-b719ba35` (composition),
`20260924T035058376907Z-d54f4758` (full-fit evaluation) and
`20260924T035102011139Z-81c91e3b` (outer-excluded evaluation).
