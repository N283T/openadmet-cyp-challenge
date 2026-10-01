# Boltz trunk regression adoption and prospective submission

The user adopted the improving CYP1A2 and CYP2D6 ensemble families and explicitly
authorized a regression upload. CYP2C9 and CYP3A4 return to their best scored
submitted vectors. Classification is unchanged.

Package run `regression/20261001T122958256323Z-1cf63fe0` completed without model
refitting. Candidate `20261001T123003757010Z-e61d651c` is prepared as submission
`20261001T123022061267Z-ad1689ae`, SHA-256
`0cc78087990bdd27ef1bdfdedc135a05d6f79dd8dc02208c14a183d977734559`.
The official validator and sender dry run both pass for 750 exact identities.

## Deployment composition

CYP1A2 selects its replacement cap using held-fold prediction error on all five
development folds, then fits the chosen greedy replacement procedure on all
development OOF. The resulting cap is two. Within the original ten-member equal
ensemble, `tabpfn35_fast_anchor59_raw_scores_logfc` is replaced by
`ligand_only--trunk--tabpfn_top500` and `jeremy_frozen_tabpfn26` by
`core_bidirectional_hem--trunk--tabpfn_top500`. Each replacement has weight 0.1;
the other eight original members retain weight 0.1. Both trunk recipes exclude
molecular feature augmentation.

CYP2D6 selects the best standalone trunk on all development OOF:
`pxr_core--trunk_molecular--lgbm_l1`, which combines protein and ligand single
representations and pocket-to-ligand pair pooling with CheMeleon/Mordred features.
Caruana selection on the original two top-level branches plus this trunk gives
weights 0.4 for the submitted two-member core, 0.4 for the moment-matched robust
branch and 0.2 for the trunk. Its step count is chosen using development folds.
This preserves the core/robust hierarchy and uses the previous single-trunk
extension family, which marginally outperformed the multiple-trunk bank.

Original members use retained full-fit blind predictions. Each selected trunk
uses the mean of its five saved fold-model blind predictions. D6 robust predictions
are moment-matched to the unplaced original core before composition. Each changed
endpoint's final prediction mean and SD are matched to the best scored control
`20260926T090548323924Z-46188564`, without another calibration search. The
unmodified blind controls reconstruct that submission within `3.55e-15` for 1A2
and `1.78e-15` for D6, validating retained atomic sources and placement.

Development ensemble selection is conditional on saved OOF, not fully nested
base-model refitting. The submitted composition is fitted on all development
OOF rather than treating the five outer evaluation compositions as a fixed
production recipe. The small reported gains have uncertainty intervals spanning
zero and do not establish improvement on the blind leaderboard.

## Changes versus the last sent payload

The last sent regression payload was the September 27 affine calibration probe,
`20260926T140601457029Z-bd3a1e40`, whose changes worsened Live ST-RAE.

| Target | Changed molecules | Final mean | Sample SD | Action |
|---|---:|---:|---:|---|
| 1A2 | 750 | 4.412 | 1.169409 | Two trunk replacements; restore best moments |
| 2C9 | 0 | 4.830 | 0.922638 | Exact best submitted vector |
| 2D6 | 750 | 3.570 | 0.883477 | Trunk extension; restore best moments |
| 3A4 | 750 | 4.880 | 1.082472 | Exact best submitted vector; undo affine probe |

C9 and A4 are exactly equal as parsed prediction columns to the best control;
their original source vectors were submitted September 17. The previous
cooldown ended September 28 at 02:13 JST. Before the live call, the prepared ID,
changes, model composition and expired cooldown were shown in the chat. The
user's explicit authorization covered this two-target adoption and rollback.

## Delivery

Attempt `20261001T123118187956Z-b567481f` captured all five pre-upload official
leaderboards and made one submission API call. It returned `sent`; the service
receipt says the Regression Prediction submission was received and is being
processed, with leaderboard appearance expected within two hours. Receipt is
delivery evidence, not a scored result. No retry was made. The next local
regression opportunity is **2026-10-02 09:30:42 JST**.

[Adopted configuration](../../configs/regression/current-boltz-adoption.json)
and [retained package outputs](../../artifacts/tracked/cyp-boltz-regression-adoption-20261001T122958256323Z-1cf63fe0/summary.json)
record exact composition and provenance. The previous best scored submission
remains the performance control until the new payload is scored.

## Reflected Live result

All five official Live boards now show N283T's **2026-10-01 12:31 UTC**
submission. Audit run `regression/20261001T125026856964Z-71f3f209` validates
the receipt, payload hash, 750 identities, exact C9/A4 control vectors and
matched submission timestamps. The unchanged C9 vector has identical values
for every leaderboard metric across all three snapshots, supporting evaluation
population continuity. A4 exactly reproduces every best-control metric.

| Target | Last sent | Previous best | New | Change vs previous best |
|---|---:|---:|---:|---:|
| 1A2 | .512345 | .508530 | **.506911** | **-.001619** |
| 2C9 | .359358 | .359358 | .359358 | .000000 |
| 2D6 | .549788 | .549341 | **.535721** | **-.013620** |
| 3A4 | .431687 | .431253 | .431253 | .000000 |
| Macro | .463295 | .462120 | **.458311** | **-.003810** |

1A2 and D6 both improve prospectively on the Live scored population, including
their Spearman correlations (.790049 to .794498 and .441151 to .450127).
Macro rank is 26; current 1A2/D6 ranks are 41/17. The result evaluates the
submitted model-and-placement combination, not isolated trunk features or
hidden-half generalization. Aggregate boards do not establish paired statistical
significance. The new payload is adopted as the scored Live control.

Submission event `f423ced2274141478574da6f45605036` attaches accepted/reflected
evidence. [Full metric comparison](../../artifacts/tracked/cyp-boltz-live-20261001T125026856964Z-71f3f209/metrics.csv)
includes MAE, R², Spearman and Kendall. The scoring monitor is now paused;
no new submission or training was performed.
