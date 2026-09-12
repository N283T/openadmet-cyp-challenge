# Alternative primary-screen featurizer

## Prespecified comparison

Replace the nested CheMeleon primary-screen log2FC model with fixed LGBM predictors trained on frozen Jeremy medium embeddings (600 columns) and the original 47 upstream predictor outputs. The downstream controlled comparison uses the same 600 frozen embeddings plus four predicted primary endpoints, with the existing fixed pIC50 LGBM readout (300 trees, 15 leaves, L1). CPU thread count is three for both controls.

Primary LGBM: squared regression, 400 trees, 31 leaves, minimum child size 20, learning rate 0.03, L2=1, deterministic seed 42. All four actual primary-log2FC endpoints are modeled independently; there is no HPO, no label imputation and no test-label use.

Ten excluded-fold pairs produce 40 saved models and predictions for all 6,896 compounds. For outer fold f, a nonanchor training recipient in fold j receives the model excluding f and j; outer-validation recipients receive the average of the four models excluding f. Anchors remain eligible for fitting and receive that same four-model average. All nonanchor outer and recipient compounds are excluded regardless of endpoint label availability. Model membership masks and independent recipient-wise routing reconstruction are saved.

The old CheMeleon auxiliary excluded anchors in its pair models as well; this alternative follows the retained-anchor policy requested for this experiment. Thus this is a comparison of complete featurizer recipes, including that small anchor-policy difference, not architecture alone.

## Run records

- Primary features: shared `20260912T082554377951Z-d272774c`.
- Old-feature pIC50 readout: regression `20260912T082657021636Z-bedc430e`; audit `20260912T082733509200Z-2a8406be`.
- Premature new-feature readout: `20260912T082735890987Z-c61c052e`, failed before model fitting because primary output was not ready. Retained. Launcher now checks upstream completion and supports `--mode new` to avoid repeating the completed old control.

Scripts: `focus_primary.py`, `launch_focus_primary_readouts.py`, `focus_primary_compare.py`. Configuration directory: `configs/focus-primary/`.

## Primary prediction accuracy

Matched nonanchor measured primary labels: 4,313 compounds per endpoint, same keys in both recipes.

| Endpoint | Old RMSE | New RMSE | Old Pearson | New Pearson |
|---|---:|---:|---:|---:|
| 1A2 | 0.557584 | 0.531817 | 0.583397 | 0.626777 |
| 2C9 | 0.313764 | 0.274769 | 0.733605 | 0.802904 |
| 2D6 | 0.878571 | 0.816909 | 0.610286 | 0.673317 |
| 3A4 | 0.932197 | 0.840917 | 0.728727 | 0.785453 |

All four improve, providing direct evidence that the previous primary featurizer had room to improve. This does not establish that a better primary predictor must always improve pIC50, so downstream controls follow.

## Controlled pIC50 LGBM result

New-feature readout `20260912T083010595992Z-1e7001de`, audit `20260912T083035807114Z-7fed236b`. Matched comparison `20260912T083048351737Z-03e8800c` evaluates all 6,384 pIC50 OOF rows and predetermined equal additions to the submitted ensemble; no weights are fitted.

| Endpoint | Old primary | New primary | Submitted ensemble | Equal six models | Equal half submitted/new |
|---|---:|---:|---:|---:|---:|
| 1A2 | 0.716652 | 0.691998 | 0.683335 | 0.677549 | 0.674183 |
| 2C9 | 0.512844 | 0.494360 | 0.463781 | 0.457919 | 0.458106 |
| 2D6 | 0.858419 | 0.842287 | 0.834741 | 0.830297 | 0.828239 |
| 3A4 | 0.463438 | 0.442934 | 0.403800 | 0.401670 | 0.407314 |

“Equal six” means `(5 × submitted five-model mean + new readout) / 6`; it improves every endpoint and reaches macro ST-RAE 0.591858. The controlled old/new readout improves every endpoint as well. These are useful but modest OOF gains, not evidence of a solved leaderboard gap.

## Controlled TabPFN check

Because actual primary prediction improved on all four endpoints, both old and new auxiliaries were also evaluated with the same frozen 600 embeddings and TabPFN 2.6 (8 estimators, 604 total dimensions). No old equivalent was available, so both controls were run.

- New: `20260912T083032513346Z-568bc526`; audit `20260912T083319485280Z-af814bd1`.
- Old: `20260912T083330269990Z-9f0909f1`; audit `20260912T083554998915Z-e22c6f22`.
- Matched comparison: `20260912T083605349220Z-70587775`.

| Endpoint | Old-primary TabPFN | New-primary TabPFN | Submitted plus new TabPFN, equal six |
|---|---:|---:|---:|
| 1A2 | 0.721231 | 0.705017 | 0.679241 |
| 2C9 | 0.526740 | 0.495867 | 0.460225 |
| 2D6 | 0.878098 | 0.874374 | 0.837052 |
| 3A4 | 0.448111 | 0.439395 | 0.403022 |

All endpoints improve in this paired comparison too. However, the new TabPFN readout is weaker than the new LGBM readout on both focus targets, and adding it to the submitted ensemble worsens 2D6. Retain the new-primary LGBM readout as the useful candidate from this experiment; the TabPFN arm is comparison evidence, not a preferred ensemble addition.

All launched fits completed except the explicitly retained premature input failure above. No submission or production pointer change was made. Ruff check and format check pass for the new scripts.
