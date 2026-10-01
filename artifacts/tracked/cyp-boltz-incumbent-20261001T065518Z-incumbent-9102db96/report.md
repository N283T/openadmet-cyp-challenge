# Frozen CYP/HEM trunk in the submitted ensembles

This replaces the old UMAP regression screen for incumbent comparisons. All four regression readouts were refitted on the exact incumbent parent-series cohort, labels, confidence bounds and five folds. Existing trunk extraction and pooling are reused.

Each incumbent is reconstructed exactly from saved atomic predictions and fold weights. D6 preserves the submitted two-member core and moment-matched robust branch at their existing equal top-level weights. Additions and replacements operate on whole members; no arbitrary fractional blend sweep is used.

Placement selection and trunk candidate selection exclude the held fold. The 2C9 bank reruns exhaustive equal-subset selection; 3A4 reruns Caruana with inner step selection; 1A2 adds one equal member; D6 also evaluates Caruana over its two preserved branches and the selected trunk member. These are conditional comparisons on saved OOF, not fully nested base-model refits. Other-fold OOF models may have seen labels from the held fold during their original base training. Postselection and multiple comparisons limit uncertainty claims.

## Regression

Lower ST-RAE is better. Standalone minimum is best-of-grid; use the outer-selected placement and bank reselection rows to assess ensemble utility.

| Target | N | Incumbent | Standalone grid minimum | Outer-selected placement | Bank reselection |
|---|---:|---:|---:|---:|---:|
| cyp1a2 | 1372 | 0.670833 | 0.747107 | 0.669032 | 0.669839 |
| cyp2c9 | 1255 | 0.454711 | 0.594674 | 0.459567 | 0.458215 |
| cyp2d6 | 1457 | 0.819829 | 0.871612 | 0.819451 | 0.815745 |
| cyp3a4 | 2300 | 0.394567 | 0.512223 | 0.395694 | 0.396542 |

Same reselection algorithm without trunk provides a separate attribution control:

| Target | Original bank reselected | Bank plus trunk reselected |
|---|---:|---:|
| cyp1a2 | 0.670833 | 0.669839 |
| cyp2c9 | 0.454711 | 0.458215 |
| cyp2d6 | 0.820524 | 0.815745 |
| cyp3a4 | 0.394567 | 0.396542 |

Conditional parent-group bootstrap intervals (2,000 draws, fixed OOF denominator):

| Target | Scenario | Delta ST-RAE | Lower 95% | Upper 95% |
|---|---|---:|---:|---:|
| cyp1a2 | bank_extension_vs_reselected_control | -0.000994 | -0.003835 | 0.001926 |
| cyp1a2 | outer_selected_placement | -0.001801 | -0.005747 | 0.002147 |
| cyp1a2 | outer_selected_trunk_bank_reselection | -0.000994 | -0.003835 | 0.001926 |
| cyp2c9 | bank_extension_vs_reselected_control | 0.003504 | -0.001694 | 0.008737 |
| cyp2c9 | outer_selected_placement | 0.004856 | -0.001619 | 0.011411 |
| cyp2c9 | outer_selected_trunk_bank_reselection | 0.003504 | -0.001694 | 0.008737 |
| cyp2d6 | bank_extension_vs_reselected_control | -0.004780 | -0.010627 | 0.000895 |
| cyp2d6 | outer_selected_placement | -0.000378 | -0.009960 | 0.009220 |
| cyp2d6 | outer_selected_trunk_bank_reselection | -0.004084 | -0.009896 | 0.001710 |
| cyp3a4 | bank_extension_vs_reselected_control | 0.001975 | 0.000240 | 0.003817 |
| cyp3a4 | outer_selected_placement | 0.001127 | -0.000145 | 0.002532 |
| cyp3a4 | outer_selected_trunk_bank_reselection | 0.001975 | 0.000240 | 0.003817 |

## TDI diagnostic

The original 500 fitted models already match the incumbent UMAP folds and labelled cohorts; no duplicate TDI fitting was performed. Reconstructed D6 core/rank/tail stages and A4 seven-member rank ensemble are retained. Trunk readouts enter the D6 selected core or replace its PCA/raw block; A4 accepts a new equal member or replaces any named incumbent atom. AP-selected placement excludes the held fold. Historical labelled-fold MCC points are diagnostic only: the full 6,146-molecule trunk pool is unavailable, so no call-rate or deployment decision follows.

| Target | Scenario | AP | AUROC |
|---|---|---:|---:|
| d6 | current | 0.296860 | 0.580764 |
| d6 | outer_selected_AP_placement | 0.306124 | 0.591360 |
| a4 | current | 0.522733 | 0.805784 |
| a4 | outer_selected_AP_placement | 0.531759 | 0.810283 |

## Provenance

Regression fits completed: 1500. Incumbent reconstruction checks: `[{"target": "cyp1a2", "n": 1372, "incumbent_max_absolute_difference": 2.6645352591003757e-15}, {"target": "cyp2c9", "n": 1255, "incumbent_max_absolute_difference": 1.7763568394002505e-15}, {"target": "cyp3a4", "n": 2300, "incumbent_max_absolute_difference": 1.7763568394002505e-15}, {"target": "cyp2d6", "n": 1457, "incumbent_max_absolute_difference": 0.0}]`.

All fixed comparisons, keyed predictions and selections remain in their recorded runs. The previous UMAP scientific results are retained. No candidate, submission or upload is created.
