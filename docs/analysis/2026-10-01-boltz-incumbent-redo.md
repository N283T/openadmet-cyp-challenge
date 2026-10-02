# Incumbent-matched frozen CYP/HEM trunk redo

The previous suite `20260930T232717Z` evaluated regression on the original UMAP
cohort. Submitted regression models use a retained parent-series cohort. Its
matched molecular controls do not answer whether trunk models help the submitted
ensembles. Preserve those runs as exploratory representation evidence, but do not
use their scores to claim incumbent improvement or ensemble utility.

The user requested a redo. Suite `20261001T065518Z-incumbent` reuses the four
completed pooling runs and refits all 1,500 regression fits. The exact cohort,
labels, confidence bounds and fold numbers come from the submitted OOF bank,
with assignments from `20260911T104703655337Z-a3923255` (`strategy=series`).
Preflight checks exact trunk coverage and fold agreement before any fitting.
Each screen saves `evaluation-cohort.parquet` and validates training/held groups.

| Target | Labelled cohort | Incumbent OOF source and recipe |
|---|---:|---|
| 1A2 | 1,372 | `46314b50`, `1A2/add_dense_and_fast_pool`, ten-member equal mean |
| 2C9 | 1,255 | `e1560701`, outer-selected equal subset |
| 2D6 | 1,457 | `0dbf57d4`, parent, delta .5, constraint .35, moment match, weight .5 |
| 3A4 | 2,300 | `dfa02730`, Caruana bank |

The 12 pooling views, trunk-only and trunk-plus-molecular modes, three readouts
and train-fold-only imputation/selection are retained. The same molecular-only
controls remain diagnostics, but the primary reference is the submitted ensemble.

## Ensemble integration

Reconstruct the incumbent exactly from saved atoms and fold weights, requiring
maximum prediction disagreement below `1e-9`. Require identical molecule IDs,
endpoint, fold, labels and confidence bounds for every trunk comparison.

- Evaluate every trunk readout as a whole member added to an equal-mean bank,
  or replacing each incumbent member at its saved fold weight.
- Preserve D6's submitted two-member core and its moment-matched robust branch
  as two explicit top-level members. The robust branch is recovered algebraically
  from its retained .5 blend and base OOF; additions/replacements preserve this
  hierarchy rather than silently flattening it.
- Select a placement on the other four OOF folds; include keeping the incumbent
  unchanged as a candidate, and score the held fold separately.
- Select a trunk candidate on other folds and rerun 2C9's exhaustive equal-subset
  selection and 3A4's Caruana selection (step counts 5/10/20/50 chosen by inner
  OOF-fold evaluation). 1A2 keeps its equal-member rule. D6 additionally tests
  Caruana on its two preserved branches plus the selected trunk readout.
- Save keyed OOF, all fixed comparisons, selections and conditional grouped
  bootstrap differences (2,000 draws). Do not report the global grid minimum as
  the estimated utility of a selected ensemble.
- Repeat the same bank selection with no trunk column to separate changing the
  ensemble selector from the incremental value of the trunk member.

These are conditional saved-OOF comparisons. Other-fold OOF predictions were
produced by base models that can include held-fold labels in their training;
selector exclusion alone does not establish fully nested base-model evaluation.
No independent confirmation or multiplicity-corrected improvement is claimed.

## TDI

The existing 500 trunk fits already use matching UMAP folds and the full labelled
cohort (D6 1,497; A4 3,584). Reuse these predictions without fitting duplicates.
Reconstruct D6's selected core, rank and tail stages, and A4's seven rank-normalized
members exactly. Trunk scores enter the D6 core as one member or replace its
PCA/raw block; A4 accepts one equal member or replacement of a named atom.
Select the AP placement on other folds, retaining the incumbent as an option.
Historical labelled-fold MCC rates are retained only for diagnosis. No operating
rate or submission decision is supported without features for the full 6,146 pool.

All substantial fitting, placement and reporting use recorded runs with explicit
aliases/configuration and `CYP_OUTPUT_DIR`. Compact report, metrics, selections
and provenance are promoted to `artifacts/tracked/`. No candidate or upload.

## Completed result

All 1,500 regression fits completed without failures. Final placement run:
`regression/20261001T111723965954Z-04ff278c`; recorded report:
`shared/20261001T111739573189Z-4e969e95`. The final analysis adds a no-trunk bank
reselection control; prior placement/report attempts remain retained. A report
promotion path error was repaired without repeating the completed recorded report.

| Target | Incumbent ST-RAE | Outer-selected placement | Extended bank reselection |
|---|---:|---:|---:|
| 1A2 | .670833 | .669032 | .669839 |
| 2C9 | .454711 | .459567 | .458215 |
| 2D6 | .819829 | .819451 | .815745 |
| 3A4 | .394567 | .395694 | .396542 |

1A2 and D6 have small conditional point improvements, with paired 95% intervals
spanning zero. D6's selected trunk member receives 20–25% Caruana weight in every
held fold. 2C9 worsens and only three folds select trunk in bank reselection.
3A4 assigns the trunk member zero final weight in every fold; including a trunk
candidate in inner step selection still changes the selected ensemble and worsens
held-fold performance. Retain the incumbent 2C9/3A4 recipes; keep 1A2/D6 results
as uncertain exploratory evidence rather than promoting a new model.
