# Direct-inhibition Emax auxiliary preparation

Run: `20260912T055345388545Z-8ab226fa` (regression). Preparation only; no model was trained.

The audited official Emax CSV provides **6,525 finite direct-inhibition Emax labels across 4,905 compounds**, suitable for four additional multitask heads. Values retain their original scale and sign. The TDI-condition Emax and `is_TDI` columns were not selected as targets.

## Interpretation and source

The [official dataset card](https://huggingface.co/datasets/openadmet/cyp-challenge-train-test) describes Emax as maximal effect and supplies direct-inhibition and TDI-condition estimates with 95% confidence intervals. The card was saved as `official-dataset-card.md` in this run, with its hash in `summary.json`.

The selected columns are named `CYP*_EmaxVsPosCtrl_direct_inhibition`. They are Emax estimates, not pIC50, and the `VsPosCtrl` naming indicates comparison to the positive control. The checked official card does not give an exact normalization equation. Consequently, this preparation does not reinterpret the negative values as concentrations, invert their sign, clip them to [-1, 0], or convert them to percentages. The native values and confidence intervals are retained.

The source CSV is the downloaded snapshot in regression run `20260911T141511080701Z-fd93c411`, compared against `data/raw/cyp-challenge-TRAIN_Emax.csv`. The downloaded version has 6,145 rows versus the original 6,146. Its only removed record is `OCNT-2329004`, SMILES `N.N.[Cl][Pt][Cl]`; there are no added records. The original source remains untouched. Although the descriptive dataset card still lists 6,146 Emax compounds, the actual audited CSV determines the prepared dataset.

## Identity and label checks

- Exact joint matching on `Molecule_Name` and raw `SMILES` against the shared aliases table mapped all 6,145 source rows.
- There are no duplicate molecule names or standardized compound IDs in the selected source.
- Every mapped compound exists in the shared compounds and current parent assignments.
- Of 24,580 potential direct Emax values, 18,055 are missing; missing observations remain in `source-rows.parquet` and are omitted from the point-target file.
- There are no censor tokens, nonnumeric values, nonfinite numeric estimates, inverted confidence intervals, or estimates outside their supplied intervals.
- Aggregation is specified as the median within standardized compound/head; no actual duplicate aggregation was necessary.

## Target distributions

| CYP | Labels | Mean | SD | Minimum | Median | Maximum |
|---|---:|---:|---:|---:|---:|---:|
| 1A2 | 1,412 | -0.98487 | 0.05839 | -1.15104 | -0.99354 | -0.45768 |
| 2C9 | 1,285 | -1.01613 | 0.08544 | -1.26977 | -1.02074 | -0.49547 |
| 2D6 | 1,493 | -1.01582 | 0.03241 | -1.08307 | -1.02680 | -0.82721 |
| 3A4 | 2,335 | -0.97723 | 0.04100 | -1.05616 | -0.97540 | -0.36766 |

These labels are concentrated near -1, especially for 2D6 and 3A4. Training-time scaling must be fitted separately per head using only the inner-training compounds. Adding Emax adds a type of supervision; the target counts do not establish additional chemical coverage beyond existing pIC50 labels.

## Outputs and integration

`artifacts/runs/regression/20260912T055345388545Z-8ab226fa/outputs/targets.parquet` has exactly the compatible columns `compound_id`, `head`, and `value`. Its four head names are the original direct-inhibition Emax column names.

The run also retains raw source values and confidence intervals, per-target provenance, fold/anchor counts, distribution quantiles, source differences, and the official dataset card. Any subsequent fine-tuning must exclude an outer compound's Emax labels alongside all its other challenge labels. No fine-tuning code or existing runs were changed here.
