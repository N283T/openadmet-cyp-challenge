# Jeremy-expanded nested Caruana diagnostic

Fit `20260912T061036581812Z-288ee828` adds frozen Jeremy600 and Jeremy600+safe nested log2fc to the original eight views. All ten views use L1, quantile 0.35 and quantile 0.65 LightGBM models. Original 24 candidate inner predictions and all membership assignments are bitwise identical to the earlier run. Each of 20 endpoint/fold fits selects 50 replacement steps using its 20% inner group holdout, retains the best inner iteration and refits the selected models on all outer training labels.

The frozen upstream encoder is a separate 47-head pretraining checkpoint. Its training molecule manifest is unavailable; these scores are diagnostic and are not adopted as the formal goal result.

All 6,384 OOF identities, endpoints, folds, labels and bounds match the fixed baseline. The two fixed pairs use equal 0.5 weights. The fixed three-way blends add Jeremy fine-tuned embedding + log2fc LGBM with equal thirds. No OOF weight optimization is performed.

Official baseline audit: `20260912T062155979037Z-daa0d3f8`. Independent comparison and unchanged-candidate verification: `20260912T062158046652Z-cebbdbe3`.

The expanded standalone improves all four endpoints, especially 2C9 and 3A4. However, substituting it into the fixed pair or three-way blend slightly worsens overall ST-RAE, while improving the neighbor subset. The standalone gain therefore does not translate into a better overall fixed ensemble; the existing original-Caruana three-way blend remains the stronger overall comparator. The two Jeremy views receive 55.1% of mean inner-selected weight, indicating that the new information is already heavily represented before adding the separate Jeremy models.

| Method | All ST-RAE | Neighbors ST-RAE |
|---|---:|---:|
| expanded | 0.628522 | 0.722938 |
| original | 0.648989 | 0.759972 |
| frozen_jeremy_pfn | 0.649557 | 0.702130 |
| expanded_jeremy_fixed_pair | 0.613236 | 0.683722 |
| original_jeremy_fixed_pair | 0.611584 | 0.689559 |
| mean_four_pubchem600 | 0.658841 | 0.787963 |
| expanded_fixed_three | 0.605946 | 0.687956 |
| original_fixed_three | 0.605039 | 0.696275 |

## Endpoint results

| CYP | Expanded | Original | Expanded + Jeremy PFN | Original + Jeremy PFN |
|---|---:|---:|---:|---:|
| CYP1A2 | 0.721577 | 0.732722 | 0.698377 | 0.693449 |
| CYP2C9 | 0.500926 | 0.536659 | 0.487526 | 0.490043 |
| CYP2D6 | 0.860072 | 0.869832 | 0.845556 | 0.845091 |
| CYP3A4 | 0.431512 | 0.456743 | 0.421482 | 0.417753 |

## Inner-selected families

Weights average across all 20 fits, including zero weights.

| View | Mean weight | Selected fits |
|---|---:|---:|
| jeremy | 0.3728 | 19 |
| jeremy_logfc | 0.1779 | 14 |
| tree | 0.1190 | 13 |
| tree_logfc | 0.0787 | 10 |
| mixed_logfc | 0.0668 | 9 |
| class | 0.0622 | 11 |
| class_logfc | 0.0408 | 6 |
| mixed | 0.0351 | 4 |
| logfc | 0.0331 | 10 |
| morgan | 0.0137 | 3 |
