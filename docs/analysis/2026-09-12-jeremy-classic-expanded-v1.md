# Jeremy + FP/Mordred expanded nested Caruana

Fit `20260912T062329198020Z-60f59445` uses 36 candidate LightGBMs with L1, quantile 0.35 and quantile 0.65 losses. The previous 30 candidate inner predictions and all membership assignments are bitwise identical to the earlier run. Each of 20 endpoint/fold fits selects 50 replacement steps using its 20% inner group holdout, retains the best inner iteration and refits the selected models on all outer training labels.

The frozen upstream encoder is a separate 47-head pretraining checkpoint. Its training molecule manifest is unavailable; these scores are diagnostic and are not adopted as the formal goal result.

All 6,384 OOF identities, endpoints, folds, labels and bounds match the fixed baseline. The two fixed pairs use equal 0.5 weights. The fixed three-way blends add Jeremy fine-tuned embedding + log2fc LGBM with equal thirds. No OOF weight optimization is performed.
For FP+Mordred variants, all 40 saved preprocessing artifacts were independently checked: exact inner/outer training IDs, finite/nonconstant descriptor filters and training medians. The original score comparator below remains the 24-candidate ensemble.

Official baseline audit: `20260912T071550737929Z-945cbf90`. Independent comparison, unchanged-candidate checks and preprocessing verification: `20260912T071553536042Z-9f26b307`.

Adding FP+Mordred barely changes the result: standalone ST-RAE worsens by 0.000204 versus 30 candidates, while the fixed three-way blend improves by only 0.000054. Both remain behind the original 24-candidate fixed three-way blend at 0.605039. The new descriptor views receive just 2.69% combined average inner-selected weight. This arm provides no reason to replace the existing best fixed blend.

| Method | All ST-RAE | Neighbors ST-RAE |
|---|---:|---:|
| expanded | 0.628726 | 0.725438 |
| original | 0.648989 | 0.759972 |
| frozen_jeremy_pfn | 0.649557 | 0.702130 |
| expanded_jeremy_fixed_pair | 0.613124 | 0.684693 |
| original_jeremy_fixed_pair | 0.611584 | 0.689559 |
| mean_four_pubchem600 | 0.658841 | 0.787963 |
| expanded_fixed_three | 0.605892 | 0.688716 |
| original_fixed_three | 0.605039 | 0.696275 |
| previous_expanded | 0.628522 | 0.722938 |

## Endpoint results

| CYP | Expanded | Original | Expanded + Jeremy PFN | Original + Jeremy PFN |
|---|---:|---:|---:|---:|
| CYP1A2 | 0.721577 | 0.732722 | 0.698377 | 0.693449 |
| CYP2C9 | 0.501186 | 0.536659 | 0.487176 | 0.490043 |
| CYP2D6 | 0.860079 | 0.869832 | 0.845584 | 0.845091 |
| CYP3A4 | 0.432062 | 0.456743 | 0.421359 | 0.417753 |

## Inner-selected families

Weights average across all 20 fits, including zero weights.

| View | Mean weight | Selected fits |
|---|---:|---:|
| jeremy | 0.3646 | 19 |
| jeremy_logfc | 0.1805 | 14 |
| tree | 0.1110 | 12 |
| tree_logfc | 0.0764 | 9 |
| mixed_logfc | 0.0643 | 9 |
| class | 0.0641 | 11 |
| class_logfc | 0.0385 | 7 |
| mixed | 0.0312 | 4 |
| logfc | 0.0305 | 9 |
| classic | 0.0229 | 5 |
| morgan | 0.0121 | 3 |
| classic_logfc | 0.0040 | 1 |
