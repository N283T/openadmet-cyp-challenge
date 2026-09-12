# Nested external Caruana ensemble

Run `20260912T054406239560Z-96fce32b` fitted 24 candidate LightGBMs per endpoint and outer fold: eight feature views crossed with L1 and quantile losses at 0.35 and 0.65. Each candidate used 300 trees, 15 leaves, min_child_samples 20, learning rate 0.05 and L2=1. Three CPU threads; no GPU.

Weights were selected on a 20% inner group holdout of outer-training compounds, retaining all anchor-containing groups in inner training. Caruana ran 50 replacement steps and retained the iteration with lowest inner official ST-RAE. Selected candidates were refitted on all outer-training observations. Outer-fold labels were not used to fit weights or choose iterations.

All 6,384 OOF predictions matched the fixed baseline identities, endpoints, folds and confidence bounds. Audit `20260912T055242888594Z-56ff81b1` evaluates the nested candidate plus the already prescribed equal blends; no OOF weight optimization was performed.

## Results

| Method | All OOF ST-RAE | Neighbor ST-RAE |
|---|---:|---:|
| nested_caruana_external_lgbm | 0.648989 | 0.759972 |
| nested_caruana_external_lgbm__four | 0.647375 | 0.777195 |
| mean_four_pubchem600 | 0.658841 | 0.787963 |
| external3_class_mt_v1_ft_tabpfn__four | 0.645787 | 0.773727 |

Four means an equal blend of the candidate, FP+Mordred LGBM, predicted-log2fc TabPFN and frozen Admetica-embedding TabPFN. The original baseline is FP+Mordred, predicted-log2fc, Admetica and PubChem600. The FT reference comes from audit `20260912T054319015342Z-f96a88a3`.

The nested ensemble improves the original baseline, especially on the neighbor subset, but its best prescribed blend does not beat the current FT-embedding TabPFN four-model blend overall. It is retained as a diverse candidate, not adopted as a new global best.

## Endpoint comparison: all OOF

| CYP | Original baseline | Nested Caruana | Caruana four | FT TabPFN four |
|---|---:|---:|---:|---:|
| CYP1A2 | 0.735039 | 0.732722 | 0.727689 | 0.723007 |
| CYP2C9 | 0.557006 | 0.536659 | 0.536316 | 0.534756 |
| CYP2D6 | 0.866655 | 0.869832 | 0.864349 | 0.863033 |
| CYP3A4 | 0.476664 | 0.456743 | 0.461146 | 0.462350 |

## Inner-selected feature families

Mean weights average across all 20 endpoint/fold fits, including zero weights. Selected fits counts a family once per fit even when multiple losses were selected. These are training-selection frequencies, not a causal feature-importance measure.

| View | Mean weight | Selected fits /20 |
|---|---:|---:|
| mixed_logfc | 0.2052 | 15 |
| tree | 0.1841 | 17 |
| class | 0.1378 | 16 |
| tree_logfc | 0.1325 | 16 |
| class_logfc | 0.1321 | 15 |
| logfc | 0.0877 | 13 |
| mixed | 0.0661 | 11 |
| morgan | 0.0547 | 10 |

| Loss | Quantile | Mean weight | Selected fits /20 |
|---|---:|---:|---:|
| quantile | 0.65 | 0.4296 | 19 |
| regression_l1 | None | 0.3520 | 20 |
| quantile | 0.35 | 0.2184 | 17 |

The main gains versus the original baseline are CYP2C9 and CYP3A4. CYP1A2 changes little, and standalone Caruana worsens CYP2D6 slightly (0.869832 versus 0.866655). Compared with the FT TabPFN four-model reference, Caruana four is better only on CYP3A4; the other three endpoints explain its slightly worse overall score. The standalone nested model has the better neighbor macro score, 0.759972 versus 0.773727 for the FT reference.

The strongest mean weight is mixed embeddings plus predicted log2fc (20.5%), followed by external tree scores (18.4%). External tree scores alone or with log2fc account for 31.7% of average weight. All three loss variants are used broadly; upper quantile 0.65 receives 43.0%, L1 35.2%, and lower quantile 0.35 21.8%. These different prediction locations complement each other inside the ensemble; the frequencies do not establish that any one loss is universally superior.

Report tables and selection-frequency calculations are recorded in `20260912T055400786169Z-e8f02247`. No new training or tuning was performed during this audit.
