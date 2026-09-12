# Independent classification Caruana audit

The 17-candidate experiment is `classification/20260912T131834422982Z-d27648b2`. Independent replay **passed** all 84 selections and 84 metric rows, covering 35,567 OOF rows, 52,500 outer-weighted fold-test rows and 10,500 deployment rows. Its source code separates isolated inner selection, outer evaluation, and full-OOF deployment fitting. Candidate probabilities are keyed by compound, endpoint and fold; a graph model's probability takes precedence over its saved binary prediction.

Independent replay `f04fde09` checks all 17 candidate tables against canonical labels and original UMAP fold assignments. Every OOF row uses its canonical outer fold; every inner row has canonical fold `(outer_fold + 1) % 5`. This is stronger than merely checking disjoint compound IDs. Base-model training isolation was separately audited in the existing LGBM, TabPFN and graph audits; this replay checks the resulting probability table alignment and selection layer.

For each endpoint, recipe and outer fold, the replay independently enumerates all eligible subsets or all 50 greedy additions, verifies the selected addition, objective value, earliest best prefix, normalized weights and final MCC threshold. Its threshold calculation uses ascending cumulative class counts, independently of the producer's descending implementation. It then reconstructs every outer prediction, outer-weighted fold-test prediction and deployment prediction, and recomputes reported outer metrics. Full-OOF fitting weights and cutoffs are also replayed, but those fitted scores are not treated as CV estimates.

`subset2` and `subset3` mean the best equal-weight pair and triple selected by inner MCC from all combinations. They do not mean averaging the two or three individually highest-scoring single candidates. MCC selection uses the same `>=` cutoff convention as prediction. Sparse Caruana permits at most three distinct members, with replacement. AUROC and Brier Caruana select weights by their named objective and subsequently fit the MCC cutoff on the same isolated inner probabilities.

The first audit `e00ecc2e` failed before numerical replay because an audit input filename resolver matched multiple `oof.parquet` suffixes. The resolver was corrected to prioritize the exact declared filename; no experiment output or model code was changed.

The current experiment's MCC values reported by the producer are D6 0.047697 and 3A4 0.317735 for unrestricted MCC Caruana, and 0.061665 / 0.331549 for the sparse-three variant. These do not establish an improvement over the previous fixed means. The full-OOF fitting scores, approximately 0.1873 / 0.3915, are selection scores, not evidence of improved held-out accuracy.

## Expanded 18-candidate experiment

Experiment `8e3be966` adds the historical stack/Mordred TabPFN candidate with newly generated isolated inner fits, retaining its historical outer OOF predictions. Audit `22321032` **passed** all 108 selections and 108 outer metric rows, 45,729 OOF predictions, 67,500 outer-weighted fold-test predictions and 13,500 deployment predictions. It also verified all separate full-OOF fitting metrics. It applies the same canonical checks and adds verification of two fixed recipes:

- `fixed_historical`: historical model weight one, cutoff fitted on the matching newly isolated inner probabilities.
- `bagged_sparse3`: 20 stratified bootstrap samples of inner rows, RNG 731 reset per fit, each using the same 50-step sparse-three MCC Caruana. Class-zero and class-one draws are concatenated without a further shuffle. The weights are averaged, the largest three retained with deterministic name-order tie handling, and renormalized. The final cutoff is fitted on the original inner probabilities. The independent replay regenerates each sample, its selected weights and its fitted MCC.

Observed outer MCC values are:

| Recipe | D6 | 3A4 |
|---|---:|---:|
| Historical fixed control, honest inner cutoff | 0.026018 | 0.365149 |
| MCC Caruana | 0.071475 | 0.337668 |
| Sparse-three MCC Caruana | 0.061777 | 0.356736 |
| Bootstrap sparse-three | 0.084431 | 0.332283 |
| AUROC Caruana, inner MCC cutoff | 0.074920 | 0.342891 |
| Brier Caruana, inner MCC cutoff | 0.038395 | 0.344491 |

The earlier D6 fixed mean achieved 0.096228; none of these 18-pool fitted ensembles exceeds it in this comparison. The fixed historical model is stronger than these ensembles on 3A4. Weight fitting or bootstrap stabilization therefore does not itself establish an improvement. These comparisons are exploratory over the existing outer predictions; no independent significance or leaderboard-improvement claim is made.
