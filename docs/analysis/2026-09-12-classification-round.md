# Classification round during regression cooldown

The user requested work on classification while waiting for the regression upload
cooldown. This round completed 17 individual recipes and seven fixed probability
means for the two official TDI endpoints. It produced saved models, OOF and test
probabilities; no classification or regression upload was made.

## Evaluation and training

Use the original five UMAP folds, shared run `edba6fc0`, throughout. CYP2D6 has
1,497 labels (324 positive); CYP3A4 has 3,584 (764 positive). Test contains 750
compounds. The explicit endpoint whitelist excludes the unscored 1A2/2C9 TDI
measurements. See the [inventory](2026-09-12-classification-inventory.md).

For every outer fold, a designated remaining fold supplies inner validation.
Neither its threshold fit nor its model training uses the outer validation labels.
Models are then fitted on all available outer-training labels. The results below
pool the five outer folds. They are screened recipe comparisons, not a nested
estimate of selecting the winning recipe from this entire round.

- LGBM: six feature views (Morgan, Mordred, their concatenation, frozen external
  classification embeddings, mixed embeddings, classification embeddings plus
  Mordred), each with plain or balanced class weights. Each of 120 outer fits
  has a separate inner model; 400 trees, 15 leaves, learning rate 0.03, minimum
  child size 20 and L2 2. Feature encoders were trained only on public data.
- TabPFN 2.6: three frozen feature views, external classification, external mixed
  and original CheMeleon. Train-fold scaling and PCA256; four estimators. Thirty
  inner and thirty outer fits. See the [report](2026-09-12-tdi-tabpfn.md).
- Chemprop: actual TDI fine-tuning from the public classification checkpoint,
  either two TDI heads alone or with ten regression auxiliaries. Inner model
  selection followed by fresh outer refits for the selected epoch counts.
  All auxiliary labels exclude held-out rows. Ten inner fits and ten outer
  refits. See [training](2026-09-12-classification-joint-finetuning.md).
- Seven equal-weight means, declared without fitted ensemble weights: classical
  LGBM, all plain LGBM, all balanced LGBM, TabPFN, graph, classical plus TabPFN,
  and all plain LGBM plus TabPFN plus graph. Thresholds are selected on the
  corresponding mean of the separately trained inner probabilities.

All individual classifiers also retain fixed-0.5 metrics. LGBM/ensemble threshold
search uses unique inner probabilities plus 0.5; the TabPFN screen uses its declared
grid. Both exclude outer labels, but the threshold search spaces differ.

## Results

Higher is better. MCC uses inner-selected thresholds unless explicitly indicated.

| Recipe | D6 MCC | D6 AUROC | A4 MCC | A4 AUROC |
|---|---:|---:|---:|---:|
| Plain LGBM six-view mean | 0.091415 | 0.565146 | 0.315875 | 0.770557 |
| Balanced LGBM six-view mean | 0.078000 | 0.555453 | **0.338112** | 0.777900 |
| TabPFN three-view mean | 0.069503 | 0.569830 | 0.296691 | 0.743326 |
| TDI-only graph | 0.049559 | 0.520171 | 0.111170 | 0.634316 |
| Graph with auxiliaries | 0.025323 | 0.529741 | 0.121850 | 0.648071 |
| Plain LGBM + TabPFN + graph mean | **0.096228** | 0.566186 | 0.314955 | 0.768413 |

The best individual LGBM A4 inner-threshold MCC is 0.321268 from balanced
external-classification embeddings plus Mordred (AUROC 0.778318). D6's best
individual inner-threshold MCC is 0.084640 from balanced external-classification
embeddings, while its AUROC is only 0.528504. That mismatch is another reason to
inspect probability ordering as well as thresholded MCC.

At the fixed 0.5 cutoff, the balanced LGBM mean has MCC 0.078473 / 0.298870
(D6 / A4); the all-family plain mean has 0 / 0.121389. Threshold selection helps
the latter, but has not created strong D6 probability ordering.

The previously saved CheMeleon + predicted-log2fc + Mordred baseline, audited in
`243a1e23`, retains A4 AUROC **0.793282** and fixed-0.5 MCC **0.258678**. Its
old all-OOF-selected threshold MCC 0.383103 is a fitting score and must not be
compared as an honest-inner score. Its submitted A4 LB MCC was 0.430653, while
D6 LB MCC was 0.001644. These are LB results, not these UMAP scores. Keep that
baseline in the candidate pool; this round has not established a superior A4
submission. The new A4 balanced-mean fixed-0.5 MCC is higher but AUROC is lower.

## Artifacts, failures and validation

| Run suffix | Role |
|---|---|
| `269e598d` | Explicit official classification label and split inventory |
| `76facc71` | Twelve LGBM recipes; 60,972 OOF and 90,000 fold-test rows |
| `8aeb9a25` | Independent saved-booster reconstruction; all 120 inner thresholds |
| `51023dac` | Three TabPFN outer-fit recipes |
| `397ddee7` | TabPFN honest inner fits and thresholds |
| `401ba682` | Independent TabPFN audit |
| `1da37202` | Two joint Chemprop recipes |
| `7e38be86` | Independent joint masks, scaling, epochs and probability audit |
| `b0a05986` | Seven fixed means, 35,567 OOF and 10,500 averaged test rows |
| `da4f1375` | Independent ensemble reconstruction, 70 thresholds and 168 metrics |

The first LGBM attempt `46d01526` was deliberately interrupted when a suffix
filter included unscored endpoints. It was replaced by the explicit endpoint
whitelist before the complete round. Ensemble attempt `a64253d1` failed at final
concatenation because label dtypes differed between families. Canonical integer
label/fold and float probability types fix this in `b0a05986`. Both failed runs
remain recorded; neither supplied reported results.

[LGBM audit](2026-09-12-classification-lgbm-audit.md) reconstructed saved-model
predictions and thresholds. [Joint audit](2026-09-12-classification-joint-audit.md)
verified all task masks, scalers, selected/refit epochs, thresholds and exports.
TabPFN was checked similarly. These checks found no split, label-alignment or
prediction-export defect explaining the weak D6 results.

The [ensemble audit](2026-09-12-classification-fixed-ensemble-audit.md) independently
reproduced all seven means, canonical inner/outer assignments, 70 thresholds,
168 metrics and the test probability means, median cutoffs and binary decisions.

Fixed-mean test artifacts are five-fold probability averages with the median of
the five inner thresholds. They are not all-label refits. Each underlying
recipe's fold-test probabilities and saved models remain available.

## Decision

Retain the historical A4 baseline and the new balanced-LGBM mean as candidates.
The broad means are useful D6 comparisons, but MCC around 0.096 and AUROC around
0.57 do not constitute a strong solution. Direct TDI fine-tuning and the tested
auxiliary loss do not justify replacing the tabular candidates. A next focused
comparison should combine the existing stronger A4 baseline with these new
probabilities using matching inner-fit predictions; do not choose a cutoff on
the reported outer evaluation labels. Regression all-endpoint submission assembly
and the cooldown status remain separately tracked in RL0052.
