# Source-faithful community Chemprop training reproduction

The user authorized reproducing the learning stage after the preceding experiments
only adapted individual ideas to the existing external600/LGBM pipeline.

## Reproduction contract

- Use the author's public preprocessed data, rather than our previous external
  curation. The historical union18 reconstructs exactly 31,670 molecules and
  161,400 labels. Challenge8 has 4,905 molecules; D6single has 1,493 labels.
- Pin the historical 18-head definition, rather than silently following HEAD's
  expanded union. D6isoform6 separately retains its later TDI/Emax head/data scope.
- Execute original Workbench native Chemprop model/training functions, including
  scratch initialization, hidden size700, depth5, FFN300×2, NormAggregation100,
  BatchNorm, dropout0.1, standardized MAE, inverse-count primary weights,
  source auxiliary weights, and Noam learning-rate schedule1e-4→1e-3→1e-4.
- Retain maximum400 epochs, patience40, batch64 and mixed16 precision, with full
  passes over the public data. Save native best-validation checkpoints.
- Evaluate on retained parent folds, with external scaffold fold assignments and
  all-connectivity exclusion. This split is explicitly different from the author.
  Add the author's exact scaffold split for the small D6single diagnostic.
- Average direct Chemprop predictions: union18 and challenge8 for other CYPs,
  those two plus D6isoform6 and D6single for CYP2D6. Test predictions follow the
  source's five-fold ensemble; the source's “100%” label does not mean a single
  model fitted on every row. No LGBM or placement correction is substituted.

See [data provenance](2026-09-12-scp-reproduction-data.md),
[independent recipe audit](2026-09-12-scp-source-recipe-audit.md), and
[metric-definition comparison](2026-09-12-scp-metric-comparison.md).

## Execution integrity

The first adapter attempt3d49be88 failed on recorded input path prefixes and is
retained. Native training then ran on full source data. An export adapter error
used bare CYP names against canonical suffixed endpoint names: this left empty
long OOF files while native training, raw per-fold predictions and checkpoints
were intact. Separate recorded export runs recover canonical outputs from those
raw predictions, with hard checks for1,457/6,384 OOF rows and750 test molecules per
endpoint. Original run files are preserved. Evaluations consume the recovered
exports, not the invalid original long OOF tables.

## Split-only control

Recorded comparison `20260912T105908385052Z-eb5c8d06` uses the same 1,457
canonical D6 labels for all three rows below; both native models use the same
published model and training defaults.

| Model | Official ST-RAE | Author-local ST-RAE | Spearman | MAE |
|---|---:|---:|---:|---:|
| Native single, parent folds | 0.906826 | 0.597431 | 0.361932 | 0.584960 |
| Native single, source scaffold folds | 0.885683 | 0.583502 | 0.400917 | 0.577057 |
| Previous global ensemble | 0.825436 | 0.543811 | 0.472752 | 0.554068 |

The source split helps this single model modestly but does not close the gap to
our previous ensemble. On its full 1,493-label cohort the source single obtains
0.886962 official ST-RAE and 0.402967 Spearman. The source split preserves two
connectivity groups split across folds; the parent runs have zero such overlap.
This control does not establish the split effect for the external multitask
ensemble, and neither single model is the blog's four-member ensemble.

## Individual completed members

Recorded early comparisons `287f0fa4` (challenge8) and `a1cb2ca5` (D6isoform6)
use the retained parent OOF keys. D6 challenge8 scores official ST-RAE0.890313
and Spearman0.411415. D6isoform6 scores0.869874 and0.415397. Public D6 tasks
improve the native single/auxiliary models, but neither beats the prior global
ensemble0.825436/0.472752. These are individual models, not the final mean.

## Completed four-member result

All 25 fits completed: four parent-CV recipes × five folds, plus five source-split
single controls. Final evaluation `20260912T112404158683Z-0c4a365a` consumes only
independently recovered canonical exports. The prescribed average has no fitted
weights or placement correction.

| Endpoint | Previous global official ST-RAE | Native fixed mean official ST-RAE | Previous Spearman | Native Spearman |
|---|---:|---:|---:|---:|
| CYP1A2 | 0.670794 | 0.695111 | 0.611303 | 0.595188 |
| CYP2C9 | 0.443453 | 0.528750 | 0.742306 | 0.704603 |
| CYP2D6 | 0.825436 | 0.823983 | 0.472752 | 0.472811 |
| CYP3A4 | 0.391582 | 0.442451 | 0.825636 | 0.810914 |
| Macro ST-RAE | 0.582816 | 0.622574 | — | — |

D6's average improves over every native individual: single0.906826,
aux8 0.890313, isoform6 0.869874, union18 0.920374. Its Spearman0.472811 also
exceeds the best native member union18's0.431863. Nevertheless, the result is
essentially tied with the previous global ensemble, not a substantive gain;
no significance claim is made for the 0.001453 ST-RAE difference. All three
other endpoints worsen. Do not adopt this as a replacement for the current model.

Under the author's local denominator, the same native mean scores macro0.378030
versus previous0.356179. Those numbers are not comparable to the official scores
in the table. The blog reports D6 ensemble Spearman0.503 vs best-member0.445;
our analogous result0.473 vs0.432 is closer in ranking quality than the raw
leaderboard error gap might suggest, but these are different CV cohorts. The
parent-to-source split effect measured for the single cannot be assumed to carry
over to the ensemble.

The unplaced D6 test mean is4.687115 with standard deviation0.360755. Thus this
training reproduction has not itself introduced a downward placement correction.
The blog's live LB0.4378 includes placement, while this run has no new LB score.
Do not treat its development OOF result as a prediction of that leaderboard score.
Source: https://supercowpowers.github.io/workbench/blogs/cyp_challenge/.

Final artifacts contain6,384 OOF rows and3,000 test rows, with IDs/endpoints/folds;
each member test output averages five trained models. Artifacts and checkpoints
remain in recorded runs. No challenge submission was made. Infrastructure
validation: Ruff passes; pytest37 passed,7 skipped (optional tests).

## Residual diversity

Independent diversity diagnostic6eb3beec finds native-member prediction Pearson
correlations0.560–0.707 but residual correlations0.876–0.942. The fixed mean's
prediction correlation with the prior ensemble is0.8732, and its residual
correlation is0.9776. Thus the average cancels some individual model error but
retains almost the same difficult compounds as the previous pipeline. D6 MAE
improves0.5541→0.5508 while RMSE worsens0.8085→0.8143; rank correlation is flat.
OOF mean bias is small (+0.0226 vs previous+0.0282). These OOF observations do not
establish the mean bias on the hidden blind labels.

Final union export audit08b26a21 passes exact native-array correspondence and
reconstructed five-fold test means, completing the five-export audit set.
