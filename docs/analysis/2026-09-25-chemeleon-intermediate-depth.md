# CheMeleon intermediate message-passing depth screen

## Question and design

Does an intermediate frozen CheMeleon state transfer better than the usual final
state to CYP regression or TDI classification? Extraction run
`20260925T065503727832Z-1c269e23` reuses the same pretrained 2048-wide,
six-step bond-message-passing checkpoint as the current frozen CheMeleon feature
table. It takes directed-bond states after steps 1, 3 and 6, applies the same
frozen final atom projection `W_o` at each step, then mean-pools atoms. This
produces a 2048-dimensional molecular vector at every depth. The final depth
matches Chemprop's forward/fingerprint path to maximum absolute difference
`2.98e-7`; it matches the previously exported CheMeleon table over all 6,896
compounds to `4.77e-7`. The 1/3 views are probes, not native Chemprop outputs.

Readout run `20260925T065626160287Z-8017d516` fits the same LightGBM recipe at
each depth, independently for each endpoint, on the fixed UMAP five-fold split.
Regression uses 200-tree L1 LightGBM and the official soft-threshold RAE;
classification uses 200-tree binary LightGBM and reports pooled AUROC/AP. No
auxiliary predictions, extra descriptors, PCA, ensemble components, test-label
feedback or model selection enter this comparison. All 6,525 regression and
5,081 classification OOF rows are present at each depth with unique keys.

## Results

Lower ST-RAE is better; higher AUROC/AP is better.

| Regression endpoint | Depth 1 ST-RAE | Depth 3 | Final depth 6 | Best |
|---|---:|---:|---:|---:|
| CYP1A2 | .83687 | **.83096** | .83912 | 3 |
| CYP2C9 | .75551 | **.73661** | .73858 | 3 |
| CYP2D6 | **.90101** | .92735 | .91246 | 1 |
| CYP3A4 | .62263 | **.62135** | .62548 | 3 |
| Four-endpoint macro | .77900 | .77907 | **.77891** | 6 |

| Classification endpoint / metric | Depth 1 | Depth 3 | Final depth 6 | Best |
|---|---:|---:|---:|---:|
| CYP2D6 AUROC | **.54004** | .52855 | .53712 | 1 |
| CYP2D6 AP | **.27309** | .25101 | .25133 | 1 |
| CYP3A4 AUROC | **.73123** | .70776 | .72705 | 1 |
| CYP3A4 AP | .40573 | .39126 | **.42105** | 6 |

Intermediate states sometimes beat the final state as *standalone* features, but
no depth wins across tasks. The regression macro values differ by less than
`.0002`, and none of these standalone readouts approaches the deployed ensemble.
Depth 1 gives the clearest exploratory classification D6 AP gain (`+.02176`
versus final), while 3A4 AP favors the final representation. This screen does
not justify replacing a submitted component or predicting leaderboard benefit.

## Similarity of held-fold predictions

The table uses pooled OOF Spearman correlation after joining on compound ID,
endpoint and fold. Pearson correlations lead to the same qualitative conclusion.

| Endpoint | Depth 1 vs 3 | Depth 1 vs 6 | Depth 3 vs 6 |
|---|---:|---:|---:|
| Regression CYP1A2 | .848 | .820 | .829 |
| Regression CYP2C9 | .843 | .817 | .844 |
| Regression CYP2D6 | .807 | .769 | .805 |
| Regression CYP3A4 | .897 | .867 | .887 |
| Classification CYP2D6 | .484 | .388 | .406 |
| Classification CYP3A4 | .709 | .652 | .666 |

Regression depths are related but not redundant; classification D6 rankings
are substantially different. Within-fold Spearman means preserve this pattern:
depth 1 versus 6 is `.756–.870` across regression endpoints, `.385` for
classification D6 and `.648` for classification 3A4. For a descriptive D6
top-call comparison at the historical `49/750` rate, the top 98 of 1,497 OOF
rows overlap by only 27 between depths 1 and 6. Pooled OOF top-call overlap
does not reproduce deployment's fold averaging or establish which calls would
improve the blinded score.

## Fixed blend follow-up

Recorded run `20260925T070604392889Z-6f643137` reads the saved, keyed OOF
predictions and evaluates two fixed arithmetic means: depths 1+6 at 50:50 and
depths 1+3+6 at one third each. It fits no blend weights and uses no blinded
labels. Because this follow-up was motivated by the earlier depth screen, its
numbers are exploratory rather than an independent evaluation of a selected
recipe.

| Endpoint / metric | Final depth 6 | Depths 1+6 | Depths 1+3+6 |
|---|---:|---:|---:|
| Regression 1A2 ST-RAE | .83912 | .82422 | **.81913** |
| Regression 2C9 ST-RAE | .73858 | .72763 | **.71969** |
| Regression 2D6 ST-RAE | .91246 | **.89681** | .90188 |
| Regression 3A4 ST-RAE | .62548 | .59976 | **.59581** |
| Regression macro ST-RAE | .77891 | .76210 | **.75913** |
| Classification 2D6 AUROC / AP | .53712 / .25133 | **.54311** / .27385 | .53516 / **.27637** |
| Classification 3A4 AUROC / AP | .72705 / .42105 | **.74908** / .44022 | .74784 / **.44356** |

The fixed depth blends improve regression ST-RAE on every endpoint and improve
classification AP on both endpoints versus the final depth alone. This is
evidence of complementary OOF errors within these simple CheMeleon readouts.
It does not measure incremental value when added to the current multi-family
submitted ensembles, nor does it establish blinded-test benefit. No candidate
or submission was changed. Compact [blend metrics](../../artifacts/tracked/chemeleon-depth-probe-20260925/blend-metrics.csv)
and [keyed blend OOF](../../artifacts/tracked/chemeleon-depth-probe-20260925/oof-blends.parquet)
are tracked with the standalone results.

## Interpretation and limits

The result supports the possibility of task-dependent transfer from earlier
message-passing states. The CYP1A2/2C9/3A4 regression gains are small and were
observed while looking across three depths and six endpoints. There is no
independent split confirmation or multiplicity-adjusted uncertainty estimate.
The frozen `W_o` was trained on the final state, which makes intermediate
extraction a practical probe but not a fully optimized intermediate readout.
Pooling, feature dimension and downstream learner were held fixed. A next
focused test, if warranted, would check whether this depth blend adds residual
value to the current submitted ensembles on matched held-out folds.

The compact [metrics](../../artifacts/tracked/chemeleon-depth-probe-20260925/metrics.csv)
and [keyed OOF predictions](../../artifacts/tracked/chemeleon-depth-probe-20260925/oof.parquet)
are tracked. Full feature tables (17/26/26 MB for depths 1/3/6) remain in the
immutable extraction run with their checksums and provenance.
