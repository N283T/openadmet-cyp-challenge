# Matched-feature TabICL comparison

## Question and scope

Does replacing TabPFN 2.6 with the official project's TabICL estimator explain the
performance gap? Use identical compound labels, five UMAP folds, fold-specific
PCA256 + nested predicted log2fc4 + imputed Mordred features. Classification and
external target augmentation remain deferred. Test predictions average five folds.

This isolates the learner on our feature pipeline. It is **not a full reproduction
of the public leaderboard baseline**: the official architecture wrapper does not
pin TabICL, and we have not established its complete submitted training recipe.

Official wrapper source commit: `4a8df806d12a7385edb6365074ce218af19002a5`.
Snapshot: `references/official/cyp/tabicl-wrapper/4a8df806/`.
Match its eight estimators, batch size one, automatic AMP/FA3/offloading and mean
prediction. Use seed 42+fold to match existing TabPFN folds, CUDA, four CPU threads.
Median and q0.55 are saved as secondary fixed readouts; they do not select the
primary mean-versus-mean comparison. No hyperparameter search or remote submission.

Resolved package: TabICL 2.2.0, recorded in pixi.lock. Checkpoint:
`tabicl-regressor-v2-20260212.ckpt`, Hugging Face `jingang/TabICL` revision
`4dcd344ece2c00be9e831fdd35bed57b5ad83e19`, SHA-256
`0db9cb538f114e79026bf08f45f41ad8dd7ad2de2aaca9a5ca8cd3bd9748ae7a`.
Acquisition run `20260911T093102177307Z-beddd36f`.
Two prior acquisition attempts failed on an incorrect local hash-helper import:
`20260911T093046008502Z-bd2f9eab`, `20260911T093056151151Z-1d4a92b0`.
All failed records remain.

The inherited launcher purpose/filename includes “interval-aware”, but these fits
use original point pIC50 labels exclusively, as the snapshotted source demonstrates.
No confidence bound is used as a training target. Source/config/input/weight snapshots
are sufficient to refit; no fitted TabICL object is serialized.

## Results

All 20 endpoint/fold fits completed. Mean ST-RAE is worse for TabICL on each endpoint.

| Fixed readout | Macro ST-RAE | MAE | Spearman |
|---|---:|---:|---:|
| TabPFN mean | 0.724704 | 0.559647 | 0.568038 |
| TabICL mean | 0.735148 | 0.567989 | 0.547448 |
| TabPFN median | 0.708706 | 0.552281 | 0.570013 |
| TabICL median | 0.729370 | 0.564036 | 0.544720 |
| TabPFN q0.55 | 0.702744 | 0.555911 | 0.570012 |
| TabICL q0.55 | 0.730685 | 0.571601 | 0.544686 |
| Existing three-member TabPFN q0.55 | 0.691191 | 0.550526 | 0.581232 |

| Endpoint | TabPFN mean ST-RAE | TabICL mean ST-RAE |
|---|---:|---:|
| 1A2 | 0.816778 | 0.819992 |
| 2C9 | 0.625471 | 0.640465 |
| 2D6 | 0.928701 | 0.942700 |
| 3A4 | 0.527865 | 0.537435 |

TabICL mean test 2D6 has mean 4.630671, SD 0.304099, and 11/750 predictions below 4.
The original TabPFN mean is 4.592663 with 10/750 below 4. Thus switching the learner
does not remove the high/narrow test prediction distribution. This observation
alone does not identify the true hidden-label distribution or justify a correction.

Conclusion: no evidence to adopt TabICL on the current feature recipe. This does
not establish that TabICL is universally worse, nor reproduce/explain the official
LB score. Next priority is the official full feature/training recipe and the
selection of screened compounds into fitted-pIC50 labels, including comparison of
log2fc and structures for missing-label compounds. Do not assign missing fits an
inactive label without assay evidence. Retain the existing prospective reference.

## Records

Fold 0–4 runs:
- `20260911T093200583715Z-a7b5212d`
- `20260911T093324497439Z-e8a75360`
- `20260911T093444498345Z-331fc5dc`
- `20260911T093604244693Z-64342e75`
- `20260911T093724505089Z-5311f2b5`

Official-score audit: `20260911T093904047313Z-d7cfa808`; 504 DB metric rows.
Outputs include endpoint/fold metrics, full keyed OOF, fold-average test and test
summary statistics. Exactly 6,525 OOF rows per method; all compound/fold keys align
with the original diagnosis. No new remote submission or local upload candidate.
`pixi run check` and 44 DB-enabled tests passed (16 existing dependency warnings).
