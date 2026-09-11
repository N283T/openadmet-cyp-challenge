# Corrected parent-retained CV with the existing full feature recipe

## Objective and correction

The earlier series-held-out experiment removed the parent hit and used only
Mordred descriptors. It therefore did not execute the requested comparison of
known-hit expansion with the existing feature pipeline. This experiment retains
parent hits, rebuilds the supervised features for the new partition, and uses
LightGBM for the trial-stage readout as requested.

## Membership

Reuse the fixed threshold-0.40 structure groups and seed selection recorded by
`20260911T100134004950Z-c7a3f2bf`; do not tune groups on the new model scores.
Every selected seed's entire InChIKey connectivity block is always in direct
readout training. Non-anchor compounds keep their shared series-fold membership.
For outer fold f, validation is `fold == f and not anchor`. Compare:

- Parent retained: train on all remaining labeled molecules.
- Parent excluded: train on `fold != f`, removing parents in the held group too.

Validation molecules are exactly identical in both arms. Neighbor-only evaluation
uses the selected hit groups after removing all anchor connectivity blocks;
all-nonanchor evaluation is also retained. Shared parent membership is applied
across all targets, with available labels used per endpoint.

Membership run `20260911T104703655337Z-a3923255` writes parent-assignments.parquet,
membership.parquet, coverage.csv, summary.json. It verifies that every anchor is
always in training and that no connectivity block has mixed roles in a fold.
Target anchor/neighbor counts are 40/78, 30/46, 36/68, 35/107 for 1A2/2C9/2D6/3A4.

## Features and controlled readouts

Use frozen CheMeleon 2048-dimensional embeddings, endpoint/outer-training PCA256,
four nested predicted primary log2fc heads, and Mordred descriptors. Drop only
training-empty/constant descriptor columns and median-impute from training rows,
matching the original descriptor preprocessing. PCA and descriptor preprocessing
are refitted for each parent-retention arm. Also evaluate the Mordred-only ablation
on exactly the same validation molecules and arms.

Refit ten CheMeleon auxiliary models excluding each unordered pair of the original
series folds, with the original learning rate, maximum epochs, early stopping and
initial checkpoint. Auxiliary training excludes complete groups, including anchor
molecules in those groups. For recipient group g under outer fold f, use the model
excluding f,g; for any recipient in group f, including a retained parent, average
the four models excluding f. Thus a retained parent's auxiliary feature does not
use its own primary target. This conservative auxiliary exclusion is common to
both direct-readout arms; only the readout gains the retained parents.

The auxiliary loader now honors `outer_strategy` (default UMAP remains unchanged).
Per-parent membership tables verify exclusion from optimization and inner model
selection. No old UMAP-supervised predictions are reused.

LightGBM settings are unchanged from the full-feature point-label baseline:
L1, 300 trees, 15 leaves, learning rate .05, min child samples20, L2=1, seed42+fold,
four CPU threads. No hyperparameter search, external targets or classification.
The readout comparison comprises 80 fits (2 retention arms × 2 feature sets ×
4 endpoints × 5 folds), in addition to ten auxiliary models.

## Planned result checks

Audit exact compound/endpoint/fold/label correspondence and equal validation keys
between arms. Report ordinary MAE, ST-RAE numerator/denominator, bias and ranking
correlation beside ST-RAE. Group-bootstrap paired retention deltas are conditional
on the fitted models. Report results even if retaining parents fails to improve.

## Results

All ten auxiliary models and 80 LightGBM readouts completed. The audit verifies
identical compound/endpoint/fold/label keys across the four readout arms. It stores
576 endpoint/fold/metric records in the DB.

| Features | Parent in training | Neighbor ST-RAE | Neighbor MAE | Neighbor Spearman |
|---|---|---:|---:|---:|
| Mordred | Excluded from held group | 0.950438 | 0.604905 | 0.382179 |
| Mordred | Retained | 0.934595 | 0.597441 | 0.398397 |
| Full original recipe | Excluded from held group | 0.873095 | 0.569082 | 0.465518 |
| Full original recipe | Retained | **0.836776** | **0.547848** | **0.510554** |

The validation population is the same in all four rows. Restoring the feature
recipe accounts for 0.077344 ST-RAE improvement under excluded parents; retaining
parents then accounts for another 0.036318. The earlier high score was not a fair
readout of the original model pipeline under the requested parent-retained design.
The Mordred/excluded result reproduces the preceding series0.40 neighbor score.

| Endpoint | n | Full excluded ST-RAE | Full retained ST-RAE | Full retained MAE |
|---|---:|---:|---:|---:|
| 1A2 | 78 | 0.885987 | 0.866457 | 0.687147 |
| 2C9 | 46 | 0.911243 | 0.901990 | 0.456718 |
| 2D6 | 68 | 1.059278 | 0.982010 | 0.538873 |
| 3A4 | 107 | 0.635870 | 0.596648 | 0.508652 |

For full-feature models, group-bootstrap intervals for the retained-minus-excluded
ST-RAE delta are [-.0618,.0236], [-.1278,.0698], [-.1369,-.0430],
[-.0795,-.0064] respectively. Thus the paired evidence is clearest for 2D6 and3A4.
For 2D6, Spearman rises .2230→.3952 and MAE falls .5692→.5389; mean prediction
bias remains around −.25 on this potency-enriched neighbor subset.

On all 6,384 non-anchor endpoint labels, full excluded/retained macro ST-RAE is
0.691585/0.687660, MAE 0.539483/0.537846. This is a different population and split
from the earlier full-UMAP/ensemble numbers and is recorded as a separate result.
The neighbor-only score remains the designated parent-expansion diagnostic.

## Decision

Use the corrected parent-retained membership with the full feature recipe for
known-hit expansion comparisons. Keep the paired excluded arm and UMAP reference
as named controls. Trial-stage readouts stay LightGBM as requested; no TabPFN fits
were added here. The remaining endpoint errors are now measured under the intended
comparison, not attributed to a different descriptor-only pipeline.

## Completed run lineage

Auxiliary fold-pair runs, in 01,02,03,04,12,13,14,23,24,34 order:
- `20260911T104624805993Z-10bc2d72`
- `20260911T104741286637Z-8500c868`
- `20260911T104849271220Z-a4988fac`
- `20260911T105008532833Z-3ed4f07c`
- `20260911T105118731416Z-b4b410b6`
- `20260911T105230072207Z-ae6b6ea0`
- `20260911T105345455281Z-367a21a9`
- `20260911T105454597743Z-2e73e197`
- `20260911T105606770114Z-f21bad3e`
- `20260911T105719482758Z-5d40310e`

Readout folds0–4:
- `20260911T105206572152Z-e998787e`
- `20260911T105552018025Z-562f006a`
- `20260911T105818167518Z-a178d271`
- `20260911T105928961261Z-66728975`
- `20260911T110021961231Z-8ec1d5bb`

Final audit: `20260911T110118749415Z-5167f1a0`. Outputs: keyed oof.parquet,
metrics.csv, macro.csv, paired-bootstrap.csv, summary.json. Readouts persist PCA,
imputers, selected descriptor columns, LightGBM boosters and fit-membership counts.
Checks and 44 DB-enabled tests pass. No remote submission or test-label use.
