# Independent source audit: SuperCowPowers CYP training recipe

## Version authority

Recorded source snapshot `20260912T103812216128Z-5ed4d83b` pins Workbench HEAD to
`a1d1b46e290e62014ee3f0dfde5667eb5e93282b`. It preserves current CYP scripts,
Chemprop implementation, GitHub tree/commit metadata, and historical revisions of
the union, mt-aux-100 and D6-specialist recipes. Support snapshot
`20260912T104213710169Z-9388c902` preserves historical/current split, task-weight,
standardization and union-data implementations. Trainer snapshots `97b9c270`
(historical) and `a7f2db9b` (current) include the non-Python Chemprop template.

The **18-head union recipe** is directly recoverable at
`4cf07a8cb4366eb29111fb42396336b20b69ada0` (2026-08-27; formatting commit following
its introduction). It has four challenge pIC50, four challenge log2fc, five ChEMBL
potency and five Veith max-response heads. Later scripts add eight challenge
TDI/Emax heads, giving 26. Current HEAD adds another five Veith potency heads,
giving **31**, although a comment still says 26. The blog's model name is not a
reliable version identifier; reproducing its named 18-head entry requires the
historical target list.

## Model recipe checklist

| Item | Source behavior | Reproduction requirement |
|---|---|---|
| Initialization | Scratch D-MPNN unless `--foundation` explicitly selected | Do not silently initialize Jeremy600 or CheMeleon |
| Scratch encoder | BondMessagePassing hidden700, depth5, dropout0.1 | Match dimensions/depth and message-passing dropout |
| Aggregation | Native NormAggregation | Preserve native norm setting |
| Predictor | RegressionFFN, hidden300, two hidden layers, dropout0.1 | Preserve native predictor and model batch normalization |
| Criterion | Native Chemprop MAE by default | Not MSE, interval loss, or BCE conversion of continuous auxiliaries |
| Task scaling | Per-target training-fold normalization, outputs unscaled | Validation labels never fit target scalers |
| Task weights | Primary inverse-label-count weights, normalized to mean1 | Compute on source recipe's selected table |
| Auxiliary weights | 0.3 × mean(primary weights) for log2fc and public heads in union-p30 | Do not interpret 0.3 as 30% gradient mass |
| Loss denominator | Number of observed labels, not weighted label mass | Use native task-weighted MAE reduction |
| Schedule | Native MPNN optimizer/scheduler; warmup2, LR1e-4→1e-3→1e-4 | Do not substitute fixed two-rate AdamW |
| Limits | Up to400 epochs, patience40, batch64,16-mixed | Any reduced budget is an explicitly different experiment |
| Checkpoint | Lowest native validation loss; reload best checkpoint | No last-epoch or pIC50-only ST-RAE checkpoint substitution |
| Deployment | Five fold-model predictions averaged | “100%” means all data participate across ensemble, not one all-row refit |
| Input | SMILES only | No appended descriptors/embeddings or tree readout in the source model |

The eight-head mt-aux-100 member uses four scored pIC50 and four log2fc targets.
The D6 isoform specialist uses six heads: direct pIC50, log2fc, TDI pIC50, direct
Emax, ChEMBL pIC50 and Veith max_response; weights are 1 for direct and 0.3 each
auxiliary. The single specialist uses only D6 pIC50. The published final ensemble
averages two members for non-D6 endpoints and all four for D6 before placement.

The public code and training-image requirements do not pin the exact deployed
checkpoint/package environment. Code-level reconstruction cannot establish
bitwise identity with an author's private saved model. Published LB placement
constants are separately board-derived and are not training hyperparameters.

## Source CV and interpretation

Default CV is five-fold Bemis-Murcko scaffold splitting, seed42. Scaffold IDs
follow first appearance in the ordered molecular table. `group_folds` visits a
seeded random permutation of groups and puts each group into the currently
smallest fold. It does not sort scaffold groups by size and is not sklearn's
GroupKFold. Empty scaffolds share a group. Each model filters to rows with at least
one of its selected targets, so public/auxiliary data change the scaffold population
and thus the challenge compounds' fold assignments.

The held-out fold is both the validation set used for best-checkpoint selection
and the source of its OOF predictions. Thus native published OOF is not the same
protocol as our inner-selected parent outer evaluation. Matching parent folds
while retaining source training settings is a controlled adaptation; a separate
native-scaffold diagnostic is needed for a closer source-CV comparison. Public
holdout rows also belong in the source validation objective, weighted across
selected heads, rather than silently validating only challenge pIC50.

## Data and SMILES requirements

`cyp_aux_features.py` left-joins measured log2fc onto the challenge inhibition
FeatureSet. It does not independently include every row of the broader TDI file.
Historical union18 combines challenge pIC50/log2fc with ChEMBL potency and Veith
max_response, keeping separate assay scales. Veith max_response is clipped to
[-150,50]. Tox21, Veith potency and activity57.14uM are not union18 targets.

The union merges challenge/public overlaps into one connectivity-keyed row.
ChEMBL duplicate connectivity keys retain the first row. Veith first averages
repeated SMILES/isoform observations through a pivot, then retains the first row
per connectivity key. This is not the same as averaging all connectivity-equivalent
stereoisomers or dropping every challenge/public overlap. Challenge identifiers
and SMILES take precedence on shared rows.

The model datapoint helper parses supplied SMILES directly; it does not itself
standardize them. Molecular standardization belongs upstream. Workbench's available
standardizer includes cleanup, fragment/charge handling and default tautomer
canonicalization, which is not automatically identical to this repository's
ChEMBL-parent preprocessing. Exact source training strings and source producer
policies must be checked before calling data reproduction exact.

The reconstruction now resolves these data requirements: public-source download
`20260912T103958592504Z-a85b8d52` obtains the author's actual ChEMBL37 table, and
`20260912T104336148000Z-ae7e37fc` reproduces the historical union's **31,670 rows
and 161,400 labels**. It preserves author SMILES and the 4,905-row inhibition
cohort. Every one of the current 6,525 direct labels and both confidence bounds
matches within 1e-9. The author's blinded structure-only file independently
matches all 750 prediction IDs and SMILES in `5cc21663`. See the companion
[data reconstruction report](2026-09-12-scp-reproduction-data.md).

One public Veith connectivity also occurs in the blinded structure set. It is
retained as public assay data, explicitly recorded, and does not carry a hidden
challenge label. Therefore “all blinded chemistry excluded” would be inaccurate.

## Independent native trainer review

Source inspection of `scripts/scp_repro_train.py` confirms that the adapter imports
the pinned original `train_chemprop_fold`, datapoint helper and scaffold splitter.
It extracts and checks the historical template's hyperparameters, calls the native
trainer and prediction helper, and adds a callback for logging. The minimal
Polars-backed frame supplies the column operations used by these source functions.
It does not reimplement MAE, optimizer, schedule, target normalization or checkpoint
selection. The historical FoldSpec signature is respected. Native best checkpoint
epochs, membership tables, per-head predictions and five-fold means are retained.

The parent-adapted mode maps every challenge connectivity, including external
aliases, to the existing held-out fold and keeps designated anchors in training.
Other external rows receive native scaffold folds, and their validation labels
participate in the native validation objective. This closes the alias exclusion
gap that a compound-ID-only override would leave. Assertions require exactly the
expected canonical scored-label membership in each parent fold.

This is a source-training reconstruction with an explicitly changed split when
`cv_policy=parent_adapted`. Its validation fold also selects the checkpoint, as
the source does; it must not be described as an independently inner-selected OOF.
`source_scaffold` retains the native grouping without parent overrides. Neither
mode establishes bitwise equivalence to a private author checkpoint because the
historical deployed package environment and checkpoint are unavailable. This
review establishes source-level correctness, not completed-run numerical results.

## Material local-score definition difference

Metric-source snapshot `20260912T104448283096Z-79a8113d` confirms the historical
and current `workbench.utils.metrics_utils.soft_threshold_rae` default to
`soft_baseline=False`. The source CYP comparison script calls that default.
Its score is interval absolute error divided by ordinary mean-predictor MAE:

`source_local = sum(distance(pred,[L,U])) / sum(abs(y - mean(y)))`.

Our saved official challenge scorer explicitly soft-thresholds the denominator:

`official = sum(distance(pred,[L,U])) / sum(distance(mean(y),[L,U]))`.

Therefore these local metrics are **not numerically comparable**, even with
identical rows and predictions. For any fixed common evaluation table their ratio
is exactly the ratio of the two denominators; it does not depend on the model.
Report both definitions when comparing source-local OOF claims, and retain the
saved tutorial definition for continuity. The deployed backend's denominator
cannot be proven from the public frontend; the completed trace below records
the available source and chronology evidence.

The source `cyp_compare.py` selects captures in this priority order: target-specific
analog holdout, generic analog holdout, target-specific CV, generic full-cross-fold.
It then joins challenge labels/confidence bounds by molecule ID, drops missing
labels/predictions, and computes Pearson, Spearman, R² and the above local score.
Consequently a quoted comparison from that script may describe the deliberate
analog holdout rather than native scaffold OOF; capture identity must accompany
reported numbers.

### Live platform trace and scorer chronology

Recorded live-source fetch `20260912T104656211077Z-06952a87` pins HF Space
`openadmet/cyp-challenge` to `63a1d372076bd1c1c8df2f3a84ed18180fead018` and the
current official tutorial to `858ae63ce79934113bccdb7fc65467de5f7b1935`.
The tutorial scorer is **byte-identical** to our stored reference (SHA-256
`c0ee287c06c0301d4b3bfd74559b90eb1bd85db68b75b7d018ae92280c0e1dd3`). Its file
history records one introduction, commit `60a779a21b43b84f569de6a24ada3bd1f796a90d`
on 2026-08-19, titled “Add eval code from backend.” This supports retaining the
soft-denominator official reference and provides no evidence of a recent tutorial
switch between denominator definitions.

The live HF app itself does not run evaluation. Public `submission.py` documents
that it writes submissions to S3 and an independent scheduled AWS Lambda
`backend/lambda_handler.py` performs scoring. That backend implementation and its
deployed package/version are absent from the public Space tree and requirements.
Therefore the deployment's exact denominator is not independently established by
this public-front-end inspection. The **confirmed issue is the source author's
local helper difference**, not a demonstrated error in our official reference or
a proven mismatch between the reference and current live leaderboard. No private
backend access or test labels were requested.
