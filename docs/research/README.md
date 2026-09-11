# Research log

The canonical research narrative lives here, in English, tracked in this private
Git repository. This carries forward the dated research-log format previously
used in PXR issues #100 and #208, without requiring an Issue for each experiment.

## Current direction

- Focus on regression; defer classification tuning until regression stalls or
  priorities change. Saved-prediction diagnosis is complete; prioritize
  source-preserving external-data curation using the existing local ChEMBL 36
  mirror. See [RL-0012](2026-09.md#rl-0012--regression-diagnosis-and-existing-chembl-inventory).
- Pooled OOF means are nearly calibrated, but CYP2D6 ranking remains weak across
  folds. Hidden-test calibration is unresolved; retain UMAP and separate future
  calibration experiments from the external-data comparison.

- Baseline-v1 uses the augmented TabPFN 2.6 fold ensemble for **both** tracks by
  user decision. First LB scores are recorded: regression ST-RAE 0.868726,
  classification MCC 0.216149; CYP2D6 is the main weakness
  ([RL-0011](2026-09.md#rl-0011--first-baseline-leaderboard-scores)).

- Use the existing UMAP five-fold assignment for routine experiments on both tracks.
  No routine random or scaffold evaluation. Reconsider the split only if observed
  leaderboard results diverge substantially from CV; record the evidence first.
- CheMeleon-initialized auxiliary pretraining is complete ([RL-0005](2026-09.md#rl-0005--chemeleon-initialized-log2fc-pretraining)).
  The first fold-matched embedding-only TabPFN readout is complete ([RL-0006](2026-09.md#rl-0006--embedding-only-tabpfn-readout)).
  Regression improves over the initial LightGBM pipelines; classification MCC does
  not. Predicted-log2fc readout inputs require nested cross-fitting. Scratch
  pretraining remains deferred.
- The PCA256 + predicted log2fc4 TabPFN experiment is complete ([RL-0007](2026-09.md#rl-0007--foundation-pca-and-predicted-log2fc-with-tabpfn)).
  Changes versus the embedding readout are marginal; CYP2D6 TDI remains weak.
- Mordred augmentation gives a small regression improvement; the first direct-FT
  recipe underperforms ([RL-0008](2026-09.md#rl-0008--mordred-augmentation-and-direct-chemeleon-ft)).
  Classification remains below the earlier LightGBM MCC.
- Retain TabPFN 2.6 for the augmented feature recipe: v3 with eight or nine
  estimators underperforms on both tracks ([RL-0009](2026-09.md#rl-0009--tabpfn-30-on-regression-and-classification)).
- See [RL-0003](2026-09.md#rl-0003--umap-only-and-local-research-log) for the current decision.

## Recording policy

Append dated entries to `YYYY-MM.md`, with stable sequential IDs (`RL-0001`, etc.).
Use the actual event date/time with an explicit timezone when available. Label
retrospective entries as such; do not invent timestamps. One entry may cover a
coherent batch of runs. Include failed, interrupted and discarded attempts, not
only improvements. Every attempt still goes through `pixi run cyp run`, even before
there is a result worth discussing. Link the log entry in a run note when useful.

Keep each entry short and evidence-backed:

- Question / hypothesis and track(s).
- Method / changes, configuration and all relevant run IDs.
- Results and limitations; link full reports or artifacts instead of copying tables.
- Interpretation, decision (keep/discard/inconclusive), and next step.

Never silently rewrite an earlier conclusion. Append a new entry referencing the
superseded one. Small factual corrections should be explicitly marked. Keep this
index's current direction up to date with links to the decisions supporting it.

PostgreSQL stores structured run, metric, split and prediction records. Run folders
store immutable source/config/input snapshots, models and outputs. This log stores
the reasoning connecting those records. Detailed reports live in `docs/analysis/`;
official source text remains separate in `references/official/`.

GitHub Issues are optional for bounded implementation tasks or actionable bugs;
they are not a second research log. Existing `--issue` run metadata stays optional.
Reports and research entries default to English; conversation language is independent.

## Entries

| Entry | Topic |
|---|---|
| [RL-0001](2026-09.md#rl-0001--preprocessing-and-descriptor-exploration) | Retrospective: preprocessing and descriptor exploration |
| [RL-0002](2026-09.md#rl-0002--initial-split-baseline-comparison) | Retrospective: LightGBM split comparison |
| [RL-0003](2026-09.md#rl-0003--umap-only-and-local-research-log) | UMAP-only evaluation and local research log |
| [RL-0004](2026-09.md#rl-0004--frozen-chemeleon-raw-2048-fingerprints) | Frozen CheMeleon features and checkpoint provenance |
| [RL-0005](2026-09.md#rl-0005--chemeleon-initialized-log2fc-pretraining) | Five-fold auxiliary pretraining from CheMeleon |
| [RL-0006](2026-09.md#rl-0006--embedding-only-tabpfn-readout) | First frozen embedding TabPFN readout |
| [RL-0007](2026-09.md#rl-0007--foundation-pca-and-predicted-log2fc-with-tabpfn) | PCA256 + nested predicted log2fc4 with TabPFN |
| [RL-0008](2026-09.md#rl-0008--mordred-augmentation-and-direct-chemeleon-ft) | Mordred-only, augmentation and direct FT |
| [RL-0009](2026-09.md#rl-0009--tabpfn-30-on-regression-and-classification) | TabPFN 3.0 with eight and nine estimators |
| [RL-0010](2026-09.md#rl-0010--tabpfn-baseline-v1-preparation-and-submission) | Baseline v1, private-profile upload automation and receipts |
| [RL-0011](2026-09.md#rl-0011--first-baseline-leaderboard-scores) | First baseline LB results, endpoint breakdown and CV comparison |
| [RL-0012](2026-09.md#rl-0012--regression-diagnosis-and-existing-chembl-inventory) | Regression diagnosis, classification deferral and local ChEMBL inventory |
