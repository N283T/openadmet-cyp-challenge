# Research log

The canonical research narrative lives here, in English, tracked in this private
Git repository. This carries forward the dated research-log format previously
used in PXR issues #100 and #208, without requiring an Issue for each experiment.

## Current direction

- **Comprehensive D6 classification audit completed:** 70 new controlled fits,
  raw-call reconciliation, actual official scores, 25-model probability audits,
  structural coverage and assay-cohort controls. Fewer labels alone do not explain
  the gap. A4 includes1,250 incomplete-potency-pair negatives versus4 forD6;
  original UMAP is structurally harsher than test even with identical neighbor
  pools. No bulk label/identity failure found; no D6-specific noise ceiling
  established. See [RL0056](2026-09.md#rl0056--comprehensive-d6-classification-explanation-audit)
  and [full audit](../analysis/2026-09-12-d6-comprehensive-audit.md).

- **Classification result confirmed:** MA-MCC0.275867 (rank45), D6 MCC0.176698
  improves from0.001644, while A4 MCC0.375036 falls from0.430653. Preserve the
  old A4 baseline for a future endpoint combination. See
  [RL0055](2026-09.md#rl0055--classification-ensemble-leaderboard-result).

- **Classification submitted:** D6 broad mean and A4 balanced LGBM mean from the
  completed round, 43/378 test positives. Official intake confirmed at September
  12 12:27 UTC; scoring pending. See [RL0054](2026-09.md#rl0054--classification-round-submission).

- **Classification round completed during regression cooldown:** 17 individual
  recipes and seven fixed probability means on the original UMAP folds. New
  inner-threshold MCC reaches 0.096228 for D6 and 0.338112 for A4, from different
  means. Keep the historical A4 baseline, whose AUROC remains higher. Models,
  inner/outer predictions and test probabilities saved; no upload. This resumes
  classification after its earlier deferral. See [RL0053](2026-09.md#rl0053--classification-round-during-regression-cooldown)
  and [completed comparison](../analysis/2026-09-12-classification-round.md).

- **Calibration comparison complete:** 30 fixed methods; native macro improves
  0.622574 to 0.601682, but D6 only 0.823983 to 0.819776 with a conditional
  uncertainty interval crossing zero. Source placement is separately applied to
  the exact previously submitted vector; D6-only CSV changes no other endpoint.
  Files validated, no upload. See [RL0051](2026-09.md#rl0051--broad-calibration-screen-and-same-vector-placement-controls)
  and [completed comparison](../analysis/2026-09-12-placement-round.md).

- **Native community learning reproduction complete:** all 25 fits use source
  data and native Chemprop defaults. Fixed mean gives D6 0.823983 / Spearman
  0.472811, essentially tied with previous 0.825436 / 0.472752; macro worsens
  to 0.622574. A separate audit identifies different local ST-RAE denominators,
  not a demonstrated live-LB scoring bug. No adoption or submission. See
  [RL0050](2026-09.md#rl0050--native-community-learning-reproduction-and-metric-definition-audit)
  and [completed comparison](../analysis/2026-09-12-scp-training-reproduction.md).

- **2D6-focused follow-up complete:** specialist scope/seed means, external efficacy
  features, actual joint public graph learning, interval-width decisions and kernel
  regressors. Fixed scope-mean addition gives D6 0.822232 versus 0.825436, but the
  bootstrap interval includes zero; expanded Caruana is slightly worse (0.825952).
  No submission. See [RL0049](2026-09.md#rl0049--cyp2d6-specialists-joint-external-responses-and-interval-decisions)
  and [completed program](../analysis/2026-09-12-d6-focus-round.md).

- **Report-derived round complete:** TDI auxiliaries, Hill-linked response heads,
  geometry/charge proxies, external substrate scores and joint descriptors.
  Expanded 126-candidate weight-crossfit macro is 0.582816 versus 0.583257.
  Small 1A2 gain; 2D6 essentially unchanged. Atom-level SoM remains untested.
  No submission. See [RL0048](2026-09.md#rl0048--report-derived-auxiliary-and-descriptor-experiments)
  and [completed program](../analysis/2026-09-12-report-experiments.md).

- **Global Caruana completed:** 107 saved candidate identities across model families.
  Weight-crossfit macro is 0.583257 versus 0.586931; gains are in 2C9 and 3A4,
  while 1A2 and 2D6 slightly worsen. Full-OOF fitted score 0.573705 is a fitting
  result, not independent CV. Weights and predictions are retained; no submission.
  See [RL0047](2026-09.md#rl0047--global-caruana-across-saved-model-families)
  and [the global comparison](../analysis/2026-09-12-global-caruana.md).

- **Jazzy comparison complete:** standard ETDG features for all 6,896 compounds
  take about two minutes; 12 net-charged structures retain missing values. Paired
  LGBM/TabPFN controls show small mixed gains, with no improvement over the latest
  ensemble. Retain features and results; latest candidate and submission remain
  unchanged. See [RL0046](2026-09.md#rl0046--dg-jazzy-descriptors-and-controlled-readouts)
  and [the controlled comparison](../analysis/2026-09-12-jazzy-readout.md).

- **Latest focused round:** prioritize1A2/2D6, allowing useful auxiliary tasks.
  Improved log2fc features and a complementary2D6 readout yield OOF-only
  candidates around0.5869; no new submission. Latest actual LB is0.721167.
  See [RL0045](2026-09.md#rl0045--focused1a22d6-experiments-and-improved-primary-features)
  and [the focused program](../analysis/2026-09-12-focus-program.md).

- **2026-09-12 update:** The regression five-component ensemble reached development
  OOF ST-RAE **0.596414** on the retained parent CV. Its full-data refit has now
  been submitted without calibration; its scored LB is 0.721167 (RL0044).
  See [RL0043](2026-09.md#rl0043--full-data-five-component-submission) and the
  [deployment report](../analysis/2026-09-12-full-ensemble-submission.md).
  This supersedes the older model/pause/split status bullets below, which retain
  links to the earlier decisions. Regression remains the priority.

- Interval-aware TabPFN/LightGBM and neighborhood trials are complete. The best
  numerical blend scores 0.688853, but its incremental gain over the simple q0.55
  ensemble is uncertain; retain alternatives without replacing the submitted baseline
  ([RL-0017](2026-09.md#rl-0017--interval-aware-readouts-and-neighborhood-alternatives)).
- Continued regression trials favor a simple three-member TabPFN 2.6 q0.55
  ensemble in exploratory OOF (ST-RAE 0.691191). Local q0.55/median candidates
  are saved, with hidden-test population shift unresolved
  ([RL-0016](2026-09.md#rl-0016--distribution-summaries-and-simple-ensembles-improve-regression-cv)).
- External activity modeling is paused. A PXR-inspired calibration screen supports
  modest interval-loss gains but does not identify the hidden test-population
  correction. Keep full/half-strength maps for prospective evaluation; preserve
  the existing submission
  ([RL-0015](2026-09.md#rl-0015--pxr-inspired-regression-calibration-without-external-activity-data)).
- Focus on regression; defer classification tuning until regression stalls or
  priorities change. Saved-prediction diagnosis is complete
  ([RL-0012](2026-09.md#rl-0012--regression-diagnosis-and-existing-chembl-inventory)).
  ChEMBL 36 and PubChem auxiliaries are now curated with source lineage and
  challenge-connectivity exclusion: 24,075 external compounds, 93,119 labels
  ([RL-0013](2026-09.md#rl-0013--external-chembl-and-pubchem-auxiliary-dataset)).
- The first matched neural comparison is complete: external auxiliaries improve
  the point estimate of MA-ST-RAE from 0.824695 to 0.797301, but do not reach the
  adopted TabPFN's 0.724704. Retain the dataset and baseline; do not adopt this
  neural recipe. The macro conditional interval slightly crosses zero
  ([RL-0014](2026-09.md#rl-0014--first-matched-external-data-regression-comparison)).
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
| [RL-0013](2026-09.md#rl-0013--external-chembl-and-pubchem-auxiliary-dataset) | ChEMBL/PubChem acquisition, curation and source lineage |
| [RL-0014](2026-09.md#rl-0014--first-matched-external-data-regression-comparison) | Matched neural external-data comparison and conditional uncertainty |
| [RL-0015](2026-09.md#rl-0015--pxr-inspired-regression-calibration-without-external-activity-data) | PXR-inspired interval calibration, label-isolation checks and unresolved test shift |
| [RL-0016](2026-09.md#rl-0016--distribution-summaries-and-simple-ensembles-improve-regression-cv) | Distribution summaries, simple ensembles and local prospective candidates |
| [RL-0017](2026-09.md#rl-0017--interval-aware-readouts-and-neighborhood-alternatives) | 120 interval/point readout fits, neighborhood alternatives and small incremental gains |
