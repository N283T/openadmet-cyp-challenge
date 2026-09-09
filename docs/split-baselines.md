# Split probe and transfer-learning plan

The previous PXR recipe is documented in
`/home/nagaet/openadmet-pxr-model-report/MODEL_REPORT.md` and its split implementation
in `/home/nagaet/pxr-iduction-challenge/track1_activity/src/splits.py`.
Its central features were predicted auxiliary log2fc values plus frozen embeddings
from log2fc-pretrained encoders, with TabPFN readout. This is the intended next
stage for CYP; this probe only fits LightGBM on molecular descriptors.

Results and decision: [2026-09-09 comparison](analysis/2026-09-09-split-baselines.md).
Current decision: [RL-0003](research/2026-09.md#rl-0003--umap-only-and-local-research-log).
`configs/splits/default.json` selects UMAP only for routine evaluation. Reconsider
only if observed LB results diverge substantially from CV; random/scaffold remain
historical comparison candidates.

## Frozen split candidates

Split run: `20260909T135130268493Z-edba6fc0`.
Configuration: `configs/splits/baseline.json`.
All non-test standardized compounds share assignments across endpoints and tracks,
including auxiliary-only compounds. An endpoint model uses only its observed labels.
Test structures are excluded from UMAP and KMeans fitting.

- Random: shuffled molecule-level KFold, seed 42, five folds.
- Scaffold: nonchiral RDKit Murcko scaffold groups, largest-first greedy fold balance.
  All acyclic molecules belong to one empty-scaffold group.
- UMAP: Morgan radius 2 / 2,048 bits / no chirality; Jaccard UMAP with 10 components,
  30 neighbors, seed 42, one thread; KMeans 50 clusters, n_init 10, seed 42;
  whole clusters assigned largest-first to five folds. This matches the previous
  PXR UMAP recipe, on this dataset and the current locked library versions.

These are unsupervised structure-based split candidates, not tuned using scores or
blinded labels. Group balance uses the complete non-test compound universe;
individual endpoints can have uneven fold size or class prevalence.

## Baseline protocol

Configuration: `configs/splits/lgbm.json`; executable: `scripts/split_baselines.py`.
Mordred's 1,613 2D descriptors versus those descriptors plus the Morgan bits above.
No global imputation, scaling, feature selection, early stopping or parameter search.
Nonfinite descriptors become missing; LightGBM handles missing values natively.
300 trees, learning rate 0.05, 15 leaves, minimum 20 samples per leaf, L2 1,
seed 42, deterministic CPU, four threads. Regression uses absolute-error objective;
TDI uses binary objective and a fixed 0.5 probability threshold.

Four regression and two TDI endpoints × three splits × two feature sets × five
folds = 180 models. Each fold model and its keyed OOF predictions are saved.
OOF predictions identify compound ID, molecule name, endpoint, split run,
strategy, fold and feature variant. No blinded-test predictions or submissions.

Regression ST-RAE uses the official credible intervals: the distance outside the
interval divided by that same error for a constant evaluation-set mean. No extra
weight is applied at pIC50 4. Pooled OOF ST-RAE is calculated per endpoint, then
macro-averaged across the four endpoints. MCC is pooled per endpoint and then
macro-averaged across the two TDI endpoints. Fold metrics are retained separately.
This is a single CV estimate, not the leaderboard's 1,000-bootstrap evaluation.

The official scorer snapshot is pinned to tutorial commit
`858ae63ce79934113bccdb7fc65467de5f7b1935` in `artifacts/metric-reference/`.
The recorded audit compares all fold and pooled ST-RAE values against that source,
recomputes MCC, and checks OOF coverage, labels and uniqueness.

Structural diagnostics use maximum Morgan Tanimoto similarity to labeled training
compounds, separately for each endpoint. Validation-to-training distances are
compared with blinded-test-to-all-labeled distances. CV reference pools have 80% of the endpoint labels, while the test comparison uses
100%, which also affects nearest-neighbor similarity. This is a covariate diagnostic;
it cannot establish which CV split predicts leaderboard performance.

## Auxiliary learning boundary

For each outer fold, exclude every held-out compound's auxiliary labels from
supervised encoder/log2fc training, across all CYP endpoints. Freeze the encoder
before its readout evaluation. Never append measured validation/test log2fc.
If predicted log2fc is used as a downstream input, generate inner out-of-fold
predictions for downstream training compounds, and outer-held-out predictions for
validation compounds, to avoid an in-sample versus out-of-sample feature mismatch.
Keep these inner folds and encoder lineage recorded separately from the outer split.

First compare Mordred alone, predicted log2fc alone, Mordred + predicted log2fc,
and frozen embedding + predicted log2fc with a common outer split. TabPFN and
encoder training are not implemented or started by this baseline probe.

## Registry

Migration `005_split_baselines.sql` registers immutable split assignments and
baseline OOF predictions. `scripts/register_split_baselines.py --baseline-run ID --audit-run AUDIT_ID`
verifies completed-run hashes and parent input hashes before an atomic append.
Recomputed `metrics` rows belong to the audit run and use split ID `SPLIT_RUN:strategy`, metric name
`features:metric`, and fold `-1` for pooled OOF. Repeated loading does not overwrite.
