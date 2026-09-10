# Foundation PCA + predicted log2fc with TabPFN

## Protocol

The requested experiment uses 256 PCA components of the unmodified CheMeleon
message-passing representation and four predicted primary-screen log2fc values,
with TabPFN-2.6 as the readout. The official preview uses TabICL; this experiment
keeps our existing TabPFN version, eight estimators and seed 42 + outer fold.
Classification remains at a fixed 0.5 threshold. There is no Mordred concatenation.

Reference inspected: OpenADMET/openadmet-models branch
`preview/chemeleon-pca-tabicl`, commit `4a8df806d12a7385edb6365074ce218af19002a5`.
The user's supplied baseline description specifies the 256 + 4 construction.
Source checkout lives at `/home/nagaet/ghq/github.com/OpenADMET/openadmet-models`.
The preview provides reusable embedding, trained-model prediction and PCA
components; it does not automatically establish our nested auxiliary CV protocol.
We use project-native scripts and the existing pixi environment.

## Auxiliary cross-fitting

Use the existing UMAP five folds. For every unordered pair of folds, train one
four-head auxiliary model excluding both folds from optimizer updates, target
scaling and early-stopping selection. This requires ten models. Hyperparameters
match the first auxiliary campaign; the inner 90/10 split is only for stopping.
Membership records mark both excluded folds as `outer_holdout`; the config records
both fold IDs. Auxiliary `outer_predictions.parquet` is diagnostic for the combined
excluded set, so its legacy single `fold` field must not be used as UMAP membership.
Actual UMAP membership always comes from the original assignment input.

For outer fold f and a readout training molecule in fold g, use the prediction
from the model excluding f and g. For validation molecules in f, average the four
models excluding f paired with each other fold. Thus every auxiliary prediction
used for learning or evaluation excludes the recipient's fold, and all four
models exclude the outer validation fold. Training features use one auxiliary
model while validation features use a four-model average, a remaining distribution
difference. Models use three of five folds before their early-stopping split;
this is less auxiliary data than a final all-data fit.

The original five single-exclusion auxiliary OOF tables are not concatenated into
the readout training data. Doing so would admit the outer holdout into auxiliary
models used by other training folds.

## PCA and outputs

Use existing raw foundation run `20260910T012630927160Z-e9576045`. For each endpoint
and outer fold, fit PCA on that endpoint's labeled training compounds only. Center
without external standardization or whitening, retain 256 components, then append
the four auxiliary predictions in the recorded head order. Persist PCA state,
explained variance, keyed 260-dimensional features, TabPFN fit states, constructor
parameters, checkpoint hashes and OOF predictions. The two foundation/readout
checkpoint provenance limitations are unchanged from the earlier reports.

Launch auxiliary fits with `scripts/launch_stacking_aux.sh`; after completion,
`scripts/launch_stacking_readout.py` builds immutable per-fold input maps and
launches readouts. Every fit is a `cyp run` with snapshotted source/config/inputs.
Configs are not overwritten on a repeat invocation. Run replay creates new runs.
`--fold` selects individual readout folds once their four auxiliary parents are
complete. The inherited `encoder_run` configuration field is unused by this
readout; actual feature lineage is `foundation_run` plus the four `auxiliary`
entries, also checked by the DB loader. Auxiliary models remain registered as
ordinary experiment runs, with explicit many-to-many parent links from readouts.

Comparisons against the prior fine-tuned-embedding readout change both the
representation and auxiliary prediction construction. They do not isolate PCA,
predicted-log2fc benefit or the effect of pretraining individually. No blinded-test
predictions or submission are part of this experiment.

## Results — 2026-09-10 UTC

All ten auxiliary fits and five readout runs completed. No failed model runs.
The audit `20260910T121725442380Z-414b10ed` verifies 11,606 unique OOF
compound/endpoint pairs, exact source-label coverage and official ST-RAE agreement
at 1e-12 tolerance. Database load `20260910T121726644068Z-21bcc734` verifies
parent input hashes and prediction fold membership before registering all outputs.
Infrastructure checks and 30 tests passed, including exclusion of both auxiliary
holdout folds from training, selection and target scaling.

| Pipeline (same UMAP split) | Macro ST-RAE ↓ | MAE ↓ | MCC ↑ | AUROC ↑ |
|---|---:|---:|---:|---:|
| LightGBM + Mordred + Morgan | 0.753365 | 0.580243 | 0.165203 | 0.647892 |
| Fine-tuned embedding 2048 + TabPFN | 0.731115 | 0.562771 | 0.112674 | 0.661075 |
| Foundation PCA256 + predicted log2fc4 + TabPFN | 0.729505 | 0.563181 | 0.124720 | 0.645634 |

These are pooled OOF endpoint metrics followed by equal-weight endpoint averaging,
not averages across folds. ST-RAE improves by only about 0.0016 versus the previous
TabPFN experiment, while MAE slightly worsens. Classification MCC increases slightly
but remains below the LightGBM baseline. There is no repeated-seed uncertainty
estimate and no evidence of a substantial improvement.

| Regression endpoint | ST-RAE | MAE | Spearman |
|---|---:|---:|---:|
| CYP1A2 | 0.810166 | 0.650819 | 0.506643 |
| CYP2C9 | 0.640403 | 0.454253 | 0.641500 |
| CYP2D6 | 0.929883 | 0.599993 | 0.354759 |
| CYP3A4 | 0.537568 | 0.547657 | 0.756292 |

| TDI endpoint | MCC | AUROC | AP | Predicted positive rate |
|---|---:|---:|---:|---:|
| CYP2D6 | 0.000000 | 0.525112 | 0.228490 | 0.00% |
| CYP3A4 | 0.249440 | 0.766157 | 0.466774 | 8.65% |

CYP2D6 TDI remains all-negative at 0.5 and has nearly chance-level ranking, so a
threshold change alone is not an established solution. Retain this experiment as
a compact alternative baseline and negative evidence against expecting a large
benefit from this feature recipe alone. Descriptor augmentation and a PCA-only
ablation are possible next experiments; neither was run here.

### Run inventory

| Excluded auxiliary folds | Run |
|---|---|
| 0, 1 | `20260910T120234802703Z-540ba590` |
| 0, 2 | `20260910T120346309537Z-45982a57` |
| 0, 3 | `20260910T120459411341Z-20bcd292` |
| 0, 4 | `20260910T120615652505Z-721bdc87` |
| 1, 2 | `20260910T120723325391Z-8b2cbaa4` |
| 1, 3 | `20260910T120827257407Z-89114984` |
| 1, 4 | `20260910T120935716146Z-38348e57` |
| 2, 3 | `20260910T121105542801Z-b3d5494f` |
| 2, 4 | `20260910T121240644581Z-370232da` |
| 3, 4 | `20260910T121356651841Z-97532ecb` |

| Readout outer fold | Run |
|---|---|
| 0 | `20260910T121031603901Z-fbc2b2ea` |
| 1 | `20260910T121149831493Z-a62f51ef` |
| 2 | `20260910T121402665787Z-08f7afb5` |
| 3 | `20260910T121548828209Z-abe37170` |
| 4 | `20260910T121624900089Z-271684c2` |
