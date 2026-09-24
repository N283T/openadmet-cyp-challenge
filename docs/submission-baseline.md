# Baseline v1 submission preparation

Both tracks use the same model family at the user's request: TabPFN 2.6 with
PCA256 of frozen foundation CheMeleon, four nested predicted log2fc values, and
Mordred. Do not substitute LightGBM for classification despite its higher CV MCC.

## Prediction provenance

Inference run `20260910T132141471698Z-343eb5ec` averages the five existing UMAP
readout models from RL-0008. This is a fold ensemble, not a new all-data refit.
Every fold uses its saved endpoint PCA, descriptor column mask, median imputer
and fitted TabPFN state. Test auxiliary features average the four matching
pair-excluded auxiliary models, exactly as in that fold's outer validation.

All inputs are snapshotted, including base weights and saved states. When loading
a `.tabpfn_fit` archive, only a scratch copy's base-weight path is remapped to the
current run input snapshot. All 30 models reproduced their saved full-fold OOF
predictions with maximum absolute difference **0.0** before test inference.
The original fitted states and completed runs remain unchanged.

Outputs include keyed predictions for each fold, mean probabilities/continuous
predictions, regression Parquet and classification Parquet at threshold 0.5.
The regression file passes the 750-compound submission contract. At 0.5, CYP2D6
TDI has zero positive predictions and CYP3A4 has 90. The all-negative CYP2D6 column
violates the official submission contract, so that classification file is not
eligible for candidate creation. A threshold decision is required separately;
do not alter the completed run or silently change its CV interpretation.

Regression candidate: `20260910T132850121465Z-7bbef880`.
Prepared regression submission: `20260910T132854718239Z-4766cd9f`.
File SHA-256: `07debb137f88ae66e9e50be043d0cc01d415ecaf0fe39e32e5e832743c679270`.

## Private profile and automated upload

`.private/submission-profile.json` contains the participant information reused
from the prior PXR submission, as authorized by the user. It is Git-ignored and
mode 0600. The report URL (`model_tag`) is empty. Publication opt-in and proprietary
data declarations retain the prior values. `open_code_checkbox` is false because
this repository is private. A blank example is tracked in
`configs/submission-profile.example.json`; never put the actual profile in configs.

```sh
# Validate the prepared file and private profile; do not contact the submit endpoint.
pixi run cyp submit SUBMISSION_ID --dry-run

# Send exactly once and persist the official response.
pixi run cyp submit SUBMISSION_ID
```

The default profile path can be replaced with `--profile PATH`. Sending acquires a
per-track database lock, snapshots the payload and profile in an ignored 0700
attempt directory, and commits an `unknown` attempt before any network operation.
It compares the live endpoint field names with the expected contract, then makes
one `/submit_predictions` call. Profiles are not echoed in CLI output. The
attempt's profile copy is mode 0600 for audit/recovery.

A recognized receipt is `sent`, not `accepted`: backend scoring has not been
confirmed. Explicit API validation/cooldown errors are `rejected`; exceptions or
unrecognized responses remain `unknown`. Raw responses and outcome records live
under `artifacts/submission-attempts/<track>/<id>/`. An uncertain payload cannot
be resubmitted until its delivery has been independently reconciled. A sent
payload can be deliberately submitted again after the competition cooldown
using a new prepared submission ID and `--resubmit-sent`; reusing the same
submission ID is blocked. This supports returning to a previously scored vector
without changing its prediction file. Rejected payloads require an explicit
`--retry-rejected` after inspecting the saved receipt. There is no automatic
retry or scheduled retry.

`cyp submission-event` still records subsequent scoring evidence. It does not
clear uncertain attempt guards: independently reconcile remote delivery before
any manual database correction/recovery. Applied migration 011 records attempts;
submission events retain the existing append-only evidence history.

Infrastructure tests mock the network and verify that intent is visible from a
second DB connection before sending, dry runs never send, all three outcomes
persist, and an implicit second upload is blocked. They use an isolated DB schema
and fictitious profiles. No production upload is performed by tests.

## User-authorized classification threshold calibration

The user authorized selecting each endpoint's maximum-MCC threshold using the
labeled OOF predictions only. Run `20260910T133042436909Z-c95bbbfc` preserves the
complete threshold search and a separate classification output. Tie-breaking is
closest to 0.5, then the higher threshold. No test outcomes are available or used.

| Endpoint | Threshold | Selection-only OOF MCC | Test positives / 750 |
|---|---:|---:|---:|
| CYP2D6 TDI | 0.13597816228866577 | 0.086555 | 686 |
| CYP3A4 TDI | 0.3193783760070801 | 0.383103 | 249 |

These MCC values are optimization results on the same OOF labels, **not** unbiased
CV scores, and do not replace RL-0008's fixed-threshold metrics. Deployment also
averages five models whereas the selection OOF has one holdout model per row.
CYP2D6's 686 positive predictions reflect a weak classifier/unstable threshold;
passing the submission validator does not establish predictive quality.

Classification candidate: `20260910T133115130724Z-48208365`.
Prepared classification submission: `20260910T133122668124Z-68969a86`.
Both tracks pass dry-run validation with the ignored private profile and empty
report URL. The named baseline manifest is `configs/submission/baseline-v1.json`.
The original 0.5-threshold predictions are retained and not submitted.

## Submission receipts

Both tracks were sent once on 2026-09-10 UTC and the official Space returned
`Submission received` for N283T. These receipts confirm intake, not completed
scoring; the response says leaderboard processing should complete within two hours.

| Track | Attempt | Local state |
|---|---|---|
| Regression | `20260910T133245035134Z-aeaa96e1` | sent |
| Classification | `20260910T133255631616Z-f2dbb957` | sent |

Raw replies, source snapshots and durable outcomes are under the corresponding
`artifacts/submission-attempts/<track>/<attempt>/` directories. Private profile
snapshots remain ignored and access-restricted. Neither attempt was retried.
The classification file SHA-256 is
`fdddaea138eb27dba4a36d10a17a557a4cc579b4fdacef5f8e08ec4485ced312`.
Infrastructure check and 37 tests passed, including mocked network/durable-state
integration tests. Scoring/leaderboard results have not yet been recorded.

## First leaderboard result — fetched 2026-09-10 13:49–13:50 UTC

Both rows identify N283T and the 2026-09-10 13:32 UTC submission, matching the
baseline receipts. Submission events now record accepted/scored evidence. Rankings
are those in this snapshot and can change as other submissions are scored.

| Track | Rank | Primary metric | Other aggregate metric |
|---|---:|---:|---:|
| Regression | 95 | MA-ST-RAE 0.868726 | MA-MAE 0.990925 |
| Classification | 61 | MA-MCC 0.216149 | MA-Accuracy 0.444873 |

| Regression endpoint | UMAP ST-RAE | LB ST-RAE | LB rank | UMAP MAE | LB MAE |
|---|---:|---:|---:|---:|---:|
| CYP1A2 | 0.816778 | 0.785188 | 90 | 0.651932 | 1.017915 |
| CYP2C9 | 0.625471 | 0.588810 | 91 | 0.450591 | 0.558855 |
| CYP2D6 | 0.928701 | 1.416195 | 92 | 0.599081 | 1.723898 |
| CYP3A4 | 0.527865 | 0.684710 | 110 | 0.536985 | 0.663030 |

| Classification endpoint | LB MCC | LB rank | LB precision | LB recall |
|---|---:|---:|---:|---:|
| CYP2D6 | 0.001644 | 75 | 0.070760 | 0.923001 |
| CYP3A4 | 0.430653 | 18 | 0.562266 | 0.662249 |

CYP2D6 is the main weakness. Its regression Spearman is 0.373798 on LB versus
0.360647 in UMAP CV, despite the much larger absolute error. This motivates a
check of value ranges, calibration and population differences; aggregate scores
do not identify a specific cause or expose the hidden labels. The observed
ST-RAE gap does not by itself establish that the UMAP split should be replaced.
CYP1A2/CYP2C9 have lower LB ST-RAE despite higher MAE, illustrating why different
populations/credible intervals and metric normalization complicate comparisons.

Classification is not directly comparable with the original 0.5-threshold CV
MCC: deployed thresholds were fitted on OOF labels. CYP2D6's 686/750 positives
and near-zero LB MCC show that threshold tuning did not solve its weak ranking.
CYP3A4 is materially stronger. Do not claim either the selection-only OOF MCC or
this first LB result as evidence for general model superiority. No retraining,
new threshold search or resubmission was done while fetching these results.

Raw CSV snapshots and live API layouts are stored under `artifacts/api/<track>/`:

| Track / endpoint | Snapshot ID |
|---|---|
| Regression MA | `20260910T134937823519Z-cc971b38` |
| Classification MA | `20260910T134938608852Z-2a9f25aa` |
| Regression CYP1A2 | `20260910T135018466919Z-66f7e785` |
| Regression CYP2C9 | `20260910T135018468641Z-dfc27da3` |
| Regression CYP2D6 | `20260910T135018466919Z-88b44c6e` |
| Regression CYP3A4 | `20260910T135018484029Z-185aaa76` |
| Classification CYP2D6 | `20260910T135018498829Z-6408f824` |
| Classification CYP3A4 | `20260910T135018507465Z-d8cffd06` |
