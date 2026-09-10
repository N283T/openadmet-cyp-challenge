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
under `artifacts/submission-attempts/<track>/<id>/`. Sent or uncertain payloads
cannot be submitted again through the command, even via another prepared ID.
Rejected payloads require an explicit `--retry-rejected` after inspecting the
saved receipt. There is no automatic retry or scheduled retry.

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
