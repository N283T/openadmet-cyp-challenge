# CYP Space API — verified 2026-09-10

Sources:
- https://openadmet-cyp-challenge.hf.space/gradio_api/info
- https://huggingface.co/spaces/openadmet/cyp-challenge/blob/main/config.py
- https://huggingface.co/spaces/openadmet/cyp-challenge/blob/main/submission.py

Live inspection captures the Space revision and revision-pinned copies under
`artifacts/api/<id>/`; public source is treated as data, never imported/executed.

## Submission

`/submit_predictions` takes these named parameters:

```
username, user_alias, anon_checkbox, participant_name, discord_username,
email, affiliation, model_tag, paper_checkbox, proprietary_data_checkbox,
open_code_checkbox, track_select, file_input
```

`track_select` currently accepts `Regression Prediction` or `Classification Prediction`.
The response is a status string, not proof of completed backend scoring. Submission
uses a separate Git-ignored profile. `cyp submit SUBMISSION_ID` now automates a
single upload with durable attempt/receipt records and duplicate guards. See
[baseline submission preparation](submission-baseline.md) for the private profile,
dry-run command, delivery states and recovery policy.

Tracks have independent 12-hour submission cooldowns in the public source. Both
require exactly 750 compounds and identifiers `SMILES` and `Molecule_Name`.
CSV and Parquet are accepted; Parquet is preferred by the official instructions.

Regression columns:
- `CYP1A2_pIC50_direct_inhibition`
- `CYP2C9_pIC50_direct_inhibition`
- `CYP2D6_pIC50_direct_inhibition`
- `CYP3A4_pIC50_direct_inhibition`

Values must be finite, nonconstant and have sample standard deviation >= 0.01.
Primary metric is macro soft-threshold RAE (MA-ST-RAE).

Classification columns:
- `CYP2D6_is_TDI`
- `CYP3A4_is_TDI`

Use boolean predictions, with both classes present per column. The UI validator
also accepts binary 0/1, but the local contract uses booleans to match the stated
file schema. Primary metric is MCC. Probabilities are retained as experiment
outputs and converted using a recorded threshold before candidate creation.

Structure is currently disabled in submission options; do not use its provisional
schema or count for infrastructure decisions.

## Leaderboard artifacts

The API exposes separate Live and Interim leaderboard tabs. Generated `/partial`
names are not stable semantic identifiers. The client resolves the selected board,
track and endpoint from the current component/tab layout and fails if resolution is
ambiguous. It saves the live contract, returned response and actual CSV under a fresh
directory via Gradio's `download_files` parameter.

Use `cyp leaderboard --track regression` or `--track classification`; `--endpoint`
defaults to `MA` and `--board` defaults to `live`. Use `--board interim` for the frozen
full-test interim board. Do not assume that a remote temporary path will remain
available later. Saved manifests contain file hashes and failed requests remain visible.
