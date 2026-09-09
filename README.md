# OpenADMET CYP Challenge

Reproducible local experiment infrastructure using pixi, PostgreSQL + RDKit, and
Polars. The two active tracks share molecular data and feature code:

- `regression`: direct-inhibition pIC50 for CYP1A2, CYP2C9, CYP2D6, CYP3A4.
- `classification`: boolean TDI predictions for CYP2D6 and CYP3A4.

Official [Space](https://huggingface.co/spaces/openadmet/cyp-challenge) and
[tutorial](https://github.com/OpenADMET/CYP-Challenge-Tutorial).
Structure support and Boltz are deferred. No Jupyter dependency is declared.
Pandas may exist transitively; project table processing uses Polars.

## Setup

```sh
mkdir -p .cache/tmp
pixi install --locked
pixi run db-init
pixi run db-start
pixi run db-migrate
pixi run cyp register-raw
pixi run check
pixi run test
```

PostgreSQL uses a project-local socket in `.runtime/socket`, port 5434, and database
`cyp_challenge`. TCP listening is disabled. It does not reuse the PXR database.
`pixi run db-stop` stops it. `pixi run db-backup` writes a durable custom-format dump.

## Record every experiment

Keep even quick probes as scripts. The wrapper saves intent before execution,
snapshots source/config/inputs/environment, then executes from the source snapshot.
All model/feature commands must read their declared input snapshots and write to
`CYP_OUTPUT_DIR`. Scratch uses the persistent run-local `TMPDIR`.

```sh
pixi run cyp run --track shared \
  --purpose "Inspect blinded dataset dimensions" \
  --input data/raw/cyp-challenge-TEST-BLINDED.csv \
  -- python scripts/smoke.py

pixi run cyp replay RUN_ID
```

Record the outcome of a discarded probe without deleting it:

```sh
pixi run cyp note RUN_ID --decision discard --notes "No improvement; retain as negative evidence"
```

Use repeatable `--input FILE`, `--config FILE`, and optional `--issue URL`.
The resolved config snapshot is exposed as `CYP_CONFIG_PATH`. Input snapshots live
in `CYP_INPUT_DIR` and are numbered in argument order. The command is executed as
an argv array without a shell. Logs are in the printed run directory's `console.log`.

Replay checks recorded hashes and the current lock, reuses the saved source and
inputs, and creates a new run. Compare output hashes or task-specific tolerances;
replay is not a guarantee of deterministic GPU numerics. Model code must record
all seeds and actual effective parameters. See [infrastructure details](docs/infrastructure.md).

## Molecular preprocessing

Shared molecular data are standardized with the ChEMBL structure pipeline and
exact duplicates are resolved with source provenance. Python RDKit's full 217
registered descriptors and Mordred's 1,613 2D descriptors are stored in PostgreSQL
and versioned Parquet outputs. See [preprocessing](docs/preprocessing.md) for the
policy, dataset/feature run IDs, exclusions, and reproducible commands.
No splits or model fitting are included in this preprocessing stage.
See [split baselines](docs/split-baselines.md) for the subsequent random/scaffold/UMAP
LightGBM comparison and the planned log2fc → frozen representation → TabPFN workflow.

## Candidates and submissions

```sh
pixi run cyp candidate --track regression \
  --file artifacts/runs/regression/RUN_ID/outputs/predictions.parquet \
  --reference data/raw/cyp-challenge-TEST-BLINDED.csv \
  --run RUN_ID --notes "Candidate rationale; transformations recorded in source run"

pixi run cyp prepare-submission CANDIDATE_ID
```

Candidate creation validates 750 exact test IDs/SMILES, track columns and values,
and references completed runs. Perform ensembling/calibration/thresholding as a
recorded run first. Classification probabilities belong to run outputs; upload
files contain booleans. A preparation creates a separate immutable file copy;
**it does not call the remote submission API**.

After an actual upload, attach a saved API response/receipt using
`pixi run cyp submission-event ID --status sent --evidence FILE --notes TEXT`.
Later append `accepted` or `rejected`; use `unknown` for uncertain delivery.
Never automatically retry a potentially successful upload.

## Public API

```sh
pixi run api-inspect
pixi run cyp leaderboard --track regression
pixi run cyp leaderboard --track classification
```

Inspection persists live schema/config plus revision-pinned public source.
Leaderboard download resolves its endpoint from the live layout and saves actual
files to that invocation's persistent directory, not `/tmp`. Use `--endpoint`
with an official endpoint name to retrieve its leaderboard. The default is `MA`.
See [API notes](docs/api.md) for the current remote submission contract.

## Official references

Read the [official source archive](references/official/README.md) for saved challenge
instructions, CYP blog posts and assay documentation, and the PXR retrospective.
Original text and figures are separated from project docs and experiment outputs;
updates retain prior timestamped versions.

## Layout

```text
src/cyp_challenge/       shared infrastructure and future model implementations
configs/{regression,classification}/  track experiment configurations
scripts/                checked-in probes and operations
db/migrations/          ordered, hash-checked SQL migrations
tests/                 focused infrastructure tests
docs/                  setup, architecture, API reference
references/official/   original official source snapshots, grouped by CYP/PXR
data/{raw,processed,features}/       shared versioned data and features
artifacts/runs/<track>/<id>/         all attempts, including failed probes
artifacts/candidates/<track>/<id>/   frozen upload candidates
artifacts/submissions/<track>/<id>/  separate upload copies and receipt events
artifacts/api/                      API contracts, responses and downloaded files
```

DB stores structured records; the local [research log](docs/research/README.md)
stores hypotheses, results and decisions using run/candidate IDs. Reports and log
entries default to English; Issues are optional for bounded implementation tasks. The private repository is https://github.com/N283T/openadmet-cyp-challenge.
Keep it private until an explicit publication decision for the model report.
Routine evaluation uses UMAP only. Revisit the split if observed LB results diverge
substantially from CV. See the research log for current decisions and next steps.
