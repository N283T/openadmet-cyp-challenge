# OpenADMET CYP Challenge

Reproducible local experiment infrastructure using pixi, PostgreSQL + RDKit, and
Polars. The two active tracks share molecular data and feature code:

- `regression`: direct-inhibition pIC50 for CYP1A2, CYP2C9, CYP2D6, CYP3A4.
- `classification`: boolean TDI predictions for CYP2D6 and CYP3A4.

Official [Space](https://huggingface.co/spaces/openadmet/cyp-challenge) and
[tutorial](https://github.com/OpenADMET/CYP-Challenge-Tutorial).
Prospective structure-aware modeling has a [four-target structure-data inventory](docs/structure-data.md):
TogoMCP-discovered experimental CYP/HEM structures, canonical AFDB references,
pinned LiteFold full-length MSAs, HEM/ligand-preserving AFDB-aligned multi-block CIFs,
and ligand-contact pocket definitions. Frozen `boltz-trunk` extraction is complete for 10,779 CYP/HEM/ligand complexes;
the first [pooling/readout experiment](docs/analysis/2026-10-01-boltz-trunk-pooling-protocol.md)
is complete. The [incumbent-matched redo](docs/analysis/2026-10-01-boltz-incumbent-redo.md)
is complete on the submitted regression cohorts and folds, including whole-member
ensemble placement and bank reselection. The user-adopted [1A2/D6 trunk regression
composition](docs/analysis/2026-10-01-boltz-regression-adoption.md) was submitted on
October 1 and improves Live ST-RAE for both changed targets against the previous
best scored submission. Boltz and Jupyter are not declared in this environment.
Pandas may exist transitively; project table processing uses Polars.

The [October 2 work summary](docs/analysis/2026-10-02-work-summary.md) consolidates
retained changes through October 2. On October 3 a D6 payload built from a
[test-like training pool](docs/analysis/2026-10-03-d6-testlike-pool.md), where
primary-screen non-hits become assigned negatives, scored .364248 MCC, the best
D6 score we have observed, above the previous .308635. A
[sizing analysis](docs/analysis/2026-10-03-d6-live-gain-analysis.md) attributes
most of that move to one scored molecule and places it inside an exchangeable
null, so the recipe is retained as the control without being established as
better. A4 is unchanged at .529911, rank 1, and
macro MCC is .447079. The [A4 diagnosis](docs/analysis/2026-10-03-a4-diagnosis.md)
explains why that endpoint is deliberately left alone.

For new regression work, use the explicit [current auxiliary configuration](configs/regression/current-log2fc.json)
and [current parent-CV readout](configs/regression/current-parent-readout.json).
These select the user-adopted Rustuna log2fc recipe; older experiment reports and
run configurations retain their original inputs. See [the adoption record](docs/analysis/2026-09-11-rustuna-adoption.md).
The [current TabICL and TabPFN models](configs/regression/current-foundation-models.json)
include fitted model assets, OOF predictions and five-fold averaged test predictions
after the [completed downstream replacement](docs/analysis/2026-09-11-foundation-replacement.md).

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

## Record substantive experiments

Keep substantive probes as scripts; small read-only inspections and preflight checks
can run directly. The wrapper saves intent before execution,
snapshots source/config/inputs/environment, then executes from the source snapshot.
All model/feature commands must read their declared input snapshots and write to
`CYP_OUTPUT_DIR`. Scratch uses the persistent run-local `TMPDIR`.

```sh
pixi run cyp run --track shared \
  --purpose "Inspect blinded dataset dimensions" \
  --input test.csv=data/raw/cyp-challenge-TEST-BLINDED.csv \
  -- python scripts/smoke.py

pixi run cyp replay RUN_ID
```

Record the outcome of a scientifically informative discarded probe:

```sh
pixi run cyp note RUN_ID --decision discard --notes "No improvement; retain as negative evidence"
```

Use repeatable `--input ALIAS=FILE`, `--config FILE`, and optional `--issue URL`.
The resolved config snapshot is exposed as `CYP_CONFIG_PATH`. Input snapshots live
in `CYP_INPUT_DIR` under the declared aliases. Legacy bare `--input FILE` arguments
are numbered in argument order. Input and config paths are checked before a run is
allocated. The command is executed as an argv array without a shell. Logs are in the
printed run directory's `console.log`.

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

Frozen [CheMeleon fingerprints](docs/chemeleon.md) are also available: 2,048 dimensions
for all standardized compounds, with keyed Parquet and PostgreSQL storage.

Five-fold [log2fc pretraining from CheMeleon](docs/log2fc-pretraining.md) is complete,
with selected checkpoints and fold-specific embeddings ready for readout experiments.
The first [embedding-only TabPFN readout](docs/tabpfn-readout.md) is complete on
the fixed UMAP folds, with saved fitted states and audited OOF predictions.
The [PCA + predicted-log2fc TabPFN comparison](docs/stacked-tabpfn.md) adds nested
auxiliary cross-fitting and a compact 260-dimensional readout.
The [Mordred and direct-FT screen](docs/model-screen.md) compares descriptor
augmentation with directly fine-tuning CheMeleon on the challenge targets.
The [TabPFN 3.0 comparison](docs/tabpfn-v3.md) evaluates both tracks and retains
2.6 for the current feature recipe.

[External CYP auxiliaries](docs/external-data.md) preserve ChEMBL 36 IC50 and
PubChem AID1851 continuous responses as separate targets, with source lineage and
all challenge connectivity blocks excluded from external training.

The [regression calibration screen](docs/analysis/2026-09-11-calibration.md)
tests PXR-inspired corrections on the existing TabPFN predictions, with interval
losses adapted to CYP scoring and explicit evaluation limitations.
The subsequent [distribution and ensemble comparison](docs/analysis/2026-09-11-distribution-ensembles.md)
tests predictive quantiles and existing model diversity. The
[interval-learning experiment](docs/analysis/2026-09-11-interval-learning.md)
fits confidence-bound targets and interval-aware tree models.

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

For classification replacements, use the
[candidate evaluation protocol](docs/classification-evaluation.md) to select the
score recipe and positive-call rate on held-out folds. The earlier
[candidate preflight](docs/classification-preflight.md) remains available as a
diagnostic for fixed-rate swaps; its thresholds do not establish the best rate
or whether a candidate is worth a prospective test.
The [AP-mass90 full-pool replay](docs/analysis/2026-09-24-d6-ap-full-pool-replay.md)
shows how labelled-only call selection and full-fit ensemble weights can
overstate D6 OOF gains. The evaluator requires complete-pool predictions for
deployment decisions; labelled-only results need an explicit diagnostic opt-in.

Use `pixi run cyp submit SUBMISSION_ID --dry-run` to validate the prepared file and
private profile and review per-target changes and cooldown timing. Obtain explicit
user approval for that specific payload before running the live
`pixi run cyp submit SUBMISSION_ID` command. It requires the full prepared ID: typed in an
interactive terminal, or passed as `--confirm "SUBMIT SUBMISSION_ID"` when the user has approved
that exact payload in conversation and no interactive terminal is available. A successful upload consumes a 12-hour opportunity for that track.
Prefer combining worthwhile changes across all targets when practical; an unchanged
target is advisory, not a blocker. The sender
reads `.private/submission-profile.json` (Git-ignored), saves the MA and every
endpoint leaderboard before upload, and saves the API receipt. If any leaderboard
snapshot fails, the command stops before calling the submission API.
See [baseline submission preparation](docs/submission-baseline.md).

For manual uploads or subsequent scoring evidence, attach a saved response using
`pixi run cyp submission-event ID --status sent --evidence FILE --notes TEXT`.
Later append `accepted` or `rejected`; use `unknown` for uncertain delivery.
Never automatically retry a potentially successful upload.

## Public API

```sh
pixi run api-inspect
pixi run cyp leaderboard --track regression
pixi run cyp leaderboard --track classification
pixi run cyp leaderboard --track regression --board interim
pixi run cyp leaderboard --track regression --board final
```

Inspection persists live schema/config plus revision-pinned public source.
Leaderboard download resolves its endpoint from the live layout and saves actual
files to that invocation's persistent directory, not `/tmp`. Use `--endpoint`
with an official endpoint name to retrieve its leaderboard. The default is `MA`.
Since the challenge closed, the Space shows Final instead of Live boards; use `--board final`.
See [API notes](docs/api.md) for the current remote submission contract.

## Official references

Read the [official source archive](references/official/README.md) for saved challenge
instructions, CYP blog posts and assay documentation, and the PXR retrospective.
Original text and figures are separated from project docs and experiment outputs;
updates retain prior timestamped versions.

Participant write-ups are saved separately in the
[community source archive](references/community/README.md), including the
SuperCowPowers CYP working notes. These preserve author claims rather than
independently validated project findings.

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
