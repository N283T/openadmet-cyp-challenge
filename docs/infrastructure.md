# Infrastructure decisions

## Durable experimentation

Substantive experiments are recorded runs; small read-only inspections and preflight
checks are not. Validate declared input and config paths before allocating a run.
Failed commands and keyboard interruption produce `result.json`. Failed or interrupted
runs with no usable scientific result may be removed after recording a concise reason;
scientifically informative negative results remain durable. `cyp reconcile` repairs
final DB states from saved results after a DB outage; records without a result are
printed for human inspection because they may still be running. SIGKILL/power failure
cannot execute finalizers; intent, process identity and console logs remain available.

`pixi run cyp prune-runs --status failed` previews failed run payload cleanup. Adding
`--execute` replaces each run directory with a compact tombstone under
`artifacts/pruned-runs/`; the database keeps the failure status, command, failure detail,
notes, reclaimed size and tombstone path. Runs without a terminal result are not pruned.

Each run contains:

- Intent, purpose, track, optional Issue, command, timestamp and process identity.
- Saved source tree, Git status/commit/diff (commit can be absent during initial setup).
- pixi manifest/lock, installed-package inventory, Python and platform information.
- Copies and SHA-256 hashes of declared input files and configuration.
- Persistent scratch, outputs, console log and terminal result with file inventory.

Executing from saved source prevents ongoing source edits from changing a long
run. Only source/config/scripts/SQL and environment definitions are copied, not
arbitrary files from the working directory. Declare extra input files explicitly.
Use `CYP_INPUT_DIR`, `CYP_CONFIG_PATH`, and `CYP_OUTPUT_DIR`; commands that bypass
these and read mutable absolute paths do not inherit reproducibility guarantees.
Do not put credentials in command arguments, configs, or checked-in source.

New or edited launchers should pass inputs as `--input alias.ext=/path/to/file` and
refer to `alias.ext` from the run config. Bare `--input /path/to/file` remains
supported for older launchers and receives the legacy numeric prefix. Input paths,
aliases and the config path are validated before a run ID is allocated, so basic
name/path mistakes do not create failed run records.

Large immutable feature assets may eventually use content-addressed storage to
avoid per-run copies; this initial version intentionally copies declared inputs.
All candidate transformations must themselves be recorded runs, including blend
weights, calibration, threshold selection and output serialization.

## Track and endpoint boundaries

Share compound identities, raw data, normalization and feature implementations.
Track is a required run/candidate/submission field. A `shared` run can supply either
track; regression and classification candidates cannot cite the opposite track.
Endpoint is independent from track: one model can serve multiple endpoints.
OOF/held-out predictions should use explicit `Molecule_Name`, `endpoint`, `split_id`,
`fold`, `prediction` columns. Never join by implicit row index. Store split assignments
as immutable input files and preserve missing-label masks per endpoint. TDI training
files also contain regression measurements: establish overlap/precedence rules
before building merged labels; do not blindly concatenate the available CSVs.

Metrics have run/endpoint/split/fold/name keys. Data loading/model training and the
competition's exact soft-threshold evaluator are subsequent modeling work; the
infrastructure does not substitute ordinary RAE for ST-RAE.

## State and provenance

- PostgreSQL: run status/config/results, metric rows, input/API inventories, candidate
  lineage and submission records/events.
- Files: model weights, keyed predictions, source/input snapshots, raw API outputs.
- [Local research log](research/README.md): hypotheses, interpretations, rejected
  approaches and decisions, with stable entry IDs and run references. English by default.
- Issues: optional bounded implementation tasks or actionable bugs.
- AGENTS.md: short durable operating rules and links only.

Candidates get a fresh ID and a read-only file copy. Submissions get another fresh
ID and file copy with matching SHA-256. Preparation is separate from remote sending.
Events append receipts/status without rewriting the original manifest or file.
Before a live send, the submission command downloads MA and every endpoint
leaderboard for the selected track. All downloads must finish successfully before
the upload intent is created or the submission API is called; their artifact paths
are stored in the attempt record.
Live sending also requires explicit user approval for the specific prepared payload
and an interactive terminal confirmation. The dry run reports changed molecule counts
per target against the most recent sent payload, when that file is available. An
unchanged target is advisory: combine worthwhile changes across targets when practical.
The sender blocks a new attempt within 12 hours of the latest locally recorded sent
attempt on the same track. Manual or unrecorded remote sends must be checked separately.
SHA-256 verification catches accidental changes; chmod is not a security boundary
against the file owner. Run output files are inventoried rather than OS-locked.

## Data and recovery

The five existing CSVs were moved intact into `data/raw/`. Their original download
revision is unknown, so the import manifest says so; it does not invent provenance.
New official data should be downloaded from revision-pinned URLs and registered
as a new snapshot. Do not overwrite raw files or an applied migration.

`pixi run db-backup` saves a PostgreSQL custom-format dump with hashes. Restore to
an empty database using `pg_restore --no-owner` after creating the target database
and installing the same RDKit extension. Back up `artifacts/`, `data/`, and Git
separately: a DB dump contains metadata, not the large files themselves. Stored
paths are absolute, so relocation requires path remapping. External/off-device
backup destination is not configured yet. Durable local storage prevents `/tmp`
cleanup losses, but does not provide disk-failure protection.

## Environment boundaries

CPU/tabular foundation installed now: Polars, Arrow, NumPy, SciPy, scikit-learn,
RDKit, PostgreSQL cartridge, psycopg, LightGBM, XGBoost, CatBoost, Optuna,
Matplotlib, Plotly, YAML, HTTPX, Hugging Face Hub, Gradio client, Rich, pytest, Ruff.
PyTorch and Chemprop are now installed through pixi for frozen CheMeleon extraction;
the extraction run records the resolved versions and GPU runtime. Other foundation
model libraries can be added when the model choices are known. Boltz and Jupyter are not declared.
The lock is the resolved environment authority; do not use ad-hoc pip installs.

### Large retained immutable inputs

`cyp run --input-mode hardlink` explicitly creates same-filesystem hardlink
snapshots instead of copying large retained immutable assets. Each input still
has an alias and SHA-256 in the run manifest. Source permissions are unchanged.
The ordinary copy mode remains the default. Hardlinks share an inode, so callers
must preserve source immutability; they protect against path removal, not in-place
writes. Completion compares re-inventoried input hashes against initial hashes
and marks a changed input as a failed run. Replay retains hardlink storage when
all original input snapshots used it. Cross-filesystem linking fails rather than
silently allocating a large copy.
