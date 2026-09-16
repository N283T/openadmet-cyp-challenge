# CYP Challenge

- Keep the GitHub repository private. Change visibility only on explicit user instruction;
  public release will be considered when preparing the model report.

- Use `pixi run`; use Polars for table processing. No Jupyter or Boltz for now.
- Run substantive experiments through `pixi run cyp run`. Small read-only inspections
  and preflight checks do not need a recorded run. Supply a purpose, configuration,
  and all input files. Use explicit input aliases for new or edited launchers instead
  of predicting recorder-generated numeric prefixes. Keep scripts under `scripts/`
  or `src/`; never leave the only copy of substantive work in `/tmp` or a terminal
  command.
- Write outputs to `CYP_OUTPUT_DIR`; scratch files to `TMPDIR` (persistent per run).
  Preserve completed runs, candidates, submissions, scientifically informative
  negative results, source data, and applied migrations. Failed or interrupted runs
  with no usable scientific result may be removed after recording a concise reason.
  Never overwrite an immutable record that is being retained.
- Shared molecular data/features; explicit `regression`/`classification` tracks.
  Predictions must carry molecule IDs, endpoints, and split/fold IDs, not row indices.
- DB = structured records. [Research log](docs/research/README.md) = hypotheses,
  results (including failures), and decisions with run IDs; write logs/reports in English.
  Use UMAP only for routine evaluation; revisit if LB substantially diverges from CV.
  Do not put evolving experiment history or best-model claims in this file.
- Track project-authored source, configs, reports, compact derived data, model weights,
  and runtime results whenever they are safe, legally shareable, and fit ordinary Git
  hosting limits. Do not commit secrets, personal information, proprietary inputs, or
  third-party copyrighted model/data payloads. Do not use Releases or Git LFS merely
  to force oversized artifacts into the repository. Promote finalized, manageable run
  outputs to `artifacts/tracked/`; keep oversized or external payloads local with
  provenance and checksums. Inspect ignored files before any targeted cleanup.
- Use feature branches. Run `pixi run check` and `pixi run test` for infrastructure edits.
- Official source snapshots: [references/official/README.md](references/official/README.md).
  Keep source text separate from interpretations; preserve previous versions on refresh.
- Read [README.md](README.md) and [docs/infrastructure.md](docs/infrastructure.md).
