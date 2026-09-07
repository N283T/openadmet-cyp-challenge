# CYP Challenge

- Keep the GitHub repository private. Change visibility only on explicit user instruction;
  public release will be considered when preparing the model report.

- Use `pixi run`; use Polars for table processing. No Jupyter or Boltz for now.
- Run every experiment, including throwaway probes, through `pixi run cyp run`.
  Supply a purpose, configuration, and all input files. Keep scripts under `scripts/`
  or `src/`; never leave the only copy of work in `/tmp` or a terminal command.
- Write outputs to `CYP_OUTPUT_DIR`; scratch files to `TMPDIR` (persistent per run).
  Failed, interrupted, and discarded experiments stay recorded. Never overwrite
  completed runs, candidates, submissions, source data, or applied migrations.
- Shared molecular data/features; explicit `regression`/`classification` tracks.
  Predictions must carry molecule IDs, endpoints, and split/fold IDs, not row indices.
- DB = structured records. Issues = hypotheses, conclusions, and decisions with IDs.
  Do not put evolving experiment history or best-model claims in this file.
- Keep model weights/data/runtime files out of Git. Never broad-clean ignored files.
- Use feature branches. Run `pixi run check` and `pixi run test` for infrastructure edits.
- Read [README.md](README.md) and [docs/infrastructure.md](docs/infrastructure.md).
