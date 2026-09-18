"""Register the scaffold secondary split with the regression anchor policy."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
config = root / "configs/regression/parent-scaffold-secondary-v1.json"
inputs = {
    "base-assignments.parquet": runs
    / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "primary-parent-assignments.parquet": runs
    / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Materialize parent-retained scaffold folds as a secondary regression ordering check",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/prepare_parent_secondary_split.py"],
    cwd=root,
    check=True,
)
