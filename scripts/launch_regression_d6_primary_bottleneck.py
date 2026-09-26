"""Record measured-versus-predicted D6 primary-response diagnostic."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
    "parent-assignments.parquet": (
        runs / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    "measurements.parquet": (
        runs / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet"
    ),
    "regression-labels.parquet": (
        runs / "regression/20260926T004020214692Z-b2f253b0/inputs/regression-labels.parquet"
    ),
    "nested-primary.parquet": (
        runs / "shared/20260912T082554377951Z-d272774c/outputs/nested-log2fc.parquet"
    ),
    "control-oof.parquet": (
        runs / "regression/20260926T013148480112Z-72cf4682/outputs/aligned-oof.parquet"
    ),
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Separate D6 primary-response information from primary prediction error",
    "--config",
    str(root / "configs/regression/d6-primary-bottleneck-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_primary_bottleneck.py"],
    cwd=root,
    check=True,
)
