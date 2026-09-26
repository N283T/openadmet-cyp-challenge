"""Record group-bootstrap and fold audit for the proposed D6 bound blend."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
    "blend-oof.parquet": (runs / "regression/20260926T041934073362Z-5e9e8cb8/outputs/oof.parquet"),
    "parent-assignments.parquet": (
        runs / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    "scaffold-assignments.parquet": (
        runs / "regression/20260918T021731984021Z-11c1a26e/outputs/parent-assignments.parquet"
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
    "Audit group stability and uncertainty of selected D6 primary-bound blend",
    "--config",
    str(root / "configs/regression/d6-primary-bound-audit-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_primary_bound_audit.py"], cwd=root, check=True
)
