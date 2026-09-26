"""Record fixed D6 upper-bound readout blends against matched OOF controls."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
    "upper-oof.parquet": (runs / "regression/20260926T041418410211Z-076c47f6/outputs/oof.parquet"),
    "parent-control.parquet": (
        runs / "regression/20260926T013148480112Z-72cf4682/inputs/control-oof.parquet"
    ),
    "scaffold-control.parquet": (
        runs / "regression/20260918T021832961989Z-c89b25fa/outputs/jeremy_anchor_all-oof.parquet"
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
    "Check fixed placement and low-weight addition of one-sided D6 model",
    "--config",
    str(root / "configs/regression/d6-primary-bound-blend-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_primary_bound_blend.py"], cwd=root, check=True
)
