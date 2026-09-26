"""Record fixed D6 robust-bound blends against matched OOF controls."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
    "upper-oof.parquet": runs / "regression/20260926T085732408866Z-0167cd1b/outputs/oof.parquet",
    "parent-control.parquet": runs
    / "regression/20260926T013148480112Z-72cf4682/inputs/control-oof.parquet",
    "scaffold-control.parquet": runs
    / "regression/20260918T021832961989Z-c89b25fa/outputs/jeremy_anchor_all-oof.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Compare fixed robust-bound blends with parent and scaffold D6 controls",
    "--config",
    str(root / "configs/regression/d6-primary-robust-blend-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_primary_robust_blend.py"],
    cwd=root,
    check=True,
)
