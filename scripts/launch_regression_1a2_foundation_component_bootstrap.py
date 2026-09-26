"""Record grouped uncertainty for the selected CYP1A2 atomic combination."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
inputs = {
    "component-oof.parquet": (
        runs / "20260926T015505095663Z-0acd3ef2/outputs/oof.parquet"
    ),
    "parent-assignments.parquet": (
        runs / "20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    "scaffold-assignments.parquet": (
        runs / "20260918T021731984021Z-11c1a26e/outputs/parent-assignments.parquet"
    ),
}
command = [
    "pixi", "run", "cyp", "run", "--track", "regression",
    "--purpose", "Quantify grouped uncertainty of CYP1A2 atomic combination on both splits",
    "--config", str(root / "configs/regression/1a2-foundation-component-bootstrap-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_1a2_foundation_component_bootstrap.py"],
    cwd=root,
    check=True,
)
