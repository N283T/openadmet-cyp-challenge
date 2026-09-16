"""Launch correlation and redundancy diagnostics for the selected D6 design."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
historical = root / "artifacts/runs/classification/20260912T121146867140Z-b0a05986"
fixed = root / "artifacts/runs/classification/20260913T123045515718Z-0df22d54/outputs"
inputs = [
    historical / "outputs/ensemble-oof.parquet",
    historical / "inputs/000-oof.parquet",
    historical / "inputs/003-oof-probabilities.parquet",
    historical / "inputs/006-classification-oof.parquet",
    fixed / "oof.parquet",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Measure prediction correlation and leave-one-out redundancy in the selected D6 design",
    "--config",
    str(root / "configs/classification-round/d6-selected-correlation-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_selected_correlation.py"],
    cwd=root,
    check=True,
)
