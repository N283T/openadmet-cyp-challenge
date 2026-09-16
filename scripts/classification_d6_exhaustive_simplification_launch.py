"""Launch exhaustive performance-first D6 simplification."""

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
    "Exhaustively simplify all D6 core members with held-fold-excluded selection",
    "--config",
    str(root / "configs/classification-round/d6-exhaustive-simplification-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_simplified_hierarchy.py"],
    cwd=root,
    check=True,
)
