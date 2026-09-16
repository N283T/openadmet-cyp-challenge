"""Launch the performance-first D6 hierarchy simplification audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
classification = root / "artifacts/runs/classification"
historical = classification / "20260912T121146867140Z-b0a05986"
fixed = classification / "20260913T123045515718Z-0df22d54/outputs"
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
    "Simplify D6 for predictive performance rather than legacy-vector reproduction",
    "--config",
    str(root / "configs/classification-round/d6-simplified-hierarchy-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_simplified_hierarchy.py"],
    cwd=root,
    check=True,
)
