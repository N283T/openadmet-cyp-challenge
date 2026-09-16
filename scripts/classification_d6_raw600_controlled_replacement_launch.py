"""Launch controlled raw600 replacement in the selected D6 architecture."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
historical = root / "artifacts/runs/classification/20260912T121146867140Z-b0a05986"
fixed = root / "artifacts/runs/classification/20260913T123045515718Z-0df22d54/outputs"
raw = root / "artifacts/runs/classification/20260915T112944344816Z-69abecf5/outputs"
selected = root / "artifacts/runs/classification/20260915T111425299454Z-bea02fa4/outputs"
inputs = [
    historical / "outputs/ensemble-oof.parquet",
    historical / "inputs/000-oof.parquet",
    historical / "inputs/003-oof-probabilities.parquet",
    historical / "inputs/006-classification-oof.parquet",
    fixed / "oof.parquet",
    raw / "oof-probabilities.parquet",
    selected / "selections.parquet",
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
    "Isolate raw600 external-mixed TabPFN by fixed replacement and block averaging",
    "--config",
    str(root / "configs/classification-round/d6-raw600-controlled-replacement-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_raw600_controlled_replacement.py"],
    cwd=root,
    check=True,
)
