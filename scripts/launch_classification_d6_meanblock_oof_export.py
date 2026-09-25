"""Reproduce and export the current D6 cross-fitted OOF prediction."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
historical = root / "artifacts/runs/classification/20260912T121146867140Z-b0a05986"
fixed = root / "artifacts/runs/classification/20260913T123045515718Z-0df22d54/outputs"
raw = root / "artifacts/runs/classification/20260915T112944344816Z-69abecf5/outputs"
selected = root / "artifacts/runs/classification/20260915T111425299454Z-bea02fa4/outputs"
sources = {
    "reference-oof.parquet": historical / "outputs/ensemble-oof.parquet",
    "lgbm-oof.parquet": historical / "inputs/000-oof.parquet",
    "pfn-oof.parquet": historical / "inputs/003-oof-probabilities.parquet",
    "graph-oof.parquet": historical / "inputs/006-classification-oof.parquet",
    "fixed-oof.parquet": fixed / "oof.parquet",
    "raw600-oof.parquet": raw / "oof-probabilities.parquet",
    "selections.parquet": selected / "selections.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Reproduce the submitted D6 mean-block OOF and export keyed predictions",
    "--config",
    str(root / "configs/classification-round/d6-raw600-controlled-replacement-export-v1.json"),
]
for alias, path in sources.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/classification_d6_raw600_controlled_replacement.py"])
subprocess.run(command, cwd=root, check=True)
