"""Launch add/swap audit for the adapted-Jeremy plus CheMeleon D6 screen."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
historical = runs / "20260912T121146867140Z-b0a05986"
inputs = {
    "reference-oof.parquet": historical / "outputs/ensemble-oof.parquet",
    "lgbm-oof.parquet": historical / "inputs/000-oof.parquet",
    "pfn-oof.parquet": historical / "inputs/003-oof-probabilities.parquet",
    "graph-oof.parquet": historical / "inputs/006-classification-oof.parquet",
    "fixed-oof.parquet": runs / "20260913T123045515718Z-0df22d54/outputs/oof.parquet",
    "raw600-oof.parquet": runs
    / "20260915T112944344816Z-69abecf5/outputs/oof-probabilities.parquet",
    "selections.parquet": runs
    / "20260915T111425299454Z-bea02fa4/outputs/selections.parquet",
    "candidate-oof.parquet": runs
    / "20260919T155849130417Z-0394b600/outputs/oof.parquet",
}
for path in inputs.values():
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
    "Audit adapted-Jeremy plus CheMeleon D6 backbones as additions and replacements",
    "--config",
    str(root / "configs/classification/d6-feature-combination-swap-audit-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_feature_combination_swap_audit.py"],
    cwd=root,
    check=True,
)
