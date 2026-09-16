"""Launch A4 raw replacement and exhaustive fixed-subset simplification."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
inputs = [
    runs / "20260913T144948028103Z-cbb9ed62/outputs/member-oof.parquet",
    runs / "20260913T144948028103Z-cbb9ed62/outputs/member-test.parquet",
    runs / "20260913T152746349812Z-c510a24e/outputs/oof.parquet",
    runs / "20260913T152746349812Z-c510a24e/outputs/test.parquet",
    runs / "20260915T123533598879Z-a459b098/outputs/outer-oof.parquet",
    runs / "20260915T123533598879Z-a459b098/outputs/test-probabilities.parquet",
    runs / "20260915T073340818290Z-d781ddd2/outputs/oof.parquet",
    runs / "20260915T073340818290Z-d781ddd2/outputs/test.parquet",
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
    "Replace A4 PCA readouts with raw controls and simplify the resulting pool",
    "--config",
    str(root / "configs/classification-round/a4-raw-simplification-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_a4_raw_simplification.py"],
    cwd=root,
    check=True,
)
