"""Launch nested safe stacking around the current fixed-ratio expert."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
classification = root / "artifacts/runs/classification"
inputs = [
    classification / "20260915T121009912691Z-73198c15/outputs/prediction-bank.parquet",
    classification / "20260914T063614433212Z-c0dc3f32/outputs/oof.parquet",
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
    "Evaluate nested anchor-preserving stackers after unrestricted stacking",
    "--config",
    str(root / "configs/classification-round/anchored-stacking-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_anchored_stacking.py"],
    cwd=root,
    check=True,
)
