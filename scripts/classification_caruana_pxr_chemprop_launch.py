"""Extend the audited Caruana pool with strict PXR Chemprop readouts."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parent = root / "artifacts/runs/classification/20260912T141520708830Z-97df4c71/inputs"
nested = root / "artifacts/runs/classification/20260913T161139263699Z-41ccf0ec/outputs"
inputs = sorted(parent.glob("*.parquet")) + [
    nested / "outer-oof.parquet",
    nested / "inner-probabilities.parquet",
    nested / "fold-test-probabilities.parquet",
]
assert len(inputs) == 24
for path in inputs:
    assert path.is_file(), path
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Rerun nested sparse and Caruana selection after adding five strict PXR "
        "Chemprop readouts as peer candidates to the existing 24-model pool"
    ),
    "--config",
    str(root / "configs/classification-round/caruana-pxr-chemprop-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_caruana.py"],
    cwd=root,
    check=True,
)
