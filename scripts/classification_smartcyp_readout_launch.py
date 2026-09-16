"""Launch controlled dual-split SMARTCyp readouts."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
base = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
inputs = [
    base / "compounds.parquet",
    base / "measurements.parquet",
    base / "test_mapping.parquet",
    root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    root
    / (
        "artifacts/runs/classification/20260914T044537363003Z-5e7f9e02/outputs/"
        "smartcyp-features.parquet"
    ),
]
config = root / "configs/classification-round/smartcyp-readout-v1.json"
for path in [config, *inputs]:
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
    "Compare SMARTCyp SoM summaries and top-SoM alert interactions on UMAP and scaffold folds",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_smartcyp_readout.py"],
    cwd=root,
    check=True,
)
