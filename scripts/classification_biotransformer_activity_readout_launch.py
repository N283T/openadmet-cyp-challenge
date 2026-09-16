"""Launch classification readout for fold-safe BioTransformer product scores."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
inputs = [
    shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    shared / "20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    root
    / (
        "artifacts/runs/classification/20260915T065421221917Z-7d41918f/outputs/"
        "product-activity-features.parquet"
    ),
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
    "Evaluate fold-safe direct-pIC50 and log2fc summaries over BioTransformer products",
    "--config",
    str(
        root
        / "configs/classification-round/biotransformer-product-activity-readout-v1.json"
    ),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_biotransformer_activity_readout.py"],
    cwd=root,
    check=True,
)
