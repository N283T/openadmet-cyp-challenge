"""Launch the recorded classification split-support audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
inputs = [
    shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    shared / "20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
]
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
        "Compare existing random, scaffold and UMAP folds plus the leave-one-out "
        "support ceiling against blind test before fixing secondary validation"
    ),
    "--config",
    str(root / "configs/classification-round/split-support-audit-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_split_support_audit.py"],
    cwd=root,
    check=True,
)
