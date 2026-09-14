"""Launch fixed-Morgan classical CYP2D6 model comparison."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = [
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    root / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    root
    / (
        "artifacts/runs/classification/20260912T121146867140Z-b0a05986/outputs/ensemble-oof.parquet"
    ),
    root
    / (
        "artifacts/runs/classification/20260912T121146867140Z-b0a05986/outputs/"
        "ensemble-test.parquet"
    ),
    root / "artifacts/runs/classification/20260913T112957931554Z-32d00f6e/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T112957931554Z-32d00f6e/outputs/test.parquet",
]
assert all(path.is_file() for path in paths)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Compare standard classifiers with fixed Morgan radius2 features for CYP2D6",
    "--config",
    str(root / "configs/classification-round/d6-fixed-fp-models-v1.json"),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_fixed_fp_models.py"],
    cwd=root,
    check=True,
)
