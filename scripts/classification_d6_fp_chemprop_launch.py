"""Launch classical fingerprint and frozen-Chemprop CYP2D6 screen."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = [
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    root / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    root
    / "artifacts/runs/shared/20260912T051543210441Z-a3ea2ec3/outputs/external-embeddings.parquet",
    root
    / "artifacts/runs/shared/20260912T052853359792Z-2c0a0953/outputs/external-embeddings.parquet",
    root
    / (
        "artifacts/runs/classification/20260912T121146867140Z-b0a05986/outputs/ensemble-oof.parquet"
    ),
    root
    / (
        "artifacts/runs/classification/20260912T121146867140Z-b0a05986/outputs/"
        "ensemble-test.parquet"
    ),
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
    "Screen classical fingerprint and fingerprint-plus-Chemprop CYP2D6 rankers",
    "--config",
    str(root / "configs/classification-round/d6-fp-chemprop-v1.json"),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_fp_chemprop.py"],
    cwd=root,
    check=True,
)
