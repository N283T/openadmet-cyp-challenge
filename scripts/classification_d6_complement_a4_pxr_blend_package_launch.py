"""Launch assembly of the selected two-endpoint classification candidate."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root
    / (
        "artifacts/runs/classification/20260913T123420289849Z-b75917f6/outputs/"
        "fixed-morgan-model-selection.parquet"
    ),
    root
    / (
        "artifacts/runs/classification/20260913T153156161201Z-1954e23e/outputs/"
        "test-blends.parquet"
    ),
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
    root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
]
assert all(path.is_file() for path in inputs)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Assemble validated D6 live-ComplementNB six-swap plus A4 80/20 "
        "Caruana-Morgan-PXR-Admetica candidate; no upload"
    ),
    "--config",
    str(
        root
        / "configs/classification-round/d6-complement-a4-pxr-blend-package-v1.json"
    ),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command
    + [
        "--",
        "python",
        "scripts/classification_d6_complement_a4_pxr_blend_package.py",
    ],
    cwd=root,
    check=True,
)
