"""Launch the recorded compact dual-split classification screen."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    root / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    root
    / "artifacts/runs/shared/20260912T051543210441Z-a3ea2ec3/outputs/external-embeddings.parquet",
    root
    / "artifacts/runs/shared/20260912T052853359792Z-2c0a0953/outputs/external-embeddings.parquet",
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
        "Compare prespecified Morgan, frozen-Chemprop and PLS-compressed readouts "
        "with honest local-residual correction on UMAP and scaffold folds"
    ),
    "--config",
    str(root / "configs/classification-round/dual-split-compact-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_dual_split_compact.py"],
    cwd=root,
    check=True,
)
