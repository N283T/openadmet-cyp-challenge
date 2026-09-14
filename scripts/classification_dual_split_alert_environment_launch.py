"""Launch the recorded alert-centered Morgan-environment experiment."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    root / "artifacts/runs/classification/20260913T141558951845Z-1bdad261/outputs/oof.parquet",
    root / "artifacts/runs/classification/20260913T141558951845Z-1bdad261/outputs/test.parquet",
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
        "Test collision-free Morgan environments centered on explicit TDI-alert "
        "atoms or their one-bond neighborhoods as complementary blocks on fixed "
        "UMAP and scaffold folds"
    ),
    "--config",
    str(root / "configs/classification-round/dual-split-alert-environment-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_dual_split_alert_environment.py"],
    cwd=root,
    check=True,
)
