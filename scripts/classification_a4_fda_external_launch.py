"""Launch the recorded A4 FDA external-task prediction-block experiment."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root
    / (
        "artifacts/runs/classification/20260913T143057530832Z-47a95438/outputs/"
        "fda-cyp3a4-tdi-records.parquet"
    ),
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
        "Evaluate non-overlapping FDA CYP3A4 TDI labels as an independent Morgan "
        "LR, ExtraTrees, or Tanimoto-SVM prediction block blended with the fixed "
        "A4 radius-Chemprop block on UMAP and scaffold folds"
    ),
    "--config",
    str(root / "configs/classification-round/a4-fda-external-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_a4_fda_external.py"], cwd=root, check=True
)
