"""Record labelled-OOF depth additions to currently submitted TDI recipes."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
sources = {
    "depth-oof.parquet": root
    / "artifacts/runs/shared/20260925T065626160287Z-8017d516/outputs/oof.parquet",
    "base-d6.parquet": (
        root
        / "artifacts/runs/classification/20260925T130752164378Z-d40a0d9c/outputs"
        / "crossfit-mean-block-oof.parquet"
    ),
    "base-a4.parquet": root
    / "artifacts/runs/classification/20260915T125001009276Z-d5adff22/outputs/oof.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Screen fixed CheMeleon depth additions to submitted TDI OOF rankings",
    "--config",
    str(root / "configs/screen/chemeleon-depth-incremental-classification-v1.json"),
]
for alias, path in sources.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/evaluate_chemeleon_depth_incremental_classification.py"])
subprocess.run(command, cwd=root, check=True)
