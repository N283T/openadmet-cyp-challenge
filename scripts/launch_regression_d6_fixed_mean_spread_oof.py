"""Record an OOF screen matching the prepared D6-only SD-0.90 sensitivity."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = {
    "blend-oof.parquet": root
    / "artifacts/runs/regression/20260926T085959697226Z-0dbf57d4/outputs/oof.parquet",
    "submitted.csv": root
    / "artifacts/submissions/regression/20260926T090548323924Z-46188564/predictions.csv",
    "variant.csv": root
    / "artifacts/submissions/regression/20260926T120124524160Z-01cea54f/predictions.csv",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Test exact fixed-mean D6 final-spread ratio on parent and scaffold held folds",
    "--config",
    str(root / "configs/regression/d6-fixed-mean-spread-oof-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_fixed_mean_spread_oof.py"],
    cwd=root,
    check=True,
)
