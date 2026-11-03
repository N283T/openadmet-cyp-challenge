"""Launch the CheMeleon depth super-ensemble analysis for CYP3A4."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
regression = root / "artifacts/runs/regression"
inputs = {
    "control-3a4.parquet": regression / "20260917T130718808881Z-dfa02730/outputs/oof.parquet",
    "analog-pairs.parquet": root
    / "artifacts/tracked/d6-analog-expansion-validation-20260924/anchor-analog-pairs.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Super-ensemble of 16 CheMeleon checkpoints x 6 depths x 2 readouts for CYP3A4",
    "--config",
    str(root / "configs/regression/depth-super-ensemble-3a4-v1.json"),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/regression_depth_super_ensemble.py"], cwd=root, check=True
)
