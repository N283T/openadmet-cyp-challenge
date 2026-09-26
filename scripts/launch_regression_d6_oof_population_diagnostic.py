"""Record a descriptive D6 OOF population diagnostic for the current ensemble."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
regression_runs = root / "artifacts/runs/regression"
inputs = {
    "control-oof.parquet": (
        root / "artifacts/tracked/d6-depth-atomic-20260925/ensemble-oof.parquet"
    ),
    "reference-oof.parquet": (
        regression_runs / "20260912T072303592078Z-07f2749a/outputs/external-oof.parquet"
    ),
    "analog-pairs.parquet": (
        root / "artifacts/tracked/d6-analog-expansion-validation-20260924"
        / "anchor-analog-pairs.parquet"
    ),
}
command = [
    "pixi", "run", "cyp", "run", "--track", "regression",
    "--purpose", "Diagnose current CYP2D6 parent-OOF errors by activity and other-CYP analog proxy",
    "--config", str(root / "configs/regression/d6-oof-population-diagnostic-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_oof_population_diagnostic.py"],
    cwd=root,
    check=True,
)
