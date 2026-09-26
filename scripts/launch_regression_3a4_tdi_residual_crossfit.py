"""Record leakage-isolated 3A4 residual transfer from classification OOF scores."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = {
    "regression-oof.parquet": root
    / "artifacts/runs/regression/20260917T130718808881Z-dfa02730/outputs/oof.parquet",
    "classification-oof.parquet": root
    / "artifacts/runs/classification/20260915T125001009276Z-d5adff22/outputs/oof.parquet",
    "parent-assignments.parquet": root
    / "artifacts/runs/regression/20260926T085732408866Z-0167cd1b/inputs/parent-assignments.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Test incremental TDI-classification signal in the scored 3A4 regression bank residuals",
    "--config",
    str(root / "configs/regression/3a4-tdi-residual-crossfit-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_3a4_tdi_residual_crossfit.py"],
    cwd=root,
    check=True,
)
