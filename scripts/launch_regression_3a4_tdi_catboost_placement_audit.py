"""Record a grouped uncertainty audit for the 3A4 TDI CatBoost blend."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = {
    "scored-bank-oof.parquet": root
    / "artifacts/runs/regression/20260917T130718808881Z-dfa02730/outputs/oof.parquet",
    "placed-oof.parquet": root
    / "artifacts/runs/regression/20260926T132707217977Z-5b4d52c8/outputs/placed-oof.parquet",
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
    "Audit fold and molecular-series uncertainty of 3A4 TDI CatBoost placement",
    "--config",
    str(root / "configs/regression/3a4-tdi-catboost-placement-audit-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_3a4_tdi_catboost_placement_audit.py"],
    cwd=root,
    check=True,
)
