"""Launch leakage-safe TDI probability readouts for CYP3A4 regression."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REG = ROOT / "artifacts/runs/regression"
CLS = ROOT / "artifacts/runs/classification"
inputs = {
    "regression-oof.parquet": REG / "20260916T062147025088Z-1d1208a1/outputs/aligned-oof.parquet",
    "tdi-oof.parquet": CLS / "20260915T125001009276Z-d5adff22/outputs/oof.parquet",
    "tdi-member-oof.parquet": CLS / "20260915T125001009276Z-d5adff22/inputs/000-member-oof.parquet",
    "scorer.py": ROOT / "artifacts/metric-reference/custom_scoring_functions.py",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Cross-fit CYP3A4 pIC50 readouts from the rank-1 TDI probability bank",
    "--config",
    str(ROOT / "configs/regression/tdi-a4-probability-readout-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/screen_tdi_a4_probability_readouts.py"],
    cwd=ROOT,
    check=True,
)
