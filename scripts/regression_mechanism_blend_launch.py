"""Launch the fixed regression mechanism-feature blend audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
config = root / "configs/regression/mechanism-representation-blend-v1.json"
inputs = [
    runs / "regression/20260914T111634758781Z-477d7e4e/outputs/mechanism-oof.parquet",
    runs
    / "regression/20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet",
    runs / "regression/20260912T072935321335Z-1e8d2643/outputs/external-oof.parquet",
    runs / "regression/20260914T111634758781Z-477d7e4e/outputs/memberships.parquet",
    runs
    / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
for path in [config, *inputs]:
    if not path.is_file():
        raise FileNotFoundError(path)

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Audit fixed blends of frozen CYP mechanism readouts with retained regression ensembles",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/regression_mechanism_blend.py"],
    cwd=root,
    check=True,
)
