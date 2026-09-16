"""Register and launch the second CYP2D6 ranking-diversity screen."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = [
    runs / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    runs / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    runs / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
    runs / "shared/20260912T061043373484Z-2e344d7a/outputs/fp-mordred-features.parquet",
    runs / "shared/20260912T082554377951Z-d272774c/outputs/nested-log2fc.parquet",
    runs / "regression/20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet",
    runs / "regression/20260912T151747998390Z-f24622ac/outputs/crossfit-predictions.parquet",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)

config = root / "configs/regression/d6-ranking-kernel-classical-screen-v1.json"
payload = json.loads(config.read_text())
config.write_text(json.dumps(payload, indent=2) + "\n")

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    (
        "Compare uncertainty-weighted pairwise, Tanimoto-kernel, ordinal, additive-spline, "
        "and robust CYP2D6 readouts under retained parent CV"
    ),
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_ranking_kernel_classical_screen.py"],
    cwd=root,
    check=True,
)
