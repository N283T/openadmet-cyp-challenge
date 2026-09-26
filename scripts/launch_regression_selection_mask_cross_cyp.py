"""Record synthetic D6-style label-selection tests on other CYP endpoints."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
    "parent-assignments.parquet": (
        runs / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    "scaffold-assignments.parquet": (
        runs / "regression/20260918T021731984021Z-11c1a26e/outputs/parent-assignments.parquet"
    ),
    "measurements.parquet": (
        runs / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet"
    ),
    "regression-labels.parquet": (
        runs / "regression/20260926T004020214692Z-b2f253b0/inputs/regression-labels.parquet"
    ),
    "scores.parquet": (
        runs / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet"
    ),
    "medium.parquet": (
        runs / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet"
    ),
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Simulate D6 selective labels on other CYPs and test cross-CYP weak supervision",
    "--config",
    str(root / "configs/regression/selection-mask-cross-cyp-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_selection_mask_cross_cyp.py"],
    cwd=root,
    check=True,
)
