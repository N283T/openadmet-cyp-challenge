"""Record full-fit D6 upper-bound inference and package a test-vector probe."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
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
    "aliases.parquet": (runs / "shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet"),
    "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
    "current-submission.csv": (
        root / "artifacts/submissions/regression/20260924T054343781258Z-db4d3cf3/predictions.csv"
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
    "Full-fit selected D6 one-sided-primary model and inspect blind vector delta",
    "--config",
    str(root / "configs/regression/d6-primary-bound-fullfit-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_primary_bound_fullfit.py"],
    cwd=root,
    check=True,
)
