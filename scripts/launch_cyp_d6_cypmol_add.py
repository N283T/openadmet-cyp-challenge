"""Launch the CYPMol-add test for the D6 TDI score."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
    "testlike-pool-oof.parquet": runs
    / "classification/20261002T215345688167Z-7920b49a/outputs/oof_scores.parquet",
    "assignments.parquet": runs
    / "classification/20261002T221918697661Z-2d8d20f5/inputs/assignments.parquet",
    "analog-pairs.parquet": root
    / "artifacts/tracked/d6-analog-expansion-validation-20260924/anchor-analog-pairs.parquet",
    "embeddings.parquet": runs
    / "shared/20261029T142347849581Z-2c6a74f4/outputs/external-embeddings.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Add frozen CYPMol D6 vectors (LightGBM and logistic pool models) to the deployed L4 "
        "D6 score "
        "and compare with paired grouped bootstraps on Metric O, Metric M and the analog pool"
    ),
    "--config",
    str(root / "configs/classification/d6-cypmol-add-v1.json"),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/evaluate_cyp_d6_cypmol_add.py"], cwd=root, check=True
)
