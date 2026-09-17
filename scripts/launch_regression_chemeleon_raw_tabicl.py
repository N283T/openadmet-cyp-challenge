"""Launch the matched raw CheMeleon2048 plus nested-log2fc TabICL comparison."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/regression/chemeleon-raw2048-log2fc-tabicl-v1.json"
inputs = {
    "assignments.parquet": root
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet",
    "measurements.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "tabicl-regressor.ckpt": root
    / "artifacts/runs/shared/20260911T093102177307Z-beddd36f/outputs"
    / "tabicl-regressor-v2-20260212.ckpt",
    "embeddings.parquet": root
    / "artifacts/runs/shared/20260912T061231617816Z-299439e8/outputs"
    / "external-embeddings.parquet",
    "nested-log2fc.parquet": root
    / "artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs"
    / "nested-log2fc.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    (
        "Compare raw CheMeleon2048 plus nested log2fc4 with the retained PCA256 "
        "TabICL component on fixed parent CV"
    ),
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/external_tabicl_readout.py"],
    cwd=root,
    check=True,
)
