"""Declare the focused kernel/tree comparison without changing evaluation folds."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = []


def add(path):
    path = root / path
    assert path.is_file()
    name = f"{len(inputs):03d}-{path.name}"
    inputs.append(path)
    return name


cfg = dict(
    endpoints=[f"{ep}_pIC50_direct_inhibition" for ep in ["CYP1A2", "CYP2D6"]],
    views=["scores_svr", "medium_svr", "combined_lgbm"],
    trials=16,
    split_run="20260911T104703655337Z-a3923255",
    assignments=add(
        "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    measurements=add(
        "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet"
    ),
    aux=add("artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs/nested-log2fc.parquet"),
    baseline=add(
        "artifacts/runs/regression/20260912T072935321335Z-1e8d2643/outputs/external-oof.parquet"
    ),
    medium=add(
        "artifacts/runs/shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet"
    ),
    scores=add(
        "artifacts/runs/shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet"
    ),
    classic=add(
        "artifacts/runs/shared/20260912T061043373484Z-2e344d7a/outputs/fp-mordred-features.parquet"
    ),
)
config = root / "configs/regression/focus-readout-screen-v1.json"
with config.open("x") as f:
    json.dump(cfg, f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Focused1A2/2D6: grouped inner-selected RBF-SVR and combined external features LGBM",
    "--config",
    str(config),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/focus_readout_screen.py"], cwd=root, check=True)
