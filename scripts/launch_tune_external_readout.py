"""Launch nested external feature LGBM tuning with immutable inputs."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = [
    root
    / "artifacts/runs/shared/20260912T051543210441Z-a3ea2ec3/outputs/external-embeddings.parquet",
    root / "artifacts/runs/shared/20260912T051543210441Z-a3ea2ec3/outputs/summary.json",
    root / "artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs/nested-log2fc.parquet",
    root / "artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs/summary.json",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    root
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
cfg = dict(
    seed=42,
    n_trials=12,
    inner_fraction=0.2,
    method="external3_class_logfc_tuned_lgbm",
    split_run="20260911T104703655337Z-a3923255",
    embedding_run="20260912T051543210441Z-a3ea2ec3",
    nested_run="20260912T052549384296Z-bac5fb3a",
    selection="inner grouped ST_RAE only; anchor groups always innertrain",
    threads=4,
)
config = root / "configs/external-tuning/class-logfc-lgbm.json"
config.parent.mkdir(parents=True, exist_ok=True)
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
    "Nested Rustuna LightGBM tuning of external classification embeddings and safe log2fc readout",
    "--config",
    str(config),
]
for p in paths:
    cmd.extend(["--input", str(p)])
subprocess.run(cmd + ["--", "python", "scripts/tune_external_readout.py"], check=True, cwd=root)
