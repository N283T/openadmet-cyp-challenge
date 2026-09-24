"""Record the three frozen Morgan control atoms on complete UMAP folds."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
original = root / "artifacts/runs/classification/20260913T123045515718Z-0df22d54"
config = root / "configs/analysis/d6-control-fixed-fp-all-pool-v1.json"
cfg = json.loads(config.read_text())
inputs = {
    "compounds.parquet": original / "inputs/000-compounds.parquet",
    "labels.parquet": original / "inputs/001-measurements.parquet",
    "assignments.parquet": original / "inputs/003-assignments.parquet",
}
for method in cfg["methods"]:
    for fold in range(5):
        inputs[f"model-{method}-fold{fold}.joblib"] = (
            original / f"outputs/{method}-fold{fold}.joblib"
        )
for path in [config, *inputs.values()]:
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Extend frozen September 15 Morgan tail models to complete UMAP folds",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_control_fixed_fp_all_pool.py"],
    cwd=root,
    check=True,
)
