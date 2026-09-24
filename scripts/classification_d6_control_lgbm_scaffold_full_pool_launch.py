"""Record September 15 control LightGBM scaffold refits on complete pools."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
original = root / "artifacts/runs/classification/20260912T115220209511Z-76facc71"
config = root / "configs/analysis/d6-control-lgbm-scaffold-full-pool-v1.json"
cfg = json.loads(config.read_text())
inputs = {
    "compounds.parquet": original / "inputs/000-compounds.parquet",
    "labels.parquet": original / "inputs/001-measurements.parquet",
    "assignments.parquet": original / "inputs/003-assignments.parquet",
    "class-embeddings.parquet": original / "inputs/004-external-embeddings.parquet",
    "mixed-embeddings.parquet": original / "inputs/005-external-embeddings.parquet",
}
for index in range(cfg["mordred_files"]):
    matches = list((original / "inputs").glob(f"*-mordred-{index:05d}.parquet"))
    if len(matches) != 1:
        raise ValueError(f"Expected one frozen Mordred chunk {index}, got {matches}")
    inputs[f"mordred-{index:05d}.parquet"] = matches[0]
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
    "Refit September 15 control LightGBM atoms on complete scaffold folds",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_control_lgbm_all_pool.py"],
    cwd=root,
    check=True,
)
