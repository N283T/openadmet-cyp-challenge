"""Record scaffold refits of the three September 15 Morgan tail atoms."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
original = root / "artifacts/runs/classification/20260913T123045515718Z-0df22d54"
config = root / "configs/analysis/d6-control-fixed-fp-scaffold-full-pool-v1.json"
inputs = {
    "compounds.parquet": original / "inputs/000-compounds.parquet",
    "labels.parquet": original / "inputs/001-measurements.parquet",
    "assignments.parquet": original / "inputs/003-assignments.parquet",
}
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
    "Refit September 15 Morgan tail atoms on complete scaffold folds",
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
