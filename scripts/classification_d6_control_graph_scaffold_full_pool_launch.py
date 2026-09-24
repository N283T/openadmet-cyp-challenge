"""Record scaffold refits of the September 15 TDI-only graph atom."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
original = root / "artifacts/runs/classification/20260912T120135164747Z-1da37202"
config = root / "configs/analysis/d6-control-graph-scaffold-full-pool-v1.json"
inputs = {
    "compounds.parquet": original / "inputs/000-compounds.parquet",
    "labels.parquet": original / "inputs/001-measurements.parquet",
    "assignments.parquet": root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "external-best.pt": original / "inputs/005-best.pt",
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
    "Refit September 15 TDI-only graph atom on complete scaffold folds",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_control_graph_scaffold_full_pool.py"],
    cwd=root,
    check=True,
)
