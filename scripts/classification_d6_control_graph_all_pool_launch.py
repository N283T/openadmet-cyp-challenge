"""Record the frozen TDI-only graph control atom on complete UMAP folds."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
original = root / "artifacts/runs/classification/20260912T120135164747Z-1da37202"
config = root / "configs/analysis/d6-control-graph-all-pool-v1.json"
cfg = json.loads(config.read_text())
inputs = {
    "compounds.parquet": original / "inputs/000-compounds.parquet",
    "labels.parquet": original / "inputs/001-measurements.parquet",
    "assignments.parquet": original / "inputs/002-assignments.parquet",
    "external-best.pt": original / "inputs/005-best.pt",
}
for fold, alias in enumerate(cfg["refits"]):
    inputs[alias] = original / f"outputs/tdi_only-fold{fold}/refit.pt"
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
    "Extend frozen September 15 TDI-only graph fits to complete UMAP folds",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_control_graph_all_pool.py"],
    cwd=root,
    check=True,
)
