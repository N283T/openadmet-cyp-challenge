"""Regenerate the two fixed CYP3A4 LightGBM recipes as separate OOF columns."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/regression/caruana-3a4-internal-pruning-v1.json"
source_inputs = root / "artifacts/runs/regression/20260912T152832952478Z-3e4de7c8/inputs"
cfg = json.loads(config.read_text())
aliases = {
    cfg["assignments"],
    cfg["measurements"],
    cfg["compounds"],
    cfg["class_emb"],
    cfg["mixed_emb"],
    cfg["tree"],
    cfg["scorer"],
}
aliases.update(item["file"] for item in cfg["components"])
for fold in cfg["folds"]:
    aliases.add(fold["selection"])
    aliases.add(fold["nested"])
    aliases.update(fold["models"].values())

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Export fixed CYP3A4 class and Morgan quantile LightGBMs as atomic OOF members",
    "--config",
    str(config),
]
for alias in sorted(aliases):
    command.extend(["--input", f"{alias}={source_inputs / alias}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_caruana_3a4_pruning.py"],
    cwd=root,
    check=True,
)
