"""Reuse exact frozen Morgan/Mordred source inputs for diagnostic training controls."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
base = root / "artifacts/runs/classification/20260912T115220209511Z-76facc71"
old = json.loads((base / "config.json").read_text())
inputs = []


def add(p):
    inputs.append(p)
    return f"{len(inputs) - 1:03d}-{p.name}"


cfg = {key: add(base / "inputs" / old[key]) for key in ["compounds", "labels", "assignments"]}
cfg["mordred"] = [add(base / "inputs" / p) for p in old["mordred"]]
cfg["baseline"] = add(base / "outputs/oof.parquet")
cfg["recipes"] = {
    "strong_regularization": {"trees": 100, "leaves": 7, "min_child": 80, "lambda": 10},
    "train_label_permutation": {"trees": 400, "leaves": 15, "min_child": 20, "lambda": 2},
}
cfg.update(seed="42+fold", learning_rate=0.03, class_weight="balanced", n_jobs=3)
config = root / "configs/classification-round/d6-sanity.json"
config.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "D6 strong regularization and train-label permutation: ten fixed sanity-control fits",
    "--config",
    str(config),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/classification_d6_sanity.py"], cwd=root, check=True)
