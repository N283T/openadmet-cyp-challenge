"""Include the historical stack with freshly reconstructed honest inner predictions."""

import argparse
import json
import subprocess
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("historical_run")
a = p.parse_args()
root = Path(__file__).resolve().parents[1]
base = root / "artifacts/runs/classification/20260912T121146867140Z-b0a05986"
history = root / "artifacts/runs/classification" / a.historical_run / "outputs"
prior = json.loads((base / "config.json").read_text())
paths = []


def add(path):
    assert path.is_file(), path
    paths.append(path)
    return f"{len(paths) - 1:03d}-{path.name}"


sources = [
    dict(
        prefix=s["prefix"],
        **{kind: add(base / "inputs" / s[kind]) for kind in ["oof", "inner", "test"]},
    )
    for s in prior["sources"]
]
sources.append(
    dict(
        prefix="history",
        oof=add(history / "outer-oof.parquet"),
        inner=add(history / "inner-probabilities.parquet"),
        test=add(history / "fold-test-probabilities.parquet"),
    )
)
cfg = dict(
    sources=sources,
    steps=50,
    bagging_seed=731,
    bagging_samples=20,
    historical_run=a.historical_run,
    recipes=[
        "fixed_historical",
        "subset1",
        "subset2",
        "subset3",
        "caruana_mcc",
        "caruana_mcc_sparse3",
        "caruana_auc",
        "caruana_brier",
        "bagged_sparse3",
    ],
)
config = root / "configs/classification-round/caruana-v2.json"
assert not config.exists()
config.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Compare18 classifier candidates including historical stack, sparse and bagged Caruana",
    "--config",
    str(config),
]
for path in paths:
    cmd += ["--input", str(path)]
subprocess.run(cmd + ["--", "python", "scripts/classification_caruana.py"], check=True, cwd=root)
