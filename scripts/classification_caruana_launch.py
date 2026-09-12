"""Declare the complete currently available matched inner/outer candidate pool."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
run = root / "artifacts/runs/classification/20260912T121146867140Z-b0a05986"
base = json.loads((run / "config.json").read_text())
paths, sources = [], []
for src in base["sources"]:
    item = dict(prefix=src["prefix"])
    for kind in ["oof", "inner", "test"]:
        p = run / "inputs" / src[kind]
        assert p.is_file()
        item[kind] = f"{len(paths):03d}-{p.name}"
        paths.append(p)
    sources.append(item)
cfg = dict(
    sources=sources,
    steps=50,
    recipes=[
        "subset1",
        "subset2",
        "subset3",
        "caruana_mcc",
        "caruana_mcc_sparse3",
        "caruana_auc",
        "caruana_brier",
    ],
)
config = root / "configs/classification-round/caruana-v1.json"
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
    "Compare single sparse subsets and Caruana selection on isolated inner predictions",
    "--config",
    str(config),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/classification_caruana.py"], cwd=root, check=True)
