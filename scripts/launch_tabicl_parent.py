"""Snapshot parent-retained inputs/preprocessing and run matched TabICL baselines."""

import argparse
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parents = re.findall(
    r"run_id=(\S+)", (root / "artifacts/launches/parent-cv-readout-20260911.log").read_text()
)
assert len(parents) == 5
parser = argparse.ArgumentParser()
parser.add_argument("--fold", type=int)
parser.add_argument("--revision", default="")
args = parser.parse_args()
for fold, r in enumerate(parents):
    if args.fold is not None and fold != args.fold:
        continue
    parent = root / f"artifacts/runs/regression/{r}"
    cfg = json.loads((parent / "config.json").read_text())
    paths = []

    def add(p, paths=paths):
        paths.append(p)
        return f"{len(paths) - 1:03d}-{p.name}"

    def remap(obj, parent=parent):
        if isinstance(obj, dict):
            return {k: remap(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [remap(v) for v in obj]
        if isinstance(obj, str) and (parent / "inputs" / obj).is_file():
            return add(parent / "inputs" / obj)
        return obj

    # Auxiliary maps are explicit; other shared inputs are matched by filename.
    cfg = remap(cfg)
    for p in sorted((parent / "inputs").iterdir()):
        if p not in paths:
            add(p)
    # Nested snapshot names still end in the canonical feature filenames.
    cfg["weights"] = add(
        root
        / "artifacts/runs/shared/20260911T093102177307Z-beddd36f/outputs"
        / "tabicl-regressor-v2-20260212.ckpt"
    )
    cfg["preprocessing"] = {}
    for ep in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]:
        e = ep + "_pIC50_direct_inhibition"
        out = parent / "outputs"
        cfg["preprocessing"][e] = {
            k: add(out / f"{e}-parent_retained-{suffix}")
            for k, suffix in [
                ("pca", "pca.joblib"),
                ("imputer", "imputer.joblib"),
                ("columns", "columns.json"),
                ("booster", "full.txt"),
            ]
        }
    add(parent / "outputs/parent-oof.parquet")
    cfg["model"] = "TabICL 2.2.0; v2 fixed checkpoint; 8 estimators; mean"
    cfg["features"] = "compact PCA256+predicted-log2fc4 and full +Mordred"
    p = root / f"configs/parent-cv/tabicl-{fold}{args.revision}.json"
    with p.open("x") as f:
        json.dump(cfg, f, indent=2)
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"TabICL compact and full, parent CV fold {fold}; verify matched features",
        "--config",
        str(p),
    ]
    for p in paths:
        cmd += ["--input", str(p)]
    subprocess.run(cmd + ["--", "python", "scripts/tabicl_parent_readout.py"], check=True, cwd=root)
