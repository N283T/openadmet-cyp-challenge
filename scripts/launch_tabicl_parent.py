"""Snapshot parent-retained inputs/preprocessing and run matched TabICL baselines."""

import argparse
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--fold", type=int)
parser.add_argument("--revision", default="")
parser.add_argument("--backend", choices=["tabicl", "tabpfn"], default="tabicl")
parser.add_argument("--update-predictions", action="store_true")
parser.add_argument("--readout-log", help="Explicit historical readout log")
parser.add_argument("--readout-registry", default="configs/regression/current-parent-readout.json")
args = parser.parse_args()
if args.readout_log:
    parents = re.findall(r"run_id=(\S+)", (root / args.readout_log).read_text())
else:
    registry = json.loads((root / args.readout_registry).read_text())
    parents = [item["run"] for item in sorted(registry["models"], key=lambda item: item["fold"])]
assert len(parents) == 5
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
    if args.backend == "tabpfn":
        cfg["weights"] = add(Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt")
    cfg["backend"] = args.backend
    cfg["update_predictions"] = args.update_predictions
    if args.update_predictions:
        cfg["feature_sets"] = ["full"]
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
    cfg["model"] = args.backend + "; fixed checkpoint; 8 estimators; mean"
    cfg["features"] = "compact PCA256+predicted-log2fc4 and full +Mordred"
    p = root / f"configs/parent-cv/{args.backend}-{fold}{args.revision}.json"
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
        f"Update {args.backend} pIC50 model on adopted auxiliary features, fold {fold}",
        "--config",
        str(p),
    ]
    for p in paths:
        cmd += ["--input", str(p)]
    subprocess.run(cmd + ["--", "python", "scripts/tabicl_parent_readout.py"], check=True, cwd=root)
