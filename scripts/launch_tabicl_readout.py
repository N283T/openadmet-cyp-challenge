"""Launch interval target models after distribution inference is finished."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--backend", choices=["tabicl"], default="tabicl")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
parent = root / "artifacts/runs/shared/20260910T132141471698Z-343eb5ec"
base = json.loads((parent / "config.json").read_text())
for fold in range(5):
    cfg = json.loads(json.dumps(base))
    cfg["backend"] = args.backend
    cfg["folds"] = [cfg["folds"][fold]]
    cfg["folds"][0]["endpoints"] = {
        k: v for k, v in cfg["folds"][0]["endpoints"].items() if not k.endswith("is_TDI")
    }
    cfg["weights"].pop("classifier")
    cfg["quantiles"] = [i / 100 for i in range(1, 100)]
    cfg["scope"] = (
        "Train interval-bound marginals and midpoint on existing fold-safe features; "
        "original TabPFN 2.6 architecture/seed"
    )
    if args.backend == "lightgbm":
        cfg["scope"] = "LightGBM L1 on point labels or half-weight interval endpoint copies"
    paths = []

    def remap(obj, paths=paths):
        if isinstance(obj, dict):
            return {k: remap(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [remap(v) for v in obj]
        if isinstance(obj, str) and (parent / "inputs" / obj).is_file():
            p = parent / "inputs" / obj
            if p not in paths:
                paths.append(p)
            return f"{paths.index(p):03d}-{p.name}"
        return obj

    cfg = remap(cfg)
    p = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet"
    cfg["measurements"] = f"{len(paths):03d}-{p.name}"
    paths.append(p)
    p = (
        root
        / "artifacts/runs/shared/20260911T093102177307Z-beddd36f/outputs"
        / "tabicl-regressor-v2-20260212.ckpt"
    )
    cfg["tabicl_weights"] = f"{len(paths):03d}-{p.name}"
    paths.append(p)
    cfg["scope"] = "TabICL matched-feature UMAP comparison; official wrapper defaults, seed 42+fold"
    cfg["official_source_commit"] = "4a8df806d12a7385edb6365074ce218af19002a5"
    prefix = "tabicl-"
    path = root / f"configs/calibration/{prefix}interval-readout-fold{fold}-v1.json"
    with path.open("x") as f:
        f.write(json.dumps(cfg, indent=2) + "\n")
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Train {args.backend} interval-aware regression fold {fold}",
        "--config",
        str(path),
    ]
    for p in paths:
        cmd += ["--input", str(p)]
    subprocess.run(cmd + ["--", "python", "scripts/tabicl_readout.py"], check=True, cwd=root)
