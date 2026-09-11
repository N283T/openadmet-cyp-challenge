"""Declare immutable assets for shared-series CV comparison."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--broad", action="store_true")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
base = root / "artifacts/runs/shared"
paths = [
    base / "20260909T132734119233Z-a7063ccd/outputs" / n
    for n in ["compounds.parquet", "aliases.parquet", "measurements.parquet"]
]
paths += [base / "20260909T135130268493Z-edba6fc0/outputs/assignments.parquet"]
paths += sorted((base / "20260909T132918829974Z-bc81af58/outputs").glob("mordred-*.parquet"))
paths += [
    Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
if args.broad:
    paths.append(
        root / "artifacts/runs/regression/20260911T095536618510Z-f97937d4/outputs/oof.parquet"
    )
cfg = root / (
    "configs/calibration/series-cv-broad-v1.json"
    if args.broad
    else "configs/calibration/series-cv-v1.json"
)
with cfg.open("x") as f:
    json.dump(
        dict(
            similarity_threshold=0.4 if args.broad else 0.5,
            strategies=["series"] if args.broad else ["umap", "series"],
            seeds_per_enzyme=25,
            folds=5,
            seed=42,
            models=["lightgbm_l1_300", "tabpfn_2.6_8"],
            features="Mordred, outer-training median imputation",
        ),
        f,
        indent=2,
    )
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Compare series and UMAP CV with matched Mordred LightGBM and TabPFN",
    "--config",
    str(cfg),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/series_cv.py"], check=True, cwd=root)
