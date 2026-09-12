"""Launch jointgraph audit with every final scaler declared separately."""

import argparse
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
b = root / "artifacts/runs"
parser = argparse.ArgumentParser()
parser.add_argument("--run", default="20260912T101831048657Z-478527aa")
parser.add_argument("--outer-refit", action="store_true")
args = parser.parse_args()
model = args.run
o = b / f"regression/{model}/outputs"
cfg = root / f"configs/d6_joint_external_audit_{model}.json"
with cfg.open("x") as f:
    json.dump(dict(model_run=model, outer_refit=args.outer_refit), f)
inputs = [
    o / "joint-oof.parquet",
    o / "memberships.parquet",
    o / "external-membership.parquet",
    o / "encoder-gradient-audit.csv",
    b / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    b / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    b / "shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    b / "shared/20260912T101212878172Z-616de5d8/outputs/efficacy-targets.parquet",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
inputs += [o / f"p{w:02d}-fold{fold}/scaling.json" for fold in range(5) for w in [0, 5, 30]]
if args.outer_refit:
    inputs += [
        o / "progress.json",
        b / "regression/20260912T101831048657Z-478527aa/outputs/progress.json",
    ]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Audit joint graph identity splits scalers labelweights gradients and OOF",
    "--config",
    str(cfg),
]
for p in inputs:
    cmd.extend(["--input", str(p)])
subprocess.run(cmd + ["--", "python", "scripts/d6_joint_external_audit.py"], cwd=root, check=True)
