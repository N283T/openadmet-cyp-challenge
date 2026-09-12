"""Audit completed single-target controls against the submitted ensemble OOF."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("run")
parser.add_argument("endpoint", choices=["CYP1A2", "CYP2D6"])
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
out = root / f"artifacts/runs/regression/{args.run}/outputs"
config = root / f"configs/focus-single-target/{args.endpoint}_audit_v2_{args.run}.json"
with config.open("x") as stream:
    json.dump(dict(endpoint=args.endpoint, n=1372 if args.endpoint == "CYP1A2" else 1457), stream)
inputs = [
    out / "finetune-oof.parquet",
    root / "artifacts/runs/regression/20260912T072935321335Z-1e8d2643/outputs/external-oof.parquet",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
] + [out / f"fold{fold}/membership.parquet" for fold in range(5)]
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Independent single-target OOF and exclusion audit: " + args.endpoint,
    "--config",
    str(config),
]
for item in inputs:
    command += ["--input", str(item)]
subprocess.run(
    command + ["--", "python", "scripts/audit_focus_single_target.py"], cwd=root, check=True
)
