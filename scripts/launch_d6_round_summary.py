"""Summarize a completed 2D6 selector from declared artifacts."""

import argparse
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--tag", required=True)
args = parser.parse_args()
m = json.loads((root / f"artifacts/launches/{args.tag}_manifest.json").read_text())
b = root / "artifacts/runs/regression"
selection = b / m["selection_run"] / "outputs"
audit = b / m["audit_run"] / "outputs"
prior = b / "20260912T094021204217Z-84c28266/outputs"
cfg = root / f"configs/global-caruana/{args.tag}_summary_v2.json"
with cfg.open("x") as f:
    json.dump(dict(selection_run=m["selection_run"], audit_run=m["audit_run"]), f, indent=2)
inputs = [
    audit / "metrics.csv",
    b / "20260912T094029976646Z-846181c4/outputs/metrics.csv",
    selection / "weights.json",
    selection / "catalog.json",
    selection / "summary.json",
    prior / "global-oof.parquet",
    selection / "aligned-pool.parquet",
    selection / "global-oof.parquet",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Summarize completed 2D6 focused pool against frozen previous ensemble",
    "--config",
    str(cfg),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/d6_round_summary.py"], cwd=root, check=True)
