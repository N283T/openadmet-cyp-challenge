"""Launch coverage or completed prediction audit for the fixed series experiment."""

import argparse
import subprocess
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--scores", action="store_true")
p.add_argument("--run", default="20260911T095536618510Z-f97937d4")
p.add_argument("--config", default="series-cv-v1.json")
args = p.parse_args()
root = Path(__file__).resolve().parents[1]
run = root / "artifacts/runs/regression" / args.run / "outputs"
paths = [
    run / "assignments.parquet",
    run / "hit-seeds.csv",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
]
if args.scores:
    paths += [run / "oof.parquet", root / "artifacts/metric-reference/custom_scoring_functions.py"]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Audit shared series fold integrity and seed-excluded evaluation coverage"
    + (" and model scores" if args.scores else ""),
    "--config",
    str(root / "configs/calibration" / args.config),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/audit_series_cv.py"], check=True, cwd=root)
