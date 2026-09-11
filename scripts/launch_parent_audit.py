"""Declare all five paired-readout runs for the final parent-retention audit."""

import argparse
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--readout-log", default="artifacts/launches/parent-cv-readout-20260911.log")
parser.add_argument("--revision", default="")
args = parser.parse_args()
paths = []


def add(p):
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


parents = re.findall(r"run_id=(\S+)", (root / args.readout_log).read_text())
assert len(parents) == 5
cfg = {
    "oof": [
        add(root / f"artifacts/runs/regression/{r}/outputs/parent-oof.parquet") for r in parents
    ]
}
cfg["database"] = json.loads((root / "configs/calibration/package-interval-v1.json").read_text())[
    "database"
]
cfg["assignments"] = add(
    root
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
)
cfg["labels"] = add(
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet"
)
cfg["scorer"] = add(root / "artifacts/metric-reference/custom_scoring_functions.py")
p = root / f"configs/parent-cv/audit{args.revision}.json"
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
    "Audit corrected parent-retained full-feature LightGBM against matched parent-excluded control",
    "--config",
    str(p),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/audit_parent_cv.py"], check=True, cwd=root)
