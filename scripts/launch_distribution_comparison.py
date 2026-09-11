"""Audit completed distribution readouts with explicit parent input maps."""

import argparse
import json
import re
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--combined", action="store_true")
parser.add_argument("--revision", default="v2")
args = parser.parse_args()
paths = []


def add(p):
    p = Path(p)
    if p not in paths:
        paths.append(p)
    return f"{paths.index(p):03d}-{p.name}"


cfg = dict(
    seed=20260911,
    repeats=5,
    estimators=["mean", "median", "q0.40", "q0.45", "q0.50", "q0.55", "q0.60", "q0.65", "q0.70"],
    baseline=[],
    additional=[],
)
for key, log, oname, tname in [
    (
        "baseline",
        "distribution-readout-20260911.log",
        "distribution-oof.parquet",
        "fold_predictions.parquet",
    ),
    ("additional", "blend-inference-20260911.log", "member-oof.parquet", "member-test.parquet"),
]:
    if key == "additional" and not args.combined:
        continue
    runs = re.findall(r"^run_id=(\S+)", Path("artifacts/launches", log).read_text(), flags=re.M)
    assert len(runs) == len(set(runs)) == 5
    for run in runs:
        parent = Path("artifacts/runs/regression", run)
        assert json.loads((parent / "result.json").read_text())["status"] == "completed"
        cfg[key].append(
            dict(run=run, oof=add(parent / "outputs" / oname), test=add(parent / "outputs" / tname))
        )
for p in [
    "artifacts/runs/regression/20260911T050944039369Z-36d00d9c/outputs/regression-oof-diagnostics.parquet",
    "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "artifacts/metric-reference/custom_scoring_functions.py",
]:
    add(p)
name = "combined" if args.combined else "baseline"
path = Path(f"configs/calibration/distribution-comparison-{name}-{args.revision}.json")
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
    f"Compare {name} distribution summaries and held-out calibration; no external labels",
    "--config",
    str(path),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/compare_distributions.py"], check=True)
