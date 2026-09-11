"""Audit the five completed interval-target model runs."""

import argparse
import json
import re
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--include-lightgbm", action="store_true")
args = parser.parse_args()
runs = re.findall(
    r"^run_id=(\S+)",
    Path("artifacts/launches/interval-readout-20260911.log").read_text(),
    flags=re.M,
)
assert len(runs) == len(set(runs)) == 5
if args.include_lightgbm:
    extra = re.findall(
        r"^run_id=(\S+)",
        Path("artifacts/launches/lightgbm-interval-readout-20260911.log").read_text(),
        flags=re.M,
    )
    assert len(extra) == len(set(extra)) == 5
    runs += extra
paths = []


def add(p):
    p = Path(p)
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


cfg = dict(seed=20260911, runs=[], neighbor_method="knn(50, 2, 'interval'):weight1.0")
for run in runs:
    parent = Path("artifacts/runs/regression", run)
    assert json.loads((parent / "result.json").read_text())["status"] == "completed"
    cfg["runs"].append(
        dict(
            run=run,
            oof=add(parent / "outputs/interval-oof.parquet"),
            test=add(parent / "outputs/interval-test.parquet"),
        )
    )
for p in [
    "artifacts/runs/regression/20260911T073458027748Z-50193fa7/outputs/neighborhood-oof.parquet",
    "artifacts/runs/regression/20260911T073458027748Z-50193fa7/outputs/neighborhood-test.parquet",
    "artifacts/runs/regression/20260911T050944039369Z-36d00d9c/outputs/regression-oof-diagnostics.parquet",
    "artifacts/runs/regression/20260911T072717606484Z-8ec7593a/outputs/distribution-comparison-oof.parquet",
    "artifacts/runs/regression/20260911T072717606484Z-8ec7593a/outputs/distribution-comparison-test.parquet",
    "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "artifacts/metric-reference/custom_scoring_functions.py",
]:
    add(p)
version = "combined-v2" if args.include_lightgbm else "v1"
path = Path(f"configs/calibration/interval-audit-{version}.json")
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
    "Audit interval-marginal and midpoint predictions plus fixed blends on all UMAP folds",
    "--config",
    str(path),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/audit_interval_readout.py"], check=True)
