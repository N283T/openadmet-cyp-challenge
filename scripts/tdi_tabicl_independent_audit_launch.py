"""Declare all frozen TabICL fit states, features, and matched PFN probabilities."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("run")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
run = root / "artifacts/runs/classification" / args.run
inputs = [
    run / "config.json",
    *sorted((run / "inputs").iterdir()),
    *sorted((run / "outputs").glob("*")),
]
pfn = (
    root
    / "artifacts/runs/classification/20260912T115242271820Z-51023dac/outputs"
    / "oof-probabilities.parquet"
)
pfn_inner = (
    root
    / "artifacts/runs/classification/20260912T115749418593Z-397ddee7/outputs"
    / "inner-validation-probabilities.parquet"
)
canonical = (
    root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet"
)
cfg = dict(
    config="000-config.json",
    run=args.run,
    pfn=f"{len(inputs):03d}-{pfn.name}",
    pfn_inner=f"{len(inputs) + 1:03d}-{pfn_inner.name}",
    canonical=f"{len(inputs) + 2:03d}-{canonical.name}",
)
inputs += [pfn, pfn_inner, canonical]
config = root / "configs/classification" / f"tabicl-independent-{args.run}.json"
config.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Independent TabICL feature and fit-state isolation audit with matched PFN "
        "threshold comparisons"
    ),
    "--config",
    str(config),
]
for p in inputs:
    assert p.is_file(), p
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/tdi_tabicl_independent_audit.py"], cwd=root, check=True
)
