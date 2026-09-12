"""Declare fixed-ensemble source and result tables for independent audit."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("run")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
run = root / "artifacts/runs/classification" / args.run
inputs = [run / "config.json", *sorted((run / "inputs").iterdir())]
for name in [
    "oof.parquet",
    "deployment.parquet",
    "outer-weighted-fold-test.parquet",
    "metrics.csv",
    "selections.json",
    "full-oof-fitting-only.csv",
]:
    inputs.append(run / "outputs" / name)
canonical = (
    root
    / "artifacts/runs/classification/20260912T115018743045Z-269e598d/outputs"
    / "classification-labels.parquet"
)
config = root / "configs/classification" / f"caruana-audit-{args.run}.json"
config.write_text(
    json.dumps(
        dict(
            config="000-config.json",
            canonical=f"{len(inputs):03d}-{canonical.name}",
            source_run=args.run,
        ),
        indent=2,
    )
)
inputs.append(canonical)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Independent reconstruction of classification Caruana selection and predictions",
    "--config",
    str(config),
]
for path in inputs:
    assert path.is_file(), path
    cmd += ["--input", str(path)]
subprocess.run(
    cmd + ["--", "python", "scripts/classification_caruana_audit.py"], check=True, cwd=root
)
