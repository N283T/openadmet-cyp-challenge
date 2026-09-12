"""Record the fixed four-member native prediction comparison."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
for name in ["single1", "aux8", "union18", "isoform6"]:
    parser.add_argument(f"--{name}", required=True, help="Canonical export run ID")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
reg = root / "artifacts/runs/regression"
paths = [
    reg / "20260912T094021204217Z-84c28266/outputs/global-oof.parquet",
    reg / "20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
members = []
for name in ["single1", "aux8", "union18", "isoform6"]:
    run = getattr(args, name)
    endpoints = (
        ["CYP2D6"] if name in {"single1", "isoform6"} else ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
    )
    members.append(
        dict(
            name=name,
            export_run=run,
            oof_input=len(paths),
            test_input=len(paths) + 1,
            endpoints=[f"{ep}_pIC50_direct_inhibition" for ep in endpoints],
        )
    )
    paths.extend(
        [
            reg / run / "outputs/scp-oof.parquet",
            reg / run / "outputs/test-predictions.parquet",
        ]
    )
cfg = dict(
    members=members,
    validation_description=(
        "Fixed parent folds with held-out connectivity excluded from all heads; "
        "native held-fold val_loss checkpoint selection; fixed equal direct means"
    ),
)
config = root / "configs/scp-repro/completed-ensemble-evaluation.json"
with config.open("x") as f:
    json.dump(cfg, f, indent=2)
    f.write("\n")
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Evaluate four source-native Chemprop members and their fixed direct mean",
    "--config",
    str(config),
]
for path in paths:
    assert path.is_file(), path
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/scp_repro_evaluate.py"], cwd=root, check=True)
