"""Record independent full-data Caruana verification."""

import argparse
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("--run", required=True)
a = p.parse_args()
inputs = []


def add(path):
    path = root / path
    inputs.append(path)
    return f"{len(inputs) - 1:03d}-{path.name}"


base = f"artifacts/runs/regression/{a.run}/outputs"
shared = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
cfg = {
    "source_selection": add(
        "artifacts/runs/regression/20260912T054406239560Z-96fce32b/outputs/summary.json"
    ),
    "actual_selection": add(f"{base}/summary.json"),
    "measurements": add(f"{shared}/measurements.parquet"),
    "membership": add(f"{base}/membership.parquet"),
    "predictions": add(f"{base}/test-predictions.parquet"),
    "test_mapping": add(f"{shared}/test_mapping.parquet"),
}
config = root / f"configs/external-ensemble/full-caruana-audit-{a.run}.json"
with config.open("x") as f:
    json.dump(cfg, f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    (
        "Independently verify complete-label membership and saved selection transfer "
        "for full Caruana deployment"
    ),
    "--config",
    str(config),
]
for path in inputs:
    cmd.extend(["--input", str(path)])
subprocess.run(cmd + ["--", "python", "scripts/audit_full_caruana_refit.py"], check=True, cwd=root)
