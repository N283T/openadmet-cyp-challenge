"""Snapshot native membership artifacts and source tables for independent audit."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
b = root / "artifacts/runs"
data = b / "shared/20260912T104336148000Z-ae7e37fc/outputs"
inputs = [data / "contract.json", data / "challenge-connectivity-folds.parquet"]
runs = []
for run, model in [
    ("20260912T105152218493Z-fec30715", "d6_single1"),
    ("20260912T105211810369Z-423d2caf", "challenge8"),
    ("20260912T105212997766Z-288cae8b", "union18"),
    ("20260912T105214185310Z-e6e93fbf", "d6_isoform6"),
    ("20260912T105256414096Z-03b86209", "d6_single1"),
]:
    output = b / f"regression/{run}/outputs"
    item = dict(run=run, model=model, memberships={})
    for key, path in [
        ("table", data / f"{model}.parquet"),
        ("assignments", output / "assignments.parquet"),
        ("effective", output / "effective-config.json"),
    ]:
        item[key] = len(inputs)
        inputs.append(path)
    for fold in range(5):
        path = output / f"fold{fold}/membership.parquet"
        if path.exists():
            item["memberships"][str(fold)] = len(inputs)
            inputs.append(path)
    runs.append(item)
cfg = root / "configs/community_reproduction/native_data_audit_v2.json"
with cfg.open("x") as f:
    json.dump(dict(runs=runs), f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Independently audit actual native SCP data cohorts and all-head connectivity masks",
    "--config",
    str(cfg),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/community_native_data_audit.py"], cwd=root, check=True
)
