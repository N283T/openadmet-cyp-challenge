"""Declare every partial/completed run record used for a complete legacy comparison."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
cfg_path = root / "configs/classification/d6-legacy-partition-merge-v1.json"
cfg = json.loads(cfg_path.read_text())
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Verify and merge all 400 predetermined legacy ensemble fits from four compute partitions",
    "--config",
    str(cfg_path),
]
for i, spec in enumerate(cfg["partitions"]):
    run = root / "artifacts/runs/classification" / spec["run_id"]
    for name in ["config.json", "result.json", "manifest.json"]:
        path = run / name
        if not path.is_file():
            raise FileNotFoundError(path)
        command += ["--input", f"p{i}-{name}={path}"]
    for name in [
        "nested-oof.parquet",
        "threshold-scale-audit.csv",
        "fit-audit.csv",
        "fold-metrics.csv",
        "leakage-and-operating-points.csv",
    ]:
        path = run / "outputs" / name
        if not path.is_file():
            raise FileNotFoundError(path)
        command += ["--input", f"p{i}-{name}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/merge_cyp_d6_legacy_partitions.py"], cwd=root, check=True
)
