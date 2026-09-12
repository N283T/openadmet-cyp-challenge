"""Audit both completed kernel comparisons through recorded input snapshots."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
b = root / "artifacts/runs"
k = b / "regression/20260912T101108734911Z-7aa9cf27/outputs"
for kind, run, filename in [
    ("svr", "20260912T101108734911Z-7aa9cf27", "d6-kernel-oof.parquet"),
    ("interval", "20260912T101356027800Z-320cdd53", "d6-interval-kernel-oof.parquet"),
]:
    d = b / "regression" / run / "outputs"
    assert json.loads((d.parent / "result.json").read_text())["status"] == "completed"
    cfg = root / f"configs/d6-kernel/audit_{kind}_v1.json"
    with cfg.open("x") as f:
        json.dump(dict(kind=kind, model_run=run), f, indent=2)
    inputs = [
        d / filename,
        d / "memberships.parquet",
        b / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
        b / "regression/20260912T094021204217Z-84c28266/outputs/global-oof.parquet",
        k / "kernel-ids.parquet",
        d / "inner-trials.csv",
        k / "kernels.npz",
        root / "artifacts/metric-reference/custom_scoring_functions.py",
    ] + sorted(d.glob("*.joblib"))
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "Reconstruct kernel predictions and audit exact inner and outer memberships",
        "--config",
        str(cfg),
    ]
    for p in inputs:
        cmd += ["--input", str(p)]
    subprocess.run(cmd + ["--", "python", "scripts/d6_kernel_audit.py"], cwd=root, check=True)
