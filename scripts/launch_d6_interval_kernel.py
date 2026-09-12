"""Launch recorded matched interval kernel regression."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
b = root / "artifacts/runs"
k = b / "regression/20260912T101108734911Z-7aa9cf27/outputs"
cfg = root / "configs/d6-kernel/interval_v2.json"
with cfg.open("x") as f:
    json.dump(
        dict(
            kernels=["fingerprint_mix", "fp_emb_mix"],
            regularization=[0.0003, 0.003, 0.03],
            split_run="20260911T104703655337Z-a3923255",
        ),
        f,
        indent=2,
    )
files = [
    k / "kernels.npz",
    k / "kernel-ids.parquet",
    b / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    b / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    b / "regression/20260912T094021204217Z-84c28266/outputs/global-oof.parquet",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "2D6 matched point versus interval kernel regression; inner regularization selection",
    "--config",
    str(cfg),
]
for p in files:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/d6_interval_kernel.py"], cwd=root, check=True)
