"""Recorded 2D6 fingerprint kernel comparison."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
r = root / "artifacts/runs"
cfg = root / "configs/d6-kernel/v1.json"
with cfg.open("x") as f:
    json.dump(
        dict(
            C=[0.3, 1.0, 3.0, 10.0, 30.0],
            epsilon=[0.05, 0.2],
            split_run="20260911T104703655337Z-a3923255",
        ),
        f,
        indent=2,
    )
inputs = [
    r / "shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    r / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    r / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    r / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    r / "regression/20260912T094021204217Z-84c28266/outputs/global-oof.parquet",
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
    "2D6 grouped inner-selected count-fingerprint kernel SVR and fixed global additions",
    "--config",
    str(cfg),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/d6_kernel.py"], check=True, cwd=root)
