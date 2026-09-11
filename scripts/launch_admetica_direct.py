"""Snapshot direct-model application inputs."""

import subprocess
from pathlib import Path

repo = Path("/home/nagaet/ghq/github.com/datagrok-ai/admetica")
base = Path("artifacts/runs")
files = [
    repo / f"ADMET/metabolism/{e}-inhibitor/{e}-inhibitor.ckpt"
    for e in ["cyp1a2", "cyp2c9", "cyp2d6", "cyp3a4"]
]
files += [
    base / "shared/20260909T132734119233Z-a7063ccd/outputs" / n
    for n in ["compounds.parquet", "measurements.parquet", "test_mapping.parquet"]
]
files += [
    base / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    base / "regression/20260911T142905405713Z-d0f742da/outputs/aligned-oof.parquet",
]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Apply unchanged Admetica CYP models and compare probabilities with pIC50",
]
for p in files:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/apply_admetica_direct.py"], check=True)
