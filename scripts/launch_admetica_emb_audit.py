"""Audit completed embedding-only fits against existing matched controls."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = re.findall(r"run_id=(\S+)", (root / "artifacts/launches/admetica-emb-only.log").read_text())
assert len(runs) == 1
r = runs[0]
files = [
    root / f"artifacts/runs/regression/{r}/outputs/admetica-emb-oof.parquet",
    root / "artifacts/runs/regression/20260911T142905405713Z-d0f742da/outputs/aligned-oof.parquet",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
cfg = dict(
    model_run=r,
    database=dict(host=str(root / ".runtime/socket"), port=5434, dbname="cyp_challenge"),
)
p = root / "configs/parent-cv/admetica-emb-audit.json"
with p.open("x") as stream:
    json.dump(cfg, stream, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Compare Admetica embedding-only TabPFN OOF and fixed blends with current models",
    "--config",
    str(p),
]
for f in files:
    cmd += ["--input", str(f)]
subprocess.run(cmd + ["--", "python", "scripts/audit_admetica_emb.py"], check=True, cwd=root)
