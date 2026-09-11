"""Audit completed embedding-only fits against existing matched controls."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = re.findall(r"run_id=(\S+)", (root / "artifacts/launches/pubchem600-readout.log").read_text())
assert len(runs) == 1
r = runs[0]
files = [
    root / f"artifacts/runs/regression/{r}/outputs/pubchem600-oof.parquet",
    root / "artifacts/runs/regression/20260911T150745346553Z-6a171ac2/outputs/aligned-oof.parquet",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
cfg = dict(
    model_run=r,
    database=dict(host=str(root / ".runtime/socket"), port=5434, dbname="cyp_challenge"),
)
p = root / "configs/parent-cv/pubchem600-audit.json"
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
    "Compare PubChem600 embedding-only TabPFN OOF and fixed blends with current models",
    "--config",
    str(p),
]
for f in files:
    cmd += ["--input", str(f)]
subprocess.run(cmd + ["--", "python", "scripts/audit_pubchem600.py"], check=True, cwd=root)
