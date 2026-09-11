"""Declare fine-tuning OOF and fixed reference ensemble inputs."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
r = re.findall(r"run_id=(\S+)", (root / "artifacts/launches/admetica-finetune.log").read_text())
assert len(r) == 1
files = [
    root / f"artifacts/runs/regression/{r[0]}/outputs/finetune-oof.parquet",
    root / "artifacts/runs/regression/20260911T150745346553Z-6a171ac2/outputs/aligned-oof.parquet",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
cfg = dict(
    model_run=r[0],
    database=dict(host=str(root / ".runtime/socket"), port=5434, dbname="cyp_challenge"),
)
p = root / "configs/parent-cv/admetica-ft-audit.json"
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
    "Audit matched direct fine-tuning and adapted embedding TabPFN predictions and fixed blends",
    "--config",
    str(p),
]
for f in files:
    cmd += ["--input", str(f)]
subprocess.run(cmd + ["--", "python", "scripts/audit_admetica_ft.py"], check=True, cwd=root)
