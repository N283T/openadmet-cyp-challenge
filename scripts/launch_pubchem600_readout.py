"""Snapshot completed PubChem600 features for matched downstream evaluation."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
r = re.findall(r"run_id=(\S+)", (root / "artifacts/launches/pubchem600-export.log").read_text())
assert len(r) == 1
files = [
    root / f"artifacts/runs/shared/{r[0]}/outputs/pubchem600-embeddings.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    root
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet",
]
weights = Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt"
cfg = dict(
    pretrain_run=r[0],
    split_run="20260911T104703655337Z-a3923255",
    weights=f"{len(files):03d}-{weights.name}",
)
files.append(weights)
p = root / "configs/parent-cv/pubchem600-readout.json"
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
    "Evaluate frozen PubChem600 embeddings alone using TabPFN2.6 on matched folds",
    "--config",
    str(p),
]
for f in files:
    cmd += ["--input", str(f)]
subprocess.run(cmd + ["--", "python", "scripts/pubchem600_readout.py"], check=True, cwd=root)
