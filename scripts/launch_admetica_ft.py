"""Declare assets for direct and embedding fine-tuning comparisons."""

import json
import subprocess
from pathlib import Path

base = Path("artifacts/runs")
files = sorted((base / "shared/20260911T144903295084Z-ac8cd6ac/outputs").glob("*-encoder.pt"))
files += [
    base / "shared/20260909T132734119233Z-a7063ccd/outputs" / n
    for n in ["compounds.parquet", "measurements.parquet"]
]
files += [base / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"]
weights = Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt"
cfg = dict(
    weights=f"{len(files):03d}-{weights.name}",
    split_run="20260911T104703655337Z-a3923255",
    encoder_lr=1e-5,
    head_lr=1e-4,
    weight_decay=1e-4,
    max_epochs=60,
    patience=10,
    loss="inner-training standardized pIC50 MSE",
    encoder="matched CYP Admetica",
    downstream="TabPFN2.6 8 estimators on adapted emb300",
)
files.append(weights)
p = Path("configs/parent-cv/admetica-finetune.json")
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
    "Fine-tune Admetica on pIC50; evaluate direct heads and adapted-embedding TabPFN",
    "--config",
    str(p),
]
for f in files:
    cmd += ["--input", str(f)]
subprocess.run(cmd + ["--", "python", "scripts/admetica_finetune.py"], check=True)
