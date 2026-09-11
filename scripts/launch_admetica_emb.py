"""Record checkpoint, fold and label sources for the embedding-only run."""

import json
import subprocess
from pathlib import Path

base = Path("artifacts/runs")
files = sorted((base / "shared/20260911T144903295084Z-ac8cd6ac/outputs").glob("*-encoder.pt"))
assert len(files) == 4
files += [
    base / "shared/20260909T132734119233Z-a7063ccd/outputs" / n
    for n in ["compounds.parquet", "measurements.parquet"]
]
files += [base / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"]
weights = Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt"
cfg = dict(
    weights=f"{len(files):03d}-{weights.name}",
    split_run="20260911T104703655337Z-a3923255",
    encoder_source="20260911T144903295084Z-ac8cd6ac",
    features="matching CYP raw mean-pool300 only",
    model="TabPFN2.6",
    n_estimators=8,
)
files.append(weights)
p = Path("configs/parent-cv/admetica-emb-only.json")
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
    "Evaluate frozen per-CYP Admetica embeddings alone with TabPFN 2.6 on matched folds",
    "--config",
    str(p),
]
for f in files:
    cmd += ["--input", str(f)]
subprocess.run(cmd + ["--", "python", "scripts/admetica_emb_tabpfn.py"], check=True)
