"""Record the full-data Caruana deployment component with immutable inputs."""

import argparse
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("--aux", required=True)
a = p.parse_args()
inputs = []


def add(path):
    path = root / path
    inputs.append(path)
    return f"{len(inputs) - 1:03d}-{path.name}"


base = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
cfg = {key: add(f"{base}/{key}.parquet") for key in ["compounds", "measurements", "test_mapping"]}
cfg.update(
    selection=add("artifacts/runs/regression/20260912T054406239560Z-96fce32b/outputs/summary.json"),
    class_emb=add(
        "artifacts/runs/shared/20260912T051543210441Z-a3ea2ec3/outputs/external-embeddings.parquet"
    ),
    mixed_emb=add(
        "artifacts/runs/shared/20260912T052853359792Z-2c0a0953/outputs/external-embeddings.parquet"
    ),
    tree=add(
        "artifacts/runs/shared/20260912T053008533415Z-45e5a1eb/outputs/external-scores.parquet"
    ),
    aux=add(a.aux),
)
config = root / "configs/external-ensemble/full-caruana-deployment-v1.json"
with config.open("x") as f:
    json.dump(cfg, f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    (
        "Refit original Caruana component on all 6525 endpoint labels "
        "using averaged saved inner-selection weights; predict 750 blinded test compounds"
    ),
    "--config",
    str(config),
]
for path in inputs:
    cmd.extend(["--input", str(path)])
subprocess.run(cmd + ["--", "python", "scripts/full_caruana_refit.py"], check=True, cwd=root)
