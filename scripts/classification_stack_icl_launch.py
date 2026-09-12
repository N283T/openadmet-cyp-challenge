"""Launch after explicit shared-GPU coordination, from frozen CPU-prepared matrices."""

import argparse
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
b = root / "artifacts/runs"
a = argparse.ArgumentParser()
a.add_argument("--prepared-run", required=True)
args = a.parse_args()
source = b / "classification" / args.prepared_run / "outputs"
inputs = []


def add(p):
    inputs.append(p)
    return f"{len(inputs) - 1:03d}-{p.name}"


cfg = dict(
    weights=add(
        b / "shared/20260912T134236392057Z-8c08c1de/outputs/tabicl-classifier-v2-20260212.ckpt"
    ),
    test=add(source / "test-identities.parquet"),
    prepared_run=args.prepared_run,
    matrices=[],
)
for item in json.loads((source / "matrices.json").read_text()):
    cfg["matrices"].append(
        dict(
            item,
            features=add(source / item["features"]),
            membership=add(source / item["membership"]),
        )
    )
config = root / "configs/classification-stack-icl/fit.json"
config.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Historical stack matched TabICL classification: only downstream classifier replacement,20fits",
    "--config",
    str(config),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/classification_stack_icl.py"], cwd=root, check=True)
