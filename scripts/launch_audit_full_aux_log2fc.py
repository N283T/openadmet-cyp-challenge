"""Launch independent verification of a completed full auxiliary run."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--run", required=True)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
base = shared / args.run / "outputs"
inputs = []


def add(path):
    inputs.append(path)
    return f"{len(inputs) - 1:03d}-{path.name}"


cfg = {
    "features": add(base / "full-aux.parquet"),
    "model": add(base / "full.pt"),
    "membership": add(base / "membership.parquet"),
    "nested": add(shared / "20260912T052549384296Z-bac5fb3a/outputs/nested-log2fc.parquet"),
    "targets": add(shared / "20260910T100337258849Z-927eac1c/outputs/targets.parquet"),
    "compounds": add(shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"),
}
path = root / f"configs/full-submission/audit-full-aux-{args.run}.json"
with path.open("x") as stream:
    json.dump(cfg, stream, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Independently reload full auxiliary model and verify all train/test feature routing",
    "--config",
    str(path),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/audit_full_aux_log2fc.py"], cwd=root, check=True)
