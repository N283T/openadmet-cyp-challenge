"""Declare the fixed data and model inputs for auxiliary HPO."""

import argparse
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--revision", default="")
args = parser.parse_args()
shared = root / "artifacts/runs/shared"
cfg = {
    "n_trials": 16,
    "base": json.loads((root / "configs/parent-cv/aux-0-1.json").read_text()),
    "selection": "Inner validation standardized MSE only; outer 0 and 1 excluded",
}
p = root / f"configs/parent-cv/rustuna-log2fc{args.revision}.json"
with p.open("x") as stream:
    json.dump(cfg, stream, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Rustuna 16-trial log2fc HPO including original control; hold folds 0 and 1 untouched",
    "--config",
    str(p),
]
for path in [
    shared / "20260910T100337258849Z-927eac1c/outputs/targets.parquet",
    shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    root / "artifacts/runs/regression/20260911T100134004950Z-c7a3f2bf/outputs/assignments.parquet",
    shared / "20260910T012612848389Z-bfeed434/outputs/chemeleon_mp.pt",
    shared / "20260910T012612848389Z-bfeed434/outputs/provenance.json",
]:
    cmd += ["--input", str(path)]
subprocess.run(cmd + ["--", "python", "scripts/tune_log2fc.py"], check=True, cwd=root)
