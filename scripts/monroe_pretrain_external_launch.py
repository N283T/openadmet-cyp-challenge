"""Launch continued MONROE pretraining on the external CYP corpus."""

import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("graph_run", help="shared run holding the external graphs-seed1.pkl")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
checkpoint = Path.home() / "ghq/github.com/blazejba/monroe/checkpoint"
inputs = {
    "graphs-seed1.pkl": shared / args.graph_run / "outputs/graphs-seed1.pkl",
    "observations.parquet": shared
    / "20260916T082302554619Z-656c63e6/outputs/observations-compact-cyp.parquet",
    "monroe-config.json": checkpoint / "config.json",
    "monroe-weights.pt": checkpoint / "weights.pt",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Continue MONROE pretraining on 32 challenge-disjoint external CYP assays (regression + calls)",
    "--config",
    str(root / "configs/shared/monroe-pretrain-external-v1.json"),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
python = str(root / ".pixi/envs/monroe/bin/python")
subprocess.run(
    command + ["--", python, "scripts/monroe_pretrain_external.py"], cwd=root, check=True
)
