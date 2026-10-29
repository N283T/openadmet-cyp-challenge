"""Launch the paired CYP3A4 TDI-proxy direct-label augmentation screen."""

import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--config", default="configs/regression/3a4-tdi-proxy-augmentation-v1.json")
parser.add_argument(
    "--weights",
    default="regression/20260917T121017505228Z-433d9acb/inputs/"
    "tabpfn-v3.5-fast-20260909.safetensors",
    help="TabPFN checkpoint relative to artifacts/runs",
)
parser.add_argument("--purpose", default=None)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
views = runs / "regression/20260917T121017505228Z-433d9acb/inputs"
inputs = {
    "measurements.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "parent-assignments.parquet": views / "parent-assignments.parquet",
    "control-3a4.parquet": runs / "regression/20260917T130718808881Z-dfa02730/outputs/oof.parquet",
    "custom_scoring_functions.py": views / "custom_scoring_functions.py",
    "tabpfn.safetensors": runs / args.weights,
    "nested-log2fc.parquet": views / "nested-log2fc.parquet",
    "anchor59-raw-scores.parquet": views / "anchor59-raw-scores.parquet",
    "jeremy-scores.parquet": views / "jeremy-scores.parquet",
    "chemeleon-raw.parquet": views / "chemeleon-raw.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    args.purpose
    or (
        "Paired CYP3A4 screen: add 1,249 shifted TDI-condition pIC50 of TDI-negative "
        "compounds as proxy direct labels to four incumbent-family recipes"
    ),
    "--config",
    str(root / args.config),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/regression_3a4_tdi_proxy.py"], cwd=root, check=True
)
