"""Record an identity-matched two-seed legacy/fresh policy comparison."""

import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("legacy_run")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
inputs = {
    "legacy.parquet": runs / args.legacy_run / "outputs/nested-oof.parquet",
    "fresh.parquet": runs / "20261002T084852768236Z-ef8d81ed/outputs/nested-oof.parquet",
    "fresh-cuts.csv": runs / "20261002T084852768236Z-ef8d81ed/outputs/threshold-scale-audit.csv",
    "legacy-result.json": runs / args.legacy_run / "result.json",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Compare legacy and fresh two-seed policies on identical untouched outer molecules",
    "--config",
    str(root / "configs/classification/d6-legacy-nested-review-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/review_cyp_d6_legacy_nested.py"], cwd=root, check=True
)
