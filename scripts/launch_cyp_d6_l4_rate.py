"""Launch the D6 L4 re-thresholding package."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
package = runs / "classification/20261002T223203072149Z-662a599c/inputs"
inputs = {
    "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
    "test_mapping.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    "sep15_blind_scores.parquet": package / "sep15_blind_scores.parquet",
    "test_scores.parquet": package / "test_scores.parquet",
    "last_sent.parquet": root
    / "artifacts/submissions/classification/20261002T223229553466Z-a08301d2/predictions.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Re-threshold the deployed L4 rank-mean D6 score at 95 calls (A4 unchanged, no refit)",
    "--config",
    str(root / "configs/classification/d6-l4-rate-v1.json"),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/package_cyp_d6_l4_rate.py"], cwd=root, check=True
)
