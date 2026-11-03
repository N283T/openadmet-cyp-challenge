"""Launch the legacy plus MONROE-pretrained D6 blind call package."""

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
    "pool-labels.parquet": runs
    / "classification/20261002T215345688167Z-7920b49a/outputs/oof_scores.parquet",
    "embeddings.parquet": runs
    / "regression/20261029T115010020552Z-7bae1535/inputs/embeddings.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Blind D6 calls from legacy Sep 15 plus a MONROE-pretrained pool LightGBM; A4 unchanged",
    "--config",
    str(root / "configs/classification/d6-swap-monroe-v1.json"),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/package_cyp_d6_swap_monroe.py"], cwd=root, check=True
)
