"""Launch the molecular-support-only D6 random split search."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
inputs = {
    "compounds.parquet": runs / "compounds.parquet",
    "measurements.parquet": runs / "measurements.parquet",
    "test-mapping.parquet": runs / "test_mapping.parquet",
}
for path in inputs.values():
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Search label-blind CYP2D6 random folds closest to blind-test molecular support",
    "--config",
    str(root / "configs/classification/d6-random-split-search-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_random_split_search.py"],
    cwd=root,
    check=True,
)
