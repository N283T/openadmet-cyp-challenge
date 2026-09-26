"""Record whether frozen features can identify weak D6 inhibition."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = {
    "parent-assignments.parquet": (
        runs / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    "scaffold-assignments.parquet": (
        runs / "regression/20260918T021731984021Z-11c1a26e/outputs/parent-assignments.parquet"
    ),
    "measurements.parquet": (
        runs / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet"
    ),
    "scores.parquet": (
        runs / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet"
    ),
    "medium.parquet": (
        runs / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet"
    ),
}
command = [
    "pixi", "run", "cyp", "run", "--track", "regression",
    "--purpose", "Test fold-safe identification of weak CYP2D6 activity on both splits",
    "--config", str(root / "configs/regression/d6-weak-activity-classifier-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_weak_activity_classifier.py"],
    cwd=root,
    check=True,
)
