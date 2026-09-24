"""Record complete-molecule D6 OOF for two static-feature readouts."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
config = root / "configs/classification/d6-full-pool-static-transfer-v1.json"
inputs = {
    "measurements.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "test-mapping.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    "assignments.parquet": runs
    / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "current-bank.parquet": runs
    / "classification/20260915T121009912691Z-73198c15/outputs/prediction-bank.parquet",
    "jeremy-raw.parquet": runs
    / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    "jeremy-scores.parquet": runs
    / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
    "chemeleon-raw.parquet": runs
    / "shared/20260912T061231617816Z-299439e8/outputs/external-embeddings.parquet",
}
declared = json.loads(config.read_text())
needed = {declared[key] for key in ("labels", "test_mapping", "assignments", "current_bank")}
needed.update(spec["file"] for spec in declared["blocks"].values())
if set(inputs) != needed:
    raise ValueError("Config and input aliases differ")
for path in [config, *inputs.values()]:
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
    "Refit two static D6 readouts and predict every UMAP/scaffold held-fold molecule",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_regression_feature_transfer.py"],
    cwd=root,
    check=True,
)
