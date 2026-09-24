"""Launch fold-safe transfer of regression feature pools to D6 classification."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
views = runs / "regression/20260917T111518065286Z-801ab060/outputs"
config = root / "configs/classification/regression-feature-transfer-v1.json"
inputs = {
    "measurements.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "test-mapping.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    "assignments.parquet": runs
    / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "current-prediction-bank.parquet": runs
    / ("classification/20260915T121009912691Z-73198c15/outputs/prediction-bank.parquet"),
    "fp-mordred.parquet": runs
    / "shared/20260912T061043373484Z-2e344d7a/outputs/fp-mordred-features.parquet",
    "jeremy-raw.parquet": runs
    / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    "jeremy-scores.parquet": runs
    / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
    "anchor-raw-scores.parquet": views / "jeremy_init_task_anchor-raw-scores.parquet",
    "pub030-raw.parquet": views / "jeremy_pubchem_anchor030-raw.parquet",
    "owned-raw-scores.parquet": views / "owned_init_hybrid-raw-scores.parquet",
    "scratch-raw.parquet": views / "scratch_teacher_only-raw.parquet",
}
declared = json.loads(config.read_text())
needed = {spec["file"] for spec in declared["blocks"].values()}
needed.update(
    {
        declared["labels"],
        declared["test_mapping"],
        declared["assignments"],
        declared["current_bank"],
    }
)
if needed != set(inputs):
    raise ValueError(
        f"Input mismatch: missing={needed - set(inputs)}, extra={set(inputs) - needed}"
    )
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
    (
        "Transfer recent regression checkpoint and heterogeneous feature pools "
        "into fold-safe CYP2D6 classifiers on UMAP and scaffold splits"
    ),
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
