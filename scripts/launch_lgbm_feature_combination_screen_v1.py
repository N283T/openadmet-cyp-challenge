"""Launch the broad heterogeneous LightGBM feature-combination screen."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
views = runs / "regression/20260917T111518065286Z-801ab060/outputs"
config = root / "configs/external-training/lgbm-feature-combination-screen-v1.json"
inputs = {
    "measurements.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "parent-assignments.parquet": runs
    / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    "custom_scoring_functions.py": root
    / "artifacts/metric-reference/custom_scoring_functions.py",
    "fp-mordred.parquet": runs
    / "shared/20260912T061043373484Z-2e344d7a/outputs/fp-mordred-features.parquet",
    "nested-log2fc.parquet": runs
    / "shared/20260912T052549384296Z-bac5fb3a/outputs/nested-log2fc.parquet",
    "jeremy-raw.parquet": runs
    / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    "jeremy-scores.parquet": runs
    / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
    "chemeleon-raw.parquet": runs
    / "shared/20260912T061231617816Z-299439e8/outputs/external-embeddings.parquet",
    "anchor-raw-scores.parquet": views / "jeremy_init_task_anchor-raw-scores.parquet",
    "pub030-raw.parquet": views / "jeremy_pubchem_anchor030-raw.parquet",
    "owned-raw-scores.parquet": views / "owned_init_hybrid-raw-scores.parquet",
    "scratch-raw.parquet": views / "scratch_teacher_only-raw.parquet",
}
declared = json.loads(config.read_text())
needed = {block["file"] for block in declared["blocks"].values()}
needed.update({"measurements.parquet", "parent-assignments.parquet", declared["scorer"]})
assert needed == set(inputs)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    (
        "Screen heterogeneous molecular, checkpoint, score and auxiliary feature "
        "combinations with fixed LightGBM readouts"
    ),
    "--config",
    str(config),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/lgbm_feature_combination_screen.py"],
    cwd=root,
    check=True,
)
