"""Refit the submitted AP-mass90 TabPFN atom for the full pool."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
views = runs / "regression/20260917T111518065286Z-801ab060/outputs"
config = root / "configs/classification/d6-ap-tabpfn-atom-full-pool-v1.json"
inputs = {
    "measurements.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "test-mapping.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    "assignments.parquet": runs
    / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "current-bank.parquet": runs
    / "classification/20260915T121009912691Z-73198c15/outputs/prediction-bank.parquet",
    "tabpfn-v3.5-fast-20260909.safetensors": Path.home()
    / ".cache/tabpfn/tabpfn-v3.5-fast-20260909.safetensors",
    "fp-mordred.parquet": runs
    / "shared/20260912T061043373484Z-2e344d7a/outputs/fp-mordred-features.parquet",
    "jeremy-raw.parquet": runs
    / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    "jeremy-scores.parquet": runs
    / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
    "anchor-raw-scores.parquet": views / "jeremy_init_task_anchor-raw-scores.parquet",
    "chemeleon-raw.parquet": runs
    / "shared/20260912T061231617816Z-299439e8/outputs/external-embeddings.parquet",
}
declared = json.loads(config.read_text())
needed = {
    declared[key]
    for key in ("labels", "test_mapping", "assignments", "current_bank", "tabpfn_weights")
}
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
    "Refit AP-mass90 TabPFN 3.5 atom with complete-molecule D6 OOF",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_tabpfn35_fast_feature_screen.py"],
    cwd=root,
    check=True,
)
