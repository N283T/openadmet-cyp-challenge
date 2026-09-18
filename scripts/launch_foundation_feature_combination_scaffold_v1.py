"""Launch scaffold confirmation for promising foundation feature readouts."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
views = runs / "regression/20260917T111518065286Z-801ab060/outputs"
config = root / "configs/external-training/foundation-feature-combination-scaffold-v1.json"
inputs = {
    "measurements.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "parent-assignments.parquet": runs
    / "regression/20260918T021731984021Z-11c1a26e/outputs/parent-assignments.parquet",
    "custom_scoring_functions.py": root / "artifacts/metric-reference/custom_scoring_functions.py",
    "tabpfn-v3.5-fast-20260909.safetensors": Path.home()
    / ".cache/tabpfn/tabpfn-v3.5-fast-20260909.safetensors",
    "tabicl-regressor-v2-20260212.ckpt": runs
    / "shared/20260911T093102177307Z-beddd36f/outputs/tabicl-regressor-v2-20260212.ckpt",
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
needed = {block["file"] for block in declared["blocks"].values()}
needed.update(
    {
        "measurements.parquet",
        "parent-assignments.parquet",
        declared["scorer"],
        declared["tabpfn_weights"],
        declared["tabicl_weights"],
    }
)
assert needed == set(inputs)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Confirm promising foundation feature readouts on the scaffold secondary split",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/foundation_feature_combination_screen.py"],
    cwd=root,
    check=True,
)
