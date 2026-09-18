"""Launch full-fit inference for new v8 feature-combination members."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
views = runs / "regression/20260917T111518065286Z-801ab060/outputs"
config = root / "configs/external-training/fullfit-feature-combination-members-v8.json"
inputs = {
    "measurements.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "aliases.parquet": runs / "shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
    "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
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
    "full-aux.parquet": runs / "shared/20260912T074339967571Z-ff194fcf/outputs/full-aux.parquet",
}
declared = json.loads(config.read_text())
needed = {block["file"] for block in declared["blocks"].values()}
needed.update(
    {
        "measurements.parquet",
        "aliases.parquet",
        "reference.csv",
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
    "Full-fit new v8 LGBM, TabPFN and TabICL atomic members for blinded-test inference",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/fullfit_feature_combination_members.py"],
    cwd=root,
    check=True,
)
