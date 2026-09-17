"""Launch honest nested raw600 LightGBM comparison of owned and Jeremy views."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
JOINT_SOURCE = ROOT / "artifacts/runs/regression/20260916T143736207347Z-cc6eb40b/outputs"
JOINT_BASE = ROOT / "artifacts/runs/regression/20260916T143051880581Z-7c464f30/outputs"
OLDER = ROOT / "artifacts/runs/regression/20260916T125409299702Z-7e362fe6/inputs"
CONFIG = ROOT / "configs/external-training/nested-raw-lgbm-owned-vs-jeremy-v1.json"

inputs = {
    "measurements.parquet": CHALLENGE / "measurements.parquet",
    "parent-assignments.parquet": ASSIGNMENTS,
    "custom_scoring_functions.py": ROOT / "artifacts/metric-reference/custom_scoring_functions.py",
    "owned-joint-continuous.parquet": JOINT_SOURCE
    / "trial-03-joint_all_no_binary/external-embeddings.parquet",
    "owned-no-polaris.parquet": OLDER / "source-no-polaris.parquet",
    "owned-joint-cyp.parquet": JOINT_BASE
    / "trial-02-joint_compact_cyp/external-embeddings.parquet",
    "owned-proportional.parquet": OLDER / "objective-proportional-mixed.parquet",
    "jeremy.parquet": ROOT
    / "artifacts/runs/regression/20260912T055849013223Z-f16a5dfd/inputs"
    / "000-external-embeddings.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Nested target-specific raw LightGBM comparison of owned and Jeremy embeddings",
    "--config",
    str(CONFIG),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/nested_raw_lgbm_readout.py"],
    cwd=ROOT,
    check=True,
)
