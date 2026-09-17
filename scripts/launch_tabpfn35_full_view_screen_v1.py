"""Launch full TabPFN-3.5 on target-winning and complementary views."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs"
VIEWS = RUNS / "regression/20260917T111518065286Z-801ab060/outputs"
inputs = {
    "measurements.parquet": RUNS
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "parent-assignments.parquet": RUNS
    / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    "custom_scoring_functions.py": ROOT
    / "artifacts/metric-reference/custom_scoring_functions.py",
    "tabpfn-v3.5-20260909.safetensors": Path.home()
    / ".cache/tabpfn/tabpfn-v3.5-20260909.safetensors",
    "nested-log2fc.parquet": RUNS
    / "shared/20260912T052549384296Z-bac5fb3a/outputs/nested-log2fc.parquet",
    "jeremy-raw.parquet": RUNS
    / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    "jeremy-scores.parquet": RUNS
    / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
    "anchor59-raw-scores.parquet": VIEWS / "jeremy_init_task_anchor-raw-scores.parquet",
    "pub030-raw.parquet": VIEWS / "jeremy_pubchem_anchor030-raw.parquet",
    "chemeleon-raw.parquet": RUNS
    / "shared/20260912T061231617816Z-299439e8/outputs/external-embeddings.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Evaluate full TabPFN-3.5 on fast-screen target winners and a raw control",
    "--config",
    str(ROOT / "configs/external-training/tabpfn35-full-view-screen-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/screen_checkpoint_views.py"],
    cwd=ROOT,
    check=True,
)
