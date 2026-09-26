"""Record the submitted A4 equal-seven full-pool composition."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/classification/a4-equal7-full-pool-20260927.json"
inputs = {
    "historical-all-pool.parquet": runs
    / "20260926T150611714864Z-d10195fa/outputs/all-pool-oof.parquet",
    "jeremy-all-pool.parquet": runs
    / "20260926T150342241560Z-da6a7d98/outputs/all-pool-oof.parquet",
    "trees-all-pool.parquet": runs / "20260926T150058222462Z-cf35fc2e/outputs/all-pool-oof.parquet",
    "equal7-selections.parquet": runs
    / "20260915T125001009276Z-d5adff22/outputs/selections.parquet",
    "equal7-labelled-oof.parquet": runs / "20260915T125001009276Z-d5adff22/outputs/oof.parquet",
    "equal7-test.parquet": runs / "20260915T125001009276Z-d5adff22/outputs/test.parquet",
    "submitted-test-scores.parquet": runs
    / "20260915T125939829893Z-1d921ffc/outputs/test-scores.parquet",
}
declared = json.loads(config.read_text())
if set(inputs) != {
    declared[key]
    for key in (
        "historical",
        "jeremy",
        "trees",
        "selection",
        "saved_labelled_oof",
        "saved_test",
        "submitted_test",
    )
}:
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
    "Compose the exact submitted A4 equal-seven on the full UMAP molecular pool",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
completed = subprocess.run(
    command + ["--", "python", "scripts/classification_a4_equal7_full_pool.py"],
    cwd=root,
    check=False,
)
if completed.returncode:
    raise SystemExit(f"A4 equal-seven composition failed with exit code {completed.returncode}")
