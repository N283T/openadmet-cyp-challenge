"""Launch CatBoost and TabPFN readouts for complementary checkpoint views."""

import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VIEWS = ROOT / "artifacts/runs/regression/20260917T111518065286Z-801ab060/outputs"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
config = ROOT / "configs/external-training/distilled-view-diverse-readout-v1.json"
inputs = {
    "measurements.parquet": CHALLENGE / "measurements.parquet",
    "parent-assignments.parquet": ASSIGNMENTS,
    "custom_scoring_functions.py": ROOT
    / "artifacts/metric-reference/custom_scoring_functions.py",
    "tabpfn-v2.6-regressor-v2.6_default.ckpt": Path.home()
    / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt",
}
for candidate in json.loads(config.read_text())["candidates"]:
    inputs[candidate["file"]] = VIEWS / candidate["file"]
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Test CatBoost and TabPFN on complementary raw and external-task predictor views",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/screen_checkpoint_views.py"],
    cwd=ROOT,
    check=True,
)
