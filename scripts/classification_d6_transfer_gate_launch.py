"""Record conservative transfer evidence decision for the AP replay."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/analysis/d6-ap-transfer-gate-v1.json"
inputs = {
    "umap-fixed-rate.csv": runs / "20260924T040346490330Z-638c65b7/outputs/fixed-rate.csv",
    "umap-rate-selection.csv": runs / "20260924T040346490330Z-638c65b7/outputs/comparison.csv",
    "scaffold-fixed-rate.csv": runs / "20260924T041524000098Z-99cb9161/outputs/fixed-rate.csv",
    "scaffold-rate-selection.csv": runs / "20260924T041524000098Z-99cb9161/outputs/comparison.csv",
}
declared = json.loads(config.read_text())
if set(inputs) != {
    declared[f"{strategy}_{mode}"]
    for strategy in declared["strategies"]
    for mode in ("fixed", "nested")
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
    "Apply complete-pool dual-split evidence gate to AP versus September 15 replay",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_transfer_gate.py"],
    cwd=root,
    check=True,
)
