"""Record the complete-pool D6 control call-rate robustness audit."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/classification/d6-control-rate-robustness-20260926.json"
tracked = root / "artifacts/tracked/d6-ap-full-pool-replay-20260924"
inputs = {
    "umap-control.parquet": tracked / "sep15-control-all-pool-oof.parquet",
    "scaffold-control.parquet": tracked / "sep15-control-scaffold-all-pool-oof.parquet",
}
declared = json.loads(config.read_text())
if set(declared["oof"].values()) != set(inputs):
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
    "Audit D6 September 15 control call rate across complete-pool UMAP and scaffold OOF",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_call_rate_robustness.py"],
    cwd=root,
    check=True,
)
