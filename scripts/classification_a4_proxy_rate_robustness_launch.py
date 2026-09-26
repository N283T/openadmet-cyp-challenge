"""Record a complete-pool A4 proxy call-rate robustness audit."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/classification/a4-proxy-rate-robustness-20260926.json"
proxy = root / (
    "artifacts/runs/classification/20260926T144248722886Z-fe484e91/outputs/all-pool-oof.parquet"
)
inputs = {"umap-proxy.parquet": proxy, "scaffold-proxy.parquet": proxy}
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
    "Audit A4 positive-call rates across complete-pool UMAP and scaffold proxy OOF",
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
