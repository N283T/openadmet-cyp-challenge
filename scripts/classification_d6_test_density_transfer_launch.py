"""Record fixed train-to-blind density-ratio transfer diagnostic."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/analysis/d6-test-density-transfer-v1.json"
inputs = {
    "compounds.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    "test-mapping.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    "ap-umap.parquet": runs / "20260924T035023873825Z-b719ba35/outputs/all-pool-oof.parquet",
    "control-umap.parquet": runs / "20260924T040301085704Z-914bb06c/outputs/all-pool-oof.parquet",
    "ap-scaffold.parquet": runs / "20260924T041348966036Z-8db19b4b/outputs/all-pool-oof.parquet",
    "control-scaffold.parquet": runs
    / "20260924T041452692950Z-e17778e3/outputs/all-pool-oof.parquet",
}
declared = json.loads(config.read_text())
if set(inputs) != {
    declared[key]
    for key in [
        "compounds",
        "test_mapping",
        "ap_umap",
        "control_umap",
        "ap_scaffold",
        "control_scaffold",
    ]
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
    "Diagnose test-domain density weighting for frozen AP versus September 15 control",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_test_density_transfer.py"],
    cwd=root,
    check=True,
)
