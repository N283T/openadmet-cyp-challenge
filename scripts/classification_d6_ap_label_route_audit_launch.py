"""Record the AP/control TDI label-route audit on both split families."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/analysis/d6-ap-label-route-audit-v1.json"
inputs = {
    "sep15-umap.parquet": runs / "20260924T040301085704Z-914bb06c/outputs/all-pool-oof.parquet",
    "ap-umap.parquet": runs / "20260924T035023873825Z-b719ba35/outputs/all-pool-oof.parquet",
    "sep15-scaffold.parquet": runs / "20260924T041452692950Z-e17778e3/outputs/all-pool-oof.parquet",
    "ap-scaffold.parquet": runs / "20260924T041348966036Z-8db19b4b/outputs/all-pool-oof.parquet",
    "assignments.parquet": root
    / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "assay-route.parquet": runs / "20260912T125901882388Z-4d006374/outputs/assay-route.parquet",
}
declared = json.loads(config.read_text())
expected = {
    *declared["predictions"]["umap"],
    *declared["predictions"]["scaffold"],
    declared["assignments"],
    declared["assay_route"],
}
if set(inputs) != expected:
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
    "Audit assay-label provenance for AP versus September 15 full-pool D6 decisions",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_ap_label_route_audit.py"],
    cwd=root,
    check=True,
)
