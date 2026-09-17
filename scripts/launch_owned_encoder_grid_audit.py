"""Audit completed project-owned encoder readouts as independent CYP2D6 members."""

import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "artifacts/launches/owned_encoder_grid_v1_manifest.json"
GRID = json.loads(
    (ROOT / "configs/external-training/owned-encoder-grid-v1.json").read_text()
)
members = json.loads(MANIFEST.read_text())
assert members
config = {
    "endpoint": GRID["endpoints"][0],
    "dataset_run": GRID["dataset_run"],
    "protocol": "Matched independent atomic scores only; no blending or weight search.",
    "scorer": "custom_scoring_functions.py",
    "variants": [
        {"name": member["name"], "file": f"{member['name']}.parquet"} for member in members
    ],
}
config_path = ROOT / "configs/owned-encoder-grid/owned_encoder_grid_d6_audit.json"
config_path.write_text(json.dumps(config, indent=2) + "\n")
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Matched CYP2D6 audit of project-owned encoder width and raw/PCA atomic readouts",
    "--config",
    str(config_path),
    "--input",
    "custom_scoring_functions.py="
    + str(ROOT / "artifacts/metric-reference/custom_scoring_functions.py"),
]
for member in members:
    command.extend(["--input", f"{member['name']}.parquet={member['oof']}"])
command.extend(["--", "python", "scripts/audit_atomic_readout_grid.py"])
subprocess.run(command, cwd=ROOT, check=True)
