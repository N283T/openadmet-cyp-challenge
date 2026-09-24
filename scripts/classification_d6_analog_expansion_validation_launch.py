"""Record a D6 pseudo analog-expansion validation subset and comparison."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
classification = root / "artifacts/runs/classification"
config = root / "configs/analysis/d6-analog-expansion-validation-v1.json"
inputs = {
    "compounds.parquet": shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    "regression-labels.parquet": root
    / "artifacts/runs/regression/20260912T074319275987Z-483dcce3/outputs/training-labels.parquet",
    "assignments.parquet": shared / "20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "sep15-umap.parquet": classification
    / "20260924T040301085704Z-914bb06c/outputs/all-pool-oof.parquet",
    "ap-umap.parquet": classification
    / "20260924T035023873825Z-b719ba35/outputs/all-pool-oof.parquet",
    "sep15-scaffold.parquet": classification
    / "20260924T041452692950Z-e17778e3/outputs/all-pool-oof.parquet",
    "ap-scaffold.parquet": classification
    / "20260924T041348966036Z-8db19b4b/outputs/all-pool-oof.parquet",
}
declared = json.loads(config.read_text())
expected = {
    declared["compounds"],
    declared["regression_labels"],
    declared["assignments"],
    *declared["predictions"]["umap"],
    *declared["predictions"]["scaffold"],
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
    "Mimic analog-expansion test selection using non-D6 potent anchors and compare AP with control",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_analog_expansion_validation.py"],
    cwd=root,
    check=True,
)
