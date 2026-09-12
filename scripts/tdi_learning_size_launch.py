"""Record fixed source features and exact full-model controls for learning-size study."""

import json
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
run = root / "artifacts/runs/classification/20260912T115220209511Z-76facc71"
original = json.loads((run / "config.json").read_text())
paths = []


def add(path):
    assert path.is_file(), path
    name = f"{len(paths):03d}-{path.name}"
    paths.append(path)
    return name


cfg = {key: add(run / "inputs" / original[key]) for key in ["compounds", "labels", "assignments"]}
cfg["mordred"] = [add(run / "inputs" / name) for name in original["mordred"]]
cfg["reference"] = add(run / "outputs/oof.parquet")
cfg["full_models"] = {
    f"{ep}-{fold}": add(run / "outputs" / f"fp_mordred_balanced-{ep}-{fold}.txt")
    for ep in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]
    for fold in range(5)
}
cfg["seeds"] = [17, 42, 73]
paired = "--paired" in sys.argv
if paired:
    cfg["route_ids"] = add(
        root
        / "artifacts/runs/classification/20260912T125901882388Z-4d006374/outputs"
        / "assay-route.parquet"
    )
    cfg["route_column"], cfg["route_value"] = "availability", "both"
config = (
    root
    / "configs/classification-round"
    / ("learning-size-paired-v1.json" if paired else "learning-size-v1.json")
)
config.write_text(json.dumps(cfg, indent=2))
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Fixed-validation A4 learning size matched to D6 and paired shared-population endpoint study",
    "--config",
    str(config),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/tdi_learning_size.py"], cwd=root, check=True)
