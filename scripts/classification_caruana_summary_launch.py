"""Launch recorded Caruana complexity and paired evaluation report."""

import argparse
import json
import subprocess
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("run")
a = p.parse_args()
root = Path(__file__).resolve().parents[1]
run = root / "artifacts/runs/classification" / a.run / "outputs"
paths = []


def add(path):
    assert path.is_file(), path
    paths.append(path)
    return f"{len(paths) - 1:03d}-{path.name}"


cfg = {
    k: add(run / f)
    for k, f in [
        ("oof", "oof.parquet"),
        ("deployment", "deployment.parquet"),
        ("selections", "selections.json"),
        ("metrics", "metrics.csv"),
    ]
}
cfg["reference"] = add(
    root
    / "artifacts/runs/classification/20260912T121146867140Z-b0a05986/outputs/ensemble-oof.parquet"
)
cfg["labels"] = add(
    root
    / "artifacts/runs/classification/20260912T115018743045Z-269e598d/outputs"
    / "classification-labels.parquet"
)
config = root / "configs/classification-round" / f"caruana-summary-{a.run}.json"
assert not config.exists()
config.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Summarize Caruana complexity and paired outer MCC versus fixed reference means",
    "--config",
    str(config),
]
for path in paths:
    cmd += ["--input", str(path)]
subprocess.run(
    cmd + ["--", "python", "scripts/classification_caruana_summary.py"], check=True, cwd=root
)
