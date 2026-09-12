"""Declare old/new primary inputs for the fixed five-model replacement comparison."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = []


def add_run(run):
    directory = root / f"artifacts/runs/regression/{run}"
    assert json.loads((directory / "result.json").read_text())["status"] == "completed"
    path = directory / "outputs/external-oof.parquet"
    assert path.is_file()
    name = f"{len(inputs):03d}-{path.name}"
    inputs.append(path)
    return name


cfg = dict(
    reference=add_run("20260912T072935321335Z-1e8d2643"),
    original_components=[
        add_run(r)
        for r in [
            "20260912T054406239560Z-96fce32b",
            "20260912T060059705105Z-2f0fbfbd",
            "20260912T072303592078Z-07f2749a",
            "20260912T072056445264Z-e6b6cbcb",
            "20260912T071619253798Z-73a6fb7e",
        ]
    ],
    replacements={
        "2": add_run("20260912T083209287272Z-e8487aeb"),
        "3": add_run("20260912T083421950497Z-036b9595"),
        "4": add_run("20260912T083210446408Z-8d59a1f4"),
    },
    addition=add_run("20260912T083010595992Z-1e7001de"),
)
scorer = root / "artifacts/metric-reference/custom_scoring_functions.py"
cfg["scorer"] = f"{len(inputs):03d}-{scorer.name}"
inputs.append(scorer)
config = root / "configs/focus-primary/replacement-comparison-v1.json"
with config.open("x") as f:
    json.dump(cfg, f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Compare fixed newprimary replacements individually and jointly; no OOF weight tuning",
    "--config",
    str(config),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/compare_focus_primary_replacements.py"], cwd=root, check=True
)
