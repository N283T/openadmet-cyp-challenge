"""Declare every input to the fixed context/blend experiment."""

import json
import subprocess
from pathlib import Path

path = Path("configs/calibration/context-blend-v1.json")
cfg = json.loads(path.read_text())
inputs = [Path(e["path"]) for e in cfg["members"]]
for i, e in enumerate(cfg["members"]):
    e["input"] = f"{i:03d}-{inputs[i].name}"
inputs += [
    Path(p)
    for p in [
        "artifacts/runs/regression/20260911T050944039369Z-36d00d9c/outputs/regression-oof-diagnostics.parquet",
        "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
        "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
        "artifacts/metric-reference/custom_scoring_functions.py",
    ]
]
resolved = path.with_name("context-blend-resolved-v1.json")
with resolved.open("x") as f:
    f.write(json.dumps(cfg, indent=2) + "\n")
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Fixed ensembles and context-dependent calibration on disjoint UMAP holdout halves",
    "--config",
    str(resolved),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/context_blend_screen.py"], check=True)
