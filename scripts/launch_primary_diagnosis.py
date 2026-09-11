"""Declare all nested auxiliary predictions and primary measurements for diagnosis."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(p):
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


cfg = {
    "base": add(
        root / "artifacts/runs/regression/20260911T110118749415Z-5167f1a0/outputs/oof.parquet"
    ),
    "primary": add(
        root / "artifacts/runs/shared/20260910T100337258849Z-927eac1c/outputs/targets.parquet"
    ),
    "aux": [],
}
for r in re.findall(
    r"run_id=(\S+)", (root / "artifacts/launches/parent-cv-aux-20260911.log").read_text()
):
    out = root / f"artifacts/runs/shared/{r}/outputs"
    meta = json.loads((out / "summary.json").read_text())
    cfg["aux"].append(
        {
            "excluded_folds": meta["config"]["excluded_folds"],
            "predictions": [add(p) for p in sorted(out.glob("predictions-*.parquet"))],
        }
    )
p = root / "configs/parent-cv/primary-diagnosis.json"
with p.open("x") as f:
    json.dump(cfg, f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Diagnose nested log2fc prediction accuracy and its association with pIC50 residuals",
    "--config",
    str(p),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/diagnose_predicted_primary.py"], check=True, cwd=root
)
