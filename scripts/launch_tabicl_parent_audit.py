"""Collect exactly five completed parent-CV TabICL runs."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(p):
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


runs = re.findall(
    r"run_id=(\S+)", (root / "artifacts/launches/tabicl-parent-20260911.log").read_text()
)
runs = [
    r
    for r in runs
    if (root / f"artifacts/runs/regression/{r}/outputs/tabicl-parent-oof.parquet").is_file()
]
assert len(runs) == 5
cfg = {
    "oof": [
        add(root / f"artifacts/runs/regression/{r}/outputs/tabicl-parent-oof.parquet") for r in runs
    ],
    "baseline": add(
        root / "artifacts/runs/regression/20260911T110118749415Z-5167f1a0/outputs/oof.parquet"
    ),
    "scorer": add(root / "artifacts/metric-reference/custom_scoring_functions.py"),
    "database": json.loads((root / "configs/parent-cv/audit.json").read_text())["database"],
}
p = root / "configs/parent-cv/tabicl-audit.json"
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
    "Compare TabICL compact/full and LightGBM on corrected parent-retained CV",
    "--config",
    str(p),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/audit_tabicl_parent.py"], check=True, cwd=root)
