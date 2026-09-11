"""Snapshot all three-component OOF runs and current controls."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = re.findall(
    r"run_id=(\S+)", (root / "artifacts/launches/ensemble-components-v2.log").read_text()
)
assert len(runs) == 5
paths = []


def add(p):
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


cfg = {
    "oofs": [
        add(root / "artifacts/runs/regression" / r / "outputs/component-oof.parquet") for r in runs
    ],
    "runs": runs,
    "controls": {},
    "database": {"host": str(root / ".runtime/socket"), "port": 5434, "dbname": "cyp_challenge"},
}
for backend in ["tabicl", "tabpfn"]:
    cfg["controls"]["current_full_" + backend] = add(
        root
        / f"artifacts/runs/regression/20260911T131826604125Z-ebb033cc/outputs/{backend}-oof.parquet"
    )
cfg["scorer"] = add(root / "artifacts/metric-reference/custom_scoring_functions.py")
p = root / "configs/parent-cv/ensemble-components-audit.json"
with p.open("x") as stream:
    json.dump(cfg, stream, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Compare matched component OOF prediction and residual correlations and equal-weight ensembles",
    "--config",
    str(p),
]
for f in paths:
    cmd += ["--input", str(f)]
subprocess.run(cmd + ["--", "python", "scripts/audit_ensemble_components.py"], check=True, cwd=root)
