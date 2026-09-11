"""Audit exactly five completed matched-feature TabICL folds."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(p):
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


cfg = {
    "oof": [],
    "test": [],
    "database": json.loads((root / "configs/calibration/package-interval-v1.json").read_text())[
        "database"
    ],
}
folds = []
for run in sorted((root / "artifacts/runs/regression").iterdir()):
    if not (run / "outputs/tabicl-oof.parquet").is_file():
        continue
    c = json.loads((run / "config.json").read_text())
    folds.append(c["folds"][0]["fold"])
    cfg["oof"].append(add(run / "outputs/tabicl-oof.parquet"))
    cfg["test"].append(add(run / "outputs/tabicl-test.parquet"))
assert sorted(folds) == list(range(5)), folds
cfg["labels"] = add(
    root
    / "artifacts/runs/regression/20260911T050944039369Z-36d00d9c/outputs"
    / "regression-oof-diagnostics.parquet"
)
cfg["distribution"] = add(
    root
    / "artifacts/runs/regression"
    / "20260911T072717606484Z-8ec7593a/outputs/distribution-comparison-oof.parquet"
)
cfg["scorer"] = add(root / "artifacts/metric-reference/custom_scoring_functions.py")
p = root / "configs/calibration/tabicl-audit-v1.json"
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
    "Audit matched-feature TabICL versus original TabPFN on all five UMAP folds",
    "--config",
    str(p),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/audit_tabicl.py"], check=True, cwd=root)
