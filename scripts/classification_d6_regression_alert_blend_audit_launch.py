"""Launch the recorded D6 regression/alert complement audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
historical = root / "artifacts/runs/classification/20260915T113718060068Z-b3a7120e"
vaishnavi = root / "artifacts/runs/classification/20260918T055037631826Z-71dcd147/outputs"
regression = root / "artifacts/runs/classification/20260918T055537966493Z-02baa35e/outputs"
deployed = root / "artifacts/runs/classification/20260915T125939829893Z-1d921ffc/outputs"
inputs = {
    "reference-oof.parquet": historical / "inputs/000-ensemble-oof.parquet",
    "lgbm-oof.parquet": historical / "inputs/001-000-oof.parquet",
    "pfn-oof.parquet": historical / "inputs/002-003-oof-probabilities.parquet",
    "graph-oof.parquet": historical / "inputs/003-006-classification-oof.parquet",
    "fixed-oof.parquet": historical / "inputs/004-oof.parquet",
    "raw600-oof.parquet": historical / "inputs/005-oof-probabilities.parquet",
    "selections.parquet": historical / "inputs/006-selections.parquet",
    "vaishnavi-oof.parquet": vaishnavi / "oof.parquet",
    "vaishnavi-test.parquet": vaishnavi / "test.parquet",
    "regression-oof.parquet": regression / "oof.parquet",
    "regression-test.parquet": regression / "test.parquet",
    "deployed-test-scores.parquet": deployed / "test-scores.parquet",
}
for path in inputs.values():
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
    "Audit fixed rank additions from alert and regression feature D6 models",
    "--config",
    str(root / "configs/classification/d6-regression-alert-blend-audit-v1.json"),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command
    + ["--", "python", "scripts/classification_d6_regression_alert_blend_audit.py"],
    cwd=root,
    check=True,
)
