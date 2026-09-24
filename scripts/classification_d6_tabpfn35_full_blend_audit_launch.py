"""Launch fixed blends of full TabPFN 3.5 into the deployed D6 score."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
historical = root / "artifacts/runs/classification/20260915T113718060068Z-b3a7120e"
tabpfn35 = root / "artifacts/runs/classification/20260918T063852992056Z-387d4e0e/outputs"
deployed = root / "artifacts/runs/classification/20260915T125939829893Z-1d921ffc/outputs"
inputs = {
    "reference-oof.parquet": historical / "inputs/000-ensemble-oof.parquet",
    "lgbm-oof.parquet": historical / "inputs/001-000-oof.parquet",
    "pfn-oof.parquet": historical / "inputs/002-003-oof-probabilities.parquet",
    "graph-oof.parquet": historical / "inputs/003-006-classification-oof.parquet",
    "fixed-oof.parquet": historical / "inputs/004-oof.parquet",
    "raw600-oof.parquet": historical / "inputs/005-oof-probabilities.parquet",
    "selections.parquet": historical / "inputs/006-selections.parquet",
    "tabpfn35-oof.parquet": tabpfn35 / "oof.parquet",
    "tabpfn35-test.parquet": tabpfn35 / "test.parquet",
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
    "Audit fixed full TabPFN 3.5 rank additions to the deployed CYP2D6 score",
    "--config",
    str(root / "configs/classification/d6-tabpfn35-full-blend-audit-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command
    + ["--", "python", "scripts/classification_d6_regression_alert_blend_audit.py"],
    cwd=root,
    check=True,
)
