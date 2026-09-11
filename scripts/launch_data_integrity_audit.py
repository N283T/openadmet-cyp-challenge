"""Declare all audit inputs and launch a durable audit run."""

import subprocess
from pathlib import Path

base = Path("artifacts/runs")
inputs = sorted(Path("data/raw").glob("*.csv"))
inputs += [
    base / "shared/20260909T132734119233Z-a7063ccd/outputs" / n
    for n in [
        "aliases.parquet",
        "compounds.parquet",
        "standardization.parquet",
        "measurements.parquet",
    ]
]
inputs += [
    base / "regression/20260911T131826604125Z-ebb033cc/outputs" / n
    for n in ["tabicl-oof.parquet", "tabicl-regression.parquet", "tabicl-test-folds.parquet"]
]
inputs += [base / "regression/20260911T132657697271Z-bb9311b7/outputs/regression.parquet"]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Audit official data freshness, SMILES identities, labels and submission alignment",
    "--config",
    "configs/regression/data-integrity-audit.json",
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/audit_data_integrity.py"], check=True)
