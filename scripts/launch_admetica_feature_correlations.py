"""Declare every feature source for the train/test correlation audit."""

import subprocess
from pathlib import Path

base = Path("artifacts/runs")
files = [
    base / "regression/20260911T145352053201Z-cb9e5d6a/outputs/direct-predictions.parquet",
    base / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
]
files += [
    base / "shared/20260911T125618165690Z-4e52f1df/outputs" / n
    for n in ["oof-log2fc.parquet", "test-log2fc.parquet"]
]
files += sorted((base / "shared/20260909T132918829974Z-bc81af58/outputs").glob("mordred-*.parquet"))
files += sorted(
    (base / "shared/20260910T012630927160Z-e9576045/outputs").glob("chemeleon-*.parquet")
)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Compare Admetica features with log2fc, Mordred and embeddings in train and test",
]
for p in files:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/admetica_feature_correlations.py"], check=True)
