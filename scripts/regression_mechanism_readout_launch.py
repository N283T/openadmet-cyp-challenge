"""Launch paired mechanism-feature additions under retained regression parent CV."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
smartcyp = runs / "classification/20260914T044537363003Z-5e7f9e02/outputs"
attnsom = runs / "classification/20260914T045809423993Z-c5c2c254/outputs"
cypmap = runs / "classification/20260914T062717289785Z-17d2a645/outputs"
config = root / "configs/regression/mechanism-representation-readout-v1.json"
inputs = [
    runs / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    runs / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    runs
    / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    runs / "shared/20260912T082554377951Z-d272774c/outputs/nested-log2fc.parquet",
    smartcyp / "smartcyp-features.parquet",
    attnsom / "attnsom-features.parquet",
    attnsom / "attnsom-embeddings.parquet",
    cypmap / "cypmap-features.parquet",
    cypmap / "cypmap-embeddings.parquet",
    runs / "regression/20260912T083010595992Z-1e7001de/outputs/external-oof.parquet",
]
for path in [config, *inputs]:
    if not path.is_file():
        raise FileNotFoundError(path)

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Paired parent-CV additions of SMARTCyp, ATTNSOM, and CYP-MAP to frozen600 plus nested log2fc",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/regression_mechanism_readout.py"],
    cwd=root,
    check=True,
)
