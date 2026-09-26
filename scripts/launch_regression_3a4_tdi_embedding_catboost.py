"""Record a CatBoost regression peer on the existing fold-specific TDI encoder."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
prior = root / "artifacts/runs/regression/20260912T093415808252Z-708fd656"
inputs = {
    **{
        f"embedding-fold-{fold}.parquet": prior / f"inputs/{fold:03d}-external-embeddings.parquet"
        for fold in range(5)
    },
    "measurements.parquet": prior / "inputs/005-measurements.parquet",
    "parent-assignments.parquet": prior / "inputs/006-parent-assignments.parquet",
    "nested-log2fc.parquet": prior / "inputs/007-nested-log2fc.parquet",
    "tdi-arm-oof.parquet": prior / "outputs/external-oof.parquet",
    "scored-bank-oof.parquet": root
    / "artifacts/runs/regression/20260917T130718808881Z-dfa02730/outputs/oof.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Test CatBoost regression readout of classification-trained 3A4 TDI continuous embeddings",
    "--config",
    str(root / "configs/regression/3a4-tdi-embedding-catboost-v1.json"),
]
for alias, path in inputs.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_3a4_tdi_embedding_catboost.py"],
    cwd=root,
    check=True,
)
