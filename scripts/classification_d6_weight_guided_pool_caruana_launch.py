"""Launch weight-guided pool adjustment with unrestricted Caruana."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/classification/d6-weight-guided-pool-caruana-v1.json"
inputs = {
    "legacy-bank.parquet": runs
    / "20260915T121009912691Z-73198c15/outputs/prediction-bank.parquet",
    "lgbm-oof.parquet": runs
    / "20260919T155849130417Z-0394b600/outputs/oof.parquet",
    "readout-oof.parquet": runs
    / "20260923T161449619218Z-a5f9f7ea/outputs/oof.parquet",
    "tabpfn-oof.parquet": runs
    / "20260923T161718040624Z-8a9676e7/outputs/oof.parquet",
    "pool-definition.json": root
    / "configs/classification/d6-flat-common-pool-caruana-v1.json",
}
declared = json.loads(config.read_text())
needed = {
    declared[key]
    for key in [
        "legacy_bank",
        "lgbm_oof",
        "readout_oof",
        "tabpfn_oof",
        "pool_definition",
    ]
}
if needed != set(inputs):
    raise ValueError(f"Input mismatch: missing={needed-set(inputs)}, extra={set(inputs)-needed}")
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
    "Adjust the CYP2D6 pool from unrestricted Caruana weights and inner-fold stability",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command
    + ["--", "python", "scripts/classification_d6_weight_guided_pool_caruana.py"],
    cwd=root,
    check=True,
)
