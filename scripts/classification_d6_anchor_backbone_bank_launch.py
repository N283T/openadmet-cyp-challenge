"""Launch the new-family-only D6 anchor/CheMeleon flat bank."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/classification/d6-anchor-backbone-bank-v1.json"
inputs = {
    "lgbm-oof.parquet": runs
    / "20260919T155849130417Z-0394b600/outputs/oof.parquet",
    "lgbm-test.parquet": runs
    / "20260919T155849130417Z-0394b600/outputs/test.parquet",
    "readout-oof.parquet": runs
    / "20260923T160616365102Z-9f963692/outputs/oof.parquet",
    "readout-test.parquet": runs
    / "20260923T160616365102Z-9f963692/outputs/test.parquet",
    "tabpfn-oof.parquet": runs
    / "20260923T160851829443Z-eb87a9de/outputs/oof.parquet",
    "tabpfn-test.parquet": runs
    / "20260923T160851829443Z-eb87a9de/outputs/test.parquet",
}
declared = json.loads(config.read_text())
needed = {
    declared[key]
    for key in [
        "lgbm_oof",
        "lgbm_test",
        "readout_oof",
        "readout_test",
        "tabpfn_oof",
        "tabpfn_test",
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
    (
        "Build and outer-cross-fit a new-family-only CYP2D6 flat bank around "
        "anchored Jeremy and CheMeleon"
    ),
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_anchor_backbone_bank.py"],
    cwd=root,
    check=True,
)
