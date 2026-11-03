"""Launch the grand-bank nested selection and payload build."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
screen = runs / "regression/20261029T145728905917Z-9d7ea648"
inputs = {
    f"control-{name}.parquet": screen / f"inputs/control-{name}.parquet"
    for name in ("1a2", "2c9", "2d6", "3a4")
}
inputs["base-predictions.csv"] = (
    root / "artifacts/submissions/regression/20261029T160903346294Z-52e12de6/predictions.csv"
)
inputs["test-mapping.parquet"] = (
    runs / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet"
)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    (
        "Nested anchored selection over every deployable atomic OOF column in all regression runs; "
        "adopt endpoints passing the prespecified rule or chosen by the user (1A2) "
        "and build a payload on 52e12de6"
    ),
    "--config",
    str(root / "configs/regression/grand-bank-package-v1.json"),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(command + ["--", "python", "scripts/regression_grand_bank.py"], cwd=root, check=True)
