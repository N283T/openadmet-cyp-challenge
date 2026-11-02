"""Launch the nested reselection over the incumbent and all frozen-embedding members."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
regression = root / "artifacts/runs/regression"
screen = regression / "20261029T145728905917Z-9d7ea648"
inputs: dict[str, Path] = {
    f"control-{name}.parquet": screen / f"inputs/control-{name}.parquet"
    for name in ("1a2", "2c9", "2d6", "3a4")
}
for alias, run in {
    "minimol-oof.parquet": "20261029T044856735303Z-16745644",
    "monroe-oof.parquet": "20261029T045719539314Z-43b3db57",
    "monroe-ext-oof.parquet": "20261029T115010020552Z-7bae1535",
    "orbmol-oof.parquet": "20261029T124544572615Z-c6b401f4",
    "cypmol1-oof.parquet": "20261029T131828203894Z-956823ad",
    "cypmol10-oof.parquet": "20261029T145728905917Z-9d7ea648",
}.items():
    inputs[alias] = regression / run / "outputs/predictions/oof.parquet"
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    (
        "Nested reselection over the incumbent and all six frozen-embedding member screens: "
        "incumbent, best single column, anchored and free greedy selection"
    ),
    "--config",
    str(root / "configs/regression/bank-selection-v1.json"),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/regression_bank_selection.py"], cwd=root, check=True
)
