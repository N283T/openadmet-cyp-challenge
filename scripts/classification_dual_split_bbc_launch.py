"""Launch the recorded group-bootstrap selection audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
oof = root / "artifacts/runs/classification/20260913T134843805812Z-0d2682af/outputs/oof.parquet"
assert oof.is_file(), oof
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Estimate compact dual-split candidate-selection optimism with grouped "
        "bootstrap out-of-bag evaluation"
    ),
    "--config",
    str(root / "configs/classification-round/dual-split-bbc-v1.json"),
    "--input",
    str(oof),
]
subprocess.run(
    command + ["--", "python", "scripts/classification_dual_split_bbc.py"],
    cwd=root,
    check=True,
)
