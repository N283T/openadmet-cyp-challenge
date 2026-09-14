"""Launch fixed-Morgan model-selection audit and packaging."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = [
    root / "artifacts/runs/classification/20260913T123045515718Z-0df22d54/outputs/oof.parquet",
    root
    / (
        "artifacts/runs/classification/20260913T123045515718Z-0df22d54/outputs/"
        "crossfit-selection.csv"
    ),
    root
    / "artifacts/runs/classification/20260913T123045515718Z-0df22d54/outputs/deployment.parquet",
    root
    / "artifacts/candidates/classification/20260913T113803083052Z-1322afeb/predictions.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
    root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
]
assert all(path.is_file() for path in paths)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Audit and package fold-selected fixed-Morgan classical D6 model; no upload",
    "--config",
    str(root / "configs/classification-round/d6-fixed-fp-package-v1.json"),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_fixed_fp_package.py"],
    cwd=root,
    check=True,
)
