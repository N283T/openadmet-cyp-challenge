"""Record fold-excluded selection from fixed robust-bound D6 candidates."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
path = root / "artifacts/runs/regression/20260926T085959697226Z-0dbf57d4/outputs/oof.parquet"
if not path.is_file():
    raise FileNotFoundError(path)
subprocess.run(
    [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "Estimate selection optimism with held-fold-excluded robust-bound recipe choice",
        "--config",
        str(root / "configs/regression/d6-robust-lofo-v1.json"),
        "--input",
        f"blend-oof.parquet={path}",
        "--",
        "python",
        "scripts/regression_d6_robust_lofo.py",
    ],
    cwd=root,
    check=True,
)
