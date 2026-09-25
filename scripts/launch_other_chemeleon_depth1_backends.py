"""Record parent-CV CheMeleon depth-1 foundation readouts for other targets."""

import argparse
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("backend", choices=["tabpfn", "tabicl"])
backend = parser.parse_args().backend
sources = {
    "depth-1.parquet": (
        root / "artifacts/runs/shared/20260925T065503727832Z-1c269e23/outputs" / "depth-1.parquet"
    ),
    "measurements.parquet": (
        root
        / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
        / "measurements.parquet"
    ),
    "parent-assignments.parquet": (
        root
        / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
        / "parent-assignments.parquet"
    ),
    "nested-log2fc.parquet": (
        root
        / "artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs"
        / "nested-log2fc.parquet"
    ),
}
if backend == "tabpfn":
    sources["tabpfn-regressor.ckpt"] = (
        root
        / "artifacts/runs/regression/20260912T072303592078Z-07f2749a/inputs"
        / "004-tabpfn-v2.6-regressor-v2.6_default.ckpt"
    )
    script = "scripts/external_embedding_readout.py"
else:
    sources["tabicl-regressor.ckpt"] = (
        root
        / "artifacts/runs/shared/20260911T093102177307Z-beddd36f/outputs"
        / "tabicl-regressor-v2-20260212.ckpt"
    )
    script = "scripts/external_tabicl_readout.py"
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    f"Parent-CV CheMeleon depth-1 {backend} readout for 1A2, 2C9 and 3A4",
    "--config",
    str(root / f"configs/screen/chemeleon-depth1-other-{backend}-log2fc-v1.json"),
]
for alias, path in sources.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", script])
subprocess.run(command, cwd=root, check=True)
