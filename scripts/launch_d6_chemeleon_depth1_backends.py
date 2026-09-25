"""Record matched depth-1 TabPFN and TabICL CYP2D6 readouts."""

import argparse
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("backend", choices=["tabpfn", "tabicl"])
parser.add_argument("--log2fc", action="store_true")
args = parser.parse_args()
backend = args.backend
depth = root / "artifacts/runs/shared/20260925T065503727832Z-1c269e23/outputs/depth-1.parquet"
data = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
assignments = (
    root
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
sources = {
    "depth-1.parquet": depth,
    "measurements.parquet": data / "measurements.parquet",
    "parent-assignments.parquet": assignments,
}
if args.log2fc:
    sources["nested-log2fc.parquet"] = (
        root
        / "artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs"
        / "nested-log2fc.parquet"
    )
if backend == "tabpfn":
    sources["tabpfn-regressor.ckpt"] = (
        root
        / "artifacts/runs/regression/20260912T072303592078Z-07f2749a/inputs"
        / "004-tabpfn-v2.6-regressor-v2.6_default.ckpt"
    )
    suffix = "-log2fc" if args.log2fc else ""
    config = root / f"configs/screen/d6-chemeleon-depth1-tabpfn{suffix}-v1.json"
    script = "scripts/external_embedding_readout.py"
else:
    sources["tabicl-regressor.ckpt"] = (
        root
        / "artifacts/runs/shared/20260911T093102177307Z-beddd36f/outputs"
        / "tabicl-regressor-v2-20260212.ckpt"
    )
    suffix = "-log2fc" if args.log2fc else ""
    config = root / f"configs/screen/d6-chemeleon-depth1-tabicl{suffix}-v1.json"
    script = "scripts/external_tabicl_readout.py"
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    f"CYP2D6 parent CV CheMeleon depth-1 {backend} atomic readout",
    "--config",
    str(config),
]
for alias, path in sources.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", script])
subprocess.run(command, cwd=root, check=True)
