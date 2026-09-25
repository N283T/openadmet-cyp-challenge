"""Record CheMeleon depth-1 TDI TabPFN and TabICL readouts."""

import argparse
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
parser = argparse.ArgumentParser()
parser.add_argument("backend", choices=["tabpfn", "tabicl"])
backend = parser.parse_args().backend
prep = runs / "classification/20260912T115018743045Z-269e598d/outputs"
sources = {
    "classification-labels.parquet": prep / "classification-labels.parquet",
    "umap-assignments.parquet": prep / "umap-assignments.parquet",
    "test_mapping.parquet": (
        runs / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet"
    ),
    "depth-1.parquet": (
        runs / "shared/20260925T065503727832Z-1c269e23/outputs/depth-1.parquet"
    ),
}
if backend == "tabpfn":
    sources["tabpfn-classifier.ckpt"] = (
        runs / "shared/20260910T121624900089Z-271684c2/inputs"
        / "004-tabpfn-v2.6-classifier-v2.6_default.ckpt"
    )
    script = "scripts/tdi_tabpfn_frozen.py"
else:
    sources["tabicl-classifier.ckpt"] = (
        runs / "shared/20260912T134236392057Z-8c08c1de/outputs"
        / "tabicl-classifier-v2-20260212.ckpt"
    )
    script = "scripts/tdi_tabicl_frozen.py"
command = [
    "pixi", "run", "cyp", "run", "--track", "classification", "--purpose",
    f"Matched UMAP CheMeleon depth-1 {backend} TDI readout on both endpoints",
    "--config", str(root / f"configs/screen/classification-chemeleon-d1-{backend}-v1.json"),
]
for alias, path in sources.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", script])
subprocess.run(command, cwd=root, check=True)
