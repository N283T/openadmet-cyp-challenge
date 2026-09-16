"""Launch pinned BioTransformer CYP metabolite feature extraction."""

import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("mode", choices=["smoke", "full"])
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
archive = (
    Path.home()
    / "ghq/bitbucket.org/wishartlab/biotransformer3.0jar/biotransformer3.0jar-6432cf8.tar.gz"
)
compounds = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"
config = root / f"configs/classification-round/biotransformer-{args.mode}-v1.json"
for path in [archive, compounds, config]:
    if not path.is_file():
        raise FileNotFoundError(path)
subprocess.run(
    [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        f"{args.mode.capitalize()} pinned BioTransformer CYP product-derived features",
        "--config",
        str(config),
        "--input",
        str(archive),
        "--input",
        str(compounds),
        "--",
        "python",
        "scripts/biotransformer_features.py",
    ],
    cwd=root,
    check=True,
)
