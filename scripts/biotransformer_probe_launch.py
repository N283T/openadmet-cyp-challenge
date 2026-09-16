"""Launch a recorded BioTransformer 3.0 CYP smoke test."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
archive = (
    Path.home()
    / "ghq/bitbucket.org/wishartlab/biotransformer3.0jar/biotransformer3.0jar-6432cf8.tar.gz"
)
config = root / "configs/classification-round/biotransformer-probe-v1.json"
for path in [archive, config]:
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
        "Verify pinned local BioTransformer CYP metabolite generation and output schema",
        "--config",
        str(config),
        "--input",
        str(archive),
        "--",
        "python",
        "scripts/biotransformer_probe.py",
    ],
    cwd=root,
    check=True,
)
