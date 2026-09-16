"""Launch full SMARTCyp atom and aggregate feature extraction."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
jar = (
    Path.home()
    / "ghq/github.com/MD-Studio/MDStudio_SMARTCyp/mdstudio_smartcyp/bin/smartcyp-2.4.2.jar"
)
compounds = (
    root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"
)
config = root / "configs/classification-round/smartcyp-features-v1.json"
for path in (jar, compounds, config):
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
        "Generate pinned SMARTCyp 3A4 and 2D6 atom-level SoM features for every challenge molecule",
        "--config",
        str(config),
        "--input",
        str(jar),
        "--input",
        str(compounds),
        "--",
        "python",
        "scripts/smartcyp_features.py",
    ],
    cwd=root,
    check=True,
)
