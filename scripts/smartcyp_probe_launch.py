"""Launch the recorded SMARTCyp compatibility probe."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
jar = (
    Path.home()
    / "ghq/github.com/MD-Studio/MDStudio_SMARTCyp/mdstudio_smartcyp/bin/smartcyp-2.4.2.jar"
)
config = root / "configs/classification-round/smartcyp-probe-v1.json"
if not jar.is_file():
    raise FileNotFoundError(jar)

subprocess.run(
    [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Verify pinned SMARTCyp 2.4.2 locally and preserve its raw atom-level output",
        "--config",
        str(config),
        "--input",
        str(jar),
        "--",
        "python",
        "scripts/smartcyp_probe.py",
    ],
    cwd=root,
    check=True,
)
