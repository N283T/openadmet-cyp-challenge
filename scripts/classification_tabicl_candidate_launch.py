"""Record and validate the target-wise TabICL classification candidate."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
base = root / "artifacts/runs/classification"
paths = []


def add(path):
    assert path.is_file(), path
    paths.append(path)
    return f"{len(paths) - 1:03d}-{path.name}"


cfg = dict(
    current=add(base / "20260912T121146867140Z-b0a05986/outputs/ensemble-test.parquet"),
    caruana=add(base / "20260912T141520708830Z-97df4c71/outputs/deployment.parquet"),
    aliases=add(
        root
        / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet"
    ),
    reference=add(root / "data/raw/cyp-challenge-TEST-BLINDED.csv"),
)
config = root / "configs/classification-round/tabicl-candidate-v1.json"
assert not config.exists()
config.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Package retained submitted D6 with audited A4 TabICL Caruana predictions; no upload",
    "--config",
    str(config),
]
for path in paths:
    cmd += ["--input", str(path)]
subprocess.run(
    cmd + ["--", "python", "scripts/classification_tabicl_candidate.py"], cwd=root, check=True
)
