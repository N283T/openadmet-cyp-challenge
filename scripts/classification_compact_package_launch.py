"""Record exact prediction sources for compact/restored classification candidates."""

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
    current_test=add(base / "20260912T121146867140Z-b0a05986/outputs/ensemble-test.parquet"),
    current_oof=add(base / "20260912T121146867140Z-b0a05986/outputs/ensemble-oof.parquet"),
    lgbm_test=add(base / "20260912T115220209511Z-76facc71/outputs/fold-test-probabilities.parquet"),
    lgbm_oof=add(base / "20260912T115220209511Z-76facc71/outputs/oof.parquet"),
    old_submission=add(base / "20260910T133042436909Z-c95bbbfc/outputs/classification.parquet"),
    aliases=add(
        root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet"
    ),
    reference=add(root / "data/raw/cyp-challenge-TEST-BLINDED.csv"),
)
config = root / "configs/classification-round/compact-package-v1.json"
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
    "Prepare exact old A4 restoration and compact D6 single-LGBM alternative; no upload",
    "--config",
    str(config),
]
for path in paths:
    cmd += ["--input", str(path)]
subprocess.run(
    cmd + ["--", "python", "scripts/classification_compact_package.py"], cwd=root, check=True
)
