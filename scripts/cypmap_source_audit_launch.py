"""Launch CYP-MAP source-data overlap and provenance audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source = Path.home() / "ghq/github.com/qwon135/CYP-MAP"
shared = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
inputs = [
    shared / "compounds.parquet",
    shared / "test_mapping.parquet",
    source / "data/cyp_map_train.sdf",
    source / "data/cyp_map_test.sdf",
]
config = root / "configs/classification-round/cypmap-source-audit-v1.json"
for path in [config, *inputs]:
    if not path.is_file():
        raise FileNotFoundError(path)

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Audit CYP-MAP source molecules for challenge connectivity overlap",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/cypmap_source_audit.py"], cwd=root, check=True)
