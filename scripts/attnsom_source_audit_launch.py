"""Launch ATTNSOM external-data overlap and provenance audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source = Path.home() / "ghq/github.com/dmis-lab/ATTNSOM"
base = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
inputs = [base / "compounds.parquet", base / "test_mapping.parquet"] + [
    source / "cyp_dataset" / f"{isoform}.sdf"
    for isoform in ["1A2", "2A6", "2B6", "2C8", "2C9", "2C19", "2D6", "2E1", "3A4"]
]
config = root / "configs/classification-round/attnsom-source-audit-v1.json"
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
    "Audit pinned ATTNSOM atom-level CYP data and exclude challenge connectivity leakage",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/attnsom_source_audit.py"],
    cwd=root,
    check=True,
)
