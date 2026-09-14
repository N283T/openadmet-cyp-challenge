"""Launch the recorded FDA/Instem CYP3A4 TDI source audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root
    / (
        "artifacts/runs/classification/20260913T142848479169Z-594c2424/outputs/"
        "source/Supplementary Material/Supplementary_TableS1_Full_Database.sdf"
    ),
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
]
for path in inputs:
    assert path.is_file(), path
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Audit the checksum-validated FDA/Instem CYP3A4 TDI SDF, duplicates, "
        "prevalence, and exact connectivity overlap with challenge train/test"
    ),
    "--config",
    str(root / "configs/external/fda-cyp3a4-tdi-audit-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/fda_tdi_source_audit.py"], cwd=root, check=True)
