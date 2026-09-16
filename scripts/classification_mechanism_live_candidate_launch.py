"""Launch the prespecified endpoint-specific live/mechanism classification candidate."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
inputs = [
    root
    / "artifacts/candidates/classification/20260914T020749578109Z-0a0bd4cc/predictions.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
    runs / "20260913T123420289849Z-b75917f6/inputs/000-oof.parquet",
    runs / "20260913T123420289849Z-b75917f6/inputs/001-crossfit-selection.csv",
    runs / "20260913T123420289849Z-b75917f6/inputs/002-deployment.parquet",
    runs / "20260914T051140518895Z-ca43af17/outputs/oof.parquet",
    runs / "20260914T051140518895Z-ca43af17/outputs/test.parquet",
    runs / "20260913T153156161201Z-1954e23e/inputs/000-member-oof.parquet",
    runs / "20260913T153156161201Z-1954e23e/inputs/002-oof.parquet",
    runs / "20260913T153156161201Z-1954e23e/outputs/test-blends.parquet",
    runs / "20260914T063614433212Z-c0dc3f32/outputs/oof.parquet",
    runs / "20260914T063614433212Z-c0dc3f32/outputs/test.parquet",
    root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
]
for path in inputs:
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
    (
        "Audit and package a prespecified 3:1 live/mechanism rank blend: "
        "D6 Morgan-SMARTCyp and A4 Morgan-SMARTCyp-ATTNSOM-CYP-MAP; no upload"
    ),
    "--config",
    str(root / "configs/classification-round/mechanism-live-blend-candidate-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_mechanism_live_candidate.py"],
    cwd=root,
    check=True,
)
