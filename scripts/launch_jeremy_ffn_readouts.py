"""Prespecified readout comparisons for frozen upstream predictor features."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
features = root / "artifacts/runs/shared/20260912T072105522912Z-53eb675c/outputs"
nested = (
    root / "artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs/nested-log2fc.parquet"
)
for filename, prefix, label, backend in [
    ("last-hidden.parquet", "emb_", "hidden", "lgbm"),
    ("external-scores.parquet", "score_", "scores", "lgbm"),
    ("external-scores.parquet", "score_", "scores", "tabpfn"),
    ("last-hidden.parquet", "emb_", "hidden", "tabpfn"),
]:
    command = [
        "pixi",
        "run",
        "python",
        "scripts/launch_external_readout.py",
        "--features",
        str(features / filename),
        "--prefix",
        prefix,
        "--method",
        f"jeremy_ffn_{label}_logfc_{backend}_v1",
        "--backend",
        backend,
        "--nested",
        str(nested),
    ]
    print("START", " ".join(command), flush=True)
    subprocess.run(command, cwd=root, check=True)
