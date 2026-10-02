"""Record a frozen-score composition check for the complete legacy refit recipe."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
ids = dict(
    lgbm="20260924T035725764059Z-9b7d6d23",
    pfn_pca="20260924T034908336839Z-ec9287a1",
    pfn_raw="20260924T035915695580Z-7dd31f09",
    graph="20260924T040130129624Z-558f74b9",
    fixed="20260924T040107573312Z-887ed73b",
    control="20260924T040301085704Z-914bb06c",
)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Verify the nested legacy comparator exactly recomposes all eleven historical atoms",
    "--config",
    str(root / "configs/classification/d6-legacy-composition-audit-v1.json"),
]
for alias, run in ids.items():
    path = runs / run / "outputs/all-pool-oof.parquet"
    if not path.is_file():
        raise FileNotFoundError(path)
    command += ["--input", f"{alias}.parquet={path}"]
subprocess.run(
    command + ["--", "python", "scripts/audit_cyp_d6_legacy_composition.py"], cwd=root, check=True
)
