"""Snapshot shared frozen features for the classification screen."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(path):
    p = root / path
    assert p.is_file(), p
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


base = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
cfg = {
    name: add(f"{base}/{file}.parquet")
    for name, file in [
        ("compounds", "compounds"),
        ("labels", "measurements"),
        ("test_mapping", "test_mapping"),
    ]
}
cfg["assignments"] = add(
    "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet"
)
cfg["split_run"] = "20260909T135130268493Z-edba6fc0"
cfg["class_emb"] = add(
    "artifacts/runs/shared/20260912T051543210441Z-a3ea2ec3/outputs/external-embeddings.parquet"
)
cfg["mixed_emb"] = add(
    "artifacts/runs/shared/20260912T052853359792Z-2c0a0953/outputs/external-embeddings.parquet"
)
cfg["mordred"] = [
    add(p.relative_to(root))
    for p in sorted(
        (root / "artifacts/runs/shared/20260909T132918829974Z-bc81af58/outputs").glob(
            "mordred-*.parquet"
        )
    )
]
assert cfg["mordred"]
cfg["endpoints"] = ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]
config = root / "configs/classification-round/lgbm-v2.json"
config.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Compare 12 frozen-feature LGBM TDI recipes with inner thresholds and saved test probabilities",
    "--config",
    str(config),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/classification_lgbm.py"], cwd=root, check=True)
