"""Launch target-prior-weighted LightGBM classification models."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(path):
    candidate = root / path
    assert candidate.is_file(), candidate
    paths.append(candidate)
    return f"{len(paths) - 1:03d}-{candidate.name}"


base = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
cfg = {
    name: add(f"{base}/{filename}.parquet")
    for name, filename in [
        ("compounds", "compounds"),
        ("labels", "measurements"),
        ("test_mapping", "test_mapping"),
    ]
}
cfg["assignments"] = add(
    "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet"
)
cfg["class_emb"] = add(
    "artifacts/runs/shared/20260912T051543210441Z-a3ea2ec3/outputs/"
    "external-embeddings.parquet"
)
cfg["mixed_emb"] = add(
    "artifacts/runs/shared/20260912T052853359792Z-2c0a0953/outputs/"
    "external-embeddings.parquet"
)
cfg["mordred"] = [
    add(path.relative_to(root))
    for path in sorted(
        (root / "artifacts/runs/shared/20260909T132918829974Z-bc81af58/outputs").glob(
            "mordred-*.parquet"
        )
    )
]
cfg["lb_summary"] = add(
    "artifacts/runs/classification/20260913T101151963949Z-892b1e27/outputs/summary.json"
)
cfg["split_run"] = "20260909T135130268493Z-edba6fc0"
cfg["endpoints"] = ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]
config = root / "configs/classification-round/lb-weighted-lgbm-v1.json"
config.write_text(json.dumps(cfg, indent=2) + "\n")
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Test six LightGBM views weighted for public-LB inferred hidden class priors",
    "--config",
    str(config),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_lb_weighted_lgbm.py"],
    cwd=root,
    check=True,
)
