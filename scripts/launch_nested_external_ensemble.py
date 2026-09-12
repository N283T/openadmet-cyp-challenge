"""Record all immutable inputs for nested Caruana external feature selection."""

import argparse
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--jeremy-expanded", action="store_true")
parser.add_argument("--classic-expanded", action="store_true")
args = parser.parse_args()
paths = []


def add(path):
    p = root / path
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


cfg = dict(
    seed=42,
    steps=50,
    method="nested_caruana_external_lgbm",
    split_run="20260911T104703655337Z-a3923255",
    assignments=add(
        "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    measurements=add(
        "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet"
    ),
    compounds=add(
        "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"
    ),
    class_emb=add(
        "artifacts/runs/shared/20260912T051543210441Z-a3ea2ec3/outputs/external-embeddings.parquet"
    ),
    mixed_emb=add(
        "artifacts/runs/shared/20260912T052853359792Z-2c0a0953/outputs/external-embeddings.parquet"
    ),
    tree=add(
        "artifacts/runs/shared/20260912T053008533415Z-45e5a1eb/outputs/external-scores.parquet"
    ),
    nested=add(
        "artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs/nested-log2fc.parquet"
    ),
    scorer=add("artifacts/metric-reference/custom_scoring_functions.py"),
    threads=3,
)
config = root / "configs/external-ensemble/nested-caruana-lgbm.json"
if args.jeremy_expanded or args.classic_expanded:
    cfg["jeremy_emb"] = add(
        "artifacts/runs/shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet"
    )
    cfg["method"] = "jeremy_expanded_v1"
    cfg["upstream_pretrained_diagnostic"] = True
    config = root / "configs/external-ensemble/jeremy_expanded_v1.json"
if args.classic_expanded:
    cfg["classic"] = add(
        "artifacts/runs/shared/20260912T061043373484Z-2e344d7a/outputs/fp-mordred-features.parquet"
    )
    cfg["method"] = "jeremy_classic_expanded_v1"
    config = root / "configs/external-ensemble/jeremy_classic_expanded_v1.json"
config.parent.mkdir(parents=True, exist_ok=True)
with config.open("x") as f:
    json.dump(cfg, f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Nested Caruana inner-only selection of external LGBMs"
    + (
        " with Jeremy and train-filtered FP+Mordred candidates"
        if args.classic_expanded
        else " with six frozen Jeremy candidates"
        if args.jeremy_expanded
        else " among 24 candidates"
    ),
    "--config",
    str(config),
]
for p in paths:
    cmd.extend(["--input", str(p)])
subprocess.run(cmd + ["--", "python", "scripts/nested_external_ensemble.py"], check=True, cwd=root)
