"""Declare fold-specific encoder lineage and launch bounded inner-only readout HPO."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(path):
    p = root / path
    paths.append(p)
    return f"{len(paths) - 1:03d}-{p.name}"


ft = "artifacts/runs/regression/20260912T060538153709Z-8d434eee/outputs"
cfg = dict(
    method="jeremy_ft_logfc_rustuna_v1",
    seed=42,
    n_trials=24,
    inner_fraction=0.2,
    threads=3,
    workers=4,
    split_run="20260911T104703655337Z-a3923255",
    embedding_run="20260912T060538153709Z-8d434eee",
    embeddings=[add(f"{ft}/fold{fold}/external-embeddings.parquet") for fold in range(5)],
    encoder_memberships=[add(f"{ft}/fold{fold}/membership.parquet") for fold in range(5)],
    encoder_summary=add(f"{ft}/summary.json"),
    selection="24 readout trials per endpoint and outer fold; only inner-group official ST-RAE",
    inner_hpo_optimism="The FT encoder used inner pIC50 labels. Inner HPO scores are optimistic. "
    "All nonanchor outer-fold compounds were excluded from encoder labels and readout HPO.",
    upstream_pretrained_diagnostic=True,
)
for p in [
    "artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs/nested-log2fc.parquet",
    "artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs/summary.json",
    "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    "artifacts/metric-reference/custom_scoring_functions.py",
]:
    add(p)
config = root / "configs/external-ft-tuning/jeremy_ft_logfc_rustuna_v1.json"
config.parent.mkdir(exist_ok=True)
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
    "Inner-only Rustuna LGBM tuning of outer-safe Jeremy FT600 and log2fc4",
    "--config",
    str(config),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/tune_finetuned_readout.py"], cwd=root, check=True)
