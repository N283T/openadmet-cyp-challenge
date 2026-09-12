"""Declare all saved classifier probabilities for comprehensive diagnostics."""

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


run = root / "artifacts/runs/classification/20260912T121146867140Z-b0a05986"
old = json.loads((run / "config.json").read_text())
sources = []
for src in old["sources"]:
    sources.append(
        dict(
            prefix=src["prefix"],
            **{k: add(run / "inputs" / src[k]) for k in ["oof", "inner", "test"]},
        )
    )
sources.append(
    dict(
        prefix="mean",
        oof=add(run / "outputs/ensemble-oof.parquet"),
        test=add(run / "outputs/ensemble-test.parquet"),
    )
)
sources.append(
    dict(
        prefix="old",
        method="baseline",
        oof=add("artifacts/runs/shared/20260910T123856868493Z-243a1e23/outputs/oof.parquet"),
    )
)
cfg = dict(
    sources=sources,
    canonical=add(
        "artifacts/runs/classification/20260912T115018743045Z-269e598d/outputs/classification-labels.parquet"
    ),
    seed=731,
    bootstrap=500,
    permutations=1000,
    focus=["mean::all_plain_pfn_graph_mean", "mean::lgbm_balanced_mean", "old::baseline"],
)
config = root / "configs/classification-round/prediction-diagnostics-v1.json"
assert not config.exists()
config.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Audit ranking, folds, calibration, oracle cutoffs and probability shift for25 classifiers",
    "--config",
    str(config),
]
for p in paths:
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/tdi_prediction_diagnostics.py"], check=True, cwd=root
)
