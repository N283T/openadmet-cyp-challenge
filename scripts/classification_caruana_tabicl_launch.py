"""Extend the audited 18-candidate pool with matched TabICL comparisons."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
prior_run = runs / "20260912T133027896889Z-8e3be966"
frozen = runs / "20260912T134440187114Z-7884042f/outputs"
jeremy = runs / "20260912T140401138774Z-0e4bcc8d/outputs"
stack = runs / "20260912T140805013061Z-181ca91b/outputs"
prior = json.loads((prior_run / "config.json").read_text())
paths = []


def add(path):
    assert path.is_file(), path
    paths.append(path)
    return f"{len(paths) - 1:03d}-{path.name}"


sources = [
    dict(
        prefix=source["prefix"],
        **{
            kind: add(prior_run / "inputs" / source[kind])
            for kind in ["oof", "inner", "test"]
        },
    )
    for source in prior["sources"]
]
for prefix, directory, names in [
    (
        "icl_frozen",
        frozen,
        ("outer-oof.parquet", "inner-probabilities.parquet", "fold-test-probabilities.parquet"),
    ),
    (
        "jeremy8",
        jeremy,
        ("outer-oof.parquet", "inner-probabilities.parquet", "fold-test-probabilities.parquet"),
    ),
    (
        "icl_stack",
        stack,
        ("oof.parquet", "inner-probabilities.parquet", "fold-test-probabilities.parquet"),
    ),
]:
    sources.append(
        dict(
            prefix=prefix,
            oof=add(directory / names[0]),
            inner=add(directory / names[1]),
            test=add(directory / names[2]),
        )
    )

cfg = dict(
    sources=sources,
    steps=50,
    bagging_seed=731,
    bagging_samples=20,
    parent_run=prior_run.name,
    tabicl_runs=[frozen.parent.name, jeremy.parent.name, stack.parent.name],
    recipes=prior["recipes"],
)
assert cfg["recipes"] == [
    "fixed_historical",
    "subset1",
    "subset2",
    "subset3",
    "caruana_mcc",
    "caruana_mcc_sparse3",
    "caruana_auc",
    "caruana_brier",
    "bagged_sparse3",
]
config = root / "configs/classification-round/caruana-tabicl-v3.json"
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
    "Compare 24 classifier candidates after adding TabICL and matched Jeremy readouts",
    "--config",
    str(config),
]
for path in paths:
    cmd += ["--input", str(path)]
subprocess.run(cmd + ["--", "python", "scripts/classification_caruana.py"], check=True, cwd=root)
