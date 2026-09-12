"""Snapshot all three-exclusion sources and unchanged historical predictions."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/shared"
registry = json.loads((root / "configs/classification-stack-inner/aux-registry.json").read_text())
assert len(registry) == 10
paths = []


def add(path):
    assert path.is_file(), path
    name = f"{len(paths):03d}-{path.name}"
    paths.append(path)
    return name


cfg = dict(
    labels=add(runs / "20260909T132734119233Z-a7063ccd/outputs/measurements.parquet"),
    assignments=add(runs / "20260909T135130268493Z-edba6fc0/outputs/assignments.parquet"),
    weights=add(
        runs / "20260910T121624900089Z-271684c2/inputs/004-tabpfn-v2.6-classifier-v2.6_default.ckpt"
    ),
    embedding=[
        add(p)
        for p in sorted(
            (runs / "20260910T012630927160Z-e9576045/outputs").glob("chemeleon-*.parquet")
        )
    ],
    mordred=[
        add(p)
        for p in sorted(
            (runs / "20260909T132918829974Z-bc81af58/outputs").glob("mordred-*.parquet")
        )
    ],
    old_outer=add(runs / "20260910T123856868493Z-243a1e23/outputs/oof.parquet"),
    old_test=add(runs / "20260910T132141471698Z-343eb5ec/outputs/fold_predictions.parquet"),
    old_mean=add(runs / "20260910T132141471698Z-343eb5ec/outputs/ensemble_probabilities.parquet"),
    outer_oof_run="20260910T123856868493Z-243a1e23",
    test_run="20260910T132141471698Z-343eb5ec",
    split_id="20260909T135130268493Z-edba6fc0",
    auxiliary=[],
)
assert len(cfg["embedding"]) == 108 and cfg["mordred"]
for source in registry:
    directory = runs / source["run"] / "outputs"
    cfg["auxiliary"].append(
        dict(
            source,
            summary=add(directory / "summary.json"),
            membership=add(directory / "membership.parquet"),
            predictions=[add(p) for p in sorted(directory.glob("predictions-*.parquet"))],
        )
    )
config = root / "configs/classification-stack-inner/readouts-v1.json"
config.write_text(json.dumps(cfg, indent=2))
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Historical TabPFN inner readouts with triple-excluded log2fc features",
    "--config",
    str(config),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_stack_inner.py"], cwd=root, check=True
)
