"""Snapshot saved inner selection artifacts and launch two fixed calibration variants."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source_run = "20260912T054406239560Z-96fce32b"
source = root / f"artifacts/runs/regression/{source_run}"
files = []


def add(path):
    assert path.is_file()
    files.append(path)
    return f"{len(files) - 1:03d}-{path.name}"


cfg = dict(
    source_run=source_run,
    scale_bounds=[0.5, 1.5],
    offset_bounds=[-0.5, 0.5],
    variants=["offset", "affine"],
    optimizer="scipy linprog highs",
    database=dict(host=str(root / ".runtime/socket"), port=5434, dbname="cyp_challenge"),
)
cfg["outer_oof"] = add(source / "outputs/external-oof.parquet")
cfg["jeremy_oof"] = add(
    root / "artifacts/runs/regression/20260912T060059705105Z-2f0fbfbd/outputs/external-oof.parquet"
)
cfg["scorer"] = add(root / "artifacts/metric-reference/custom_scoring_functions.py")
cfg["original_fit_source"] = add(source / "source/scripts/nested_external_ensemble.py")
cfg["groups"] = []
for fold in range(5):
    for directory in sorted((source / f"outputs/fold-{fold}").iterdir()):
        if not directory.is_dir():
            continue
        cfg["groups"].append(
            dict(
                inner=add(directory / "inner-predictions.parquet"),
                membership=add(directory / "membership.parquet"),
                selected=add(directory / "selected.json"),
                history=add(directory / "caruana-history.json"),
            )
        )
assert len(cfg["groups"]) == 20
configs = root / "configs/external-calibration"
configs.mkdir(exist_ok=True)
config = configs / "nested24_offset_affine_v1.json"
with config.open("x") as stream:
    json.dump(cfg, stream, indent=2)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Inner-only nestedCaruana24 offset/affine calibration and fixed Jeremy blend audit",
    "--config",
    str(config),
]
for path in files:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/calibrate_nested_ensemble_run.py"], cwd=root, check=True
)
