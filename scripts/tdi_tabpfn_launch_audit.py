"""Record finished TDI probability, label and 60-fit membership audit."""

import argparse
import json
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("inner_run")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
outer = runs / "20260912T115242271820Z-51023dac/outputs"
inner = runs / args.inner_run / "outputs"
paths = []


def add(path):
    assert path.is_file(), path
    name = f"{len(paths):03d}-{path.name}"
    paths.append(path)
    return name


cfg = dict(
    inner_validation=add(inner / "inner-validation-probabilities.parquet"),
    oof=add(outer / "oof-probabilities.parquet"),
    inner=add(inner / "oof-inner-threshold.parquet"),
    test=add(outer / "test-probabilities.parquet"),
    fold_test=add(outer / "test-fold-probabilities.parquet"),
    labels=add(runs / "20260912T115018743045Z-269e598d/outputs/classification-labels.parquet"),
    memberships=[],
)
for kind, directory, marker in [("outer", outer, "fold"), ("inner", inner, "outer")]:
    files = sorted(directory.glob("*-membership.parquet"))
    assert len(files) == 30
    for path in files:
        fold = int(path.name.split(marker)[-1].split("-")[0])
        cfg["memberships"].append(dict(kind=kind, fold=fold, file=add(path)))
config = root / "configs/classification-tabpfn/audit-v1.json"
config.write_text(json.dumps(cfg, indent=2))
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Audit TDI TabPFN probabilities, labels, outer and inner membership; summarize AP and MCC",
    "--config",
    str(config),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/tdi_tabpfn_audit.py"], cwd=root, check=True)
