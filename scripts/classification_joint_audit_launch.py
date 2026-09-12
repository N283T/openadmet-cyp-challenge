"""Snapshot all completed joint-classifier fit lineage and prediction audit inputs."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runid = "20260912T120135164747Z-1da37202"
run = root / "artifacts/runs/classification" / runid
paths = []


def add(path):
    assert path.is_file(), path
    name = f"{len(paths):03d}-{path.name}"
    paths.append(path)
    return name


cfg = dict(model_run=runid, fits=[])
for key, index in [("compounds", 0), ("labels", 1), ("assignments", 2), ("testmap", 3), ("aux", 4)]:
    cfg[key] = add(next((run / "inputs").glob(f"{index:03d}-*")))
for key, name in [
    ("contract", "target-contract.json"),
    ("oof", "classification-oof.parquet"),
    ("inner", "inner-probabilities.parquet"),
    ("test", "test-probabilities.parquet"),
    ("foldtest", "fold-test-probabilities.parquet"),
]:
    cfg[key] = add(run / "outputs" / name)
for recipe in ["tdi_only", "tdi_aux"]:
    for fold in range(5):
        directory = run / "outputs" / f"{recipe}-fold{fold}"
        item = dict(recipe=recipe, fold=fold)
        for phase in ["inner", "refit"]:
            for key, suffix in [
                ("membership", "membership.parquet"),
                ("scaling", "scaling.json"),
                ("history", "history.json"),
            ]:
                item[f"{phase}_{key}"] = add(directory / f"{phase}-{suffix}")
            item[f"{phase}_checkpoint"] = add(
                directory / ("inner-best.pt" if phase == "inner" else "refit.pt")
            )
        cfg["fits"].append(item)
cfg["audited_source"] = add(run / "source/scripts/classification_joint.py")
config = root / "configs/classification-round/joint-audit-v1.json"
config.write_text(json.dumps(cfg, indent=2))
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Audit TDI exclusions, scalers, epoch lineage and canonical probabilities",
    "--config",
    str(config),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_joint_audit.py"], cwd=root, check=True
)
