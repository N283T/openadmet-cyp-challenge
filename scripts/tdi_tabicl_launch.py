"""Record all inputs for the sixty matched frozen TabICL classifier fits."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
old = root / "artifacts/runs/classification/20260912T115242271820Z-51023dac"
cfg = json.loads((root / "configs/classification-tabpfn/frozen-v1.json").read_text())
paths = []


def add(path):
    assert path.is_file(), path
    name = f"{len(paths):03d}-{path.name}"
    paths.append(path)
    return name


for key in ["labels", "assignments", "test_mapping"]:
    cfg[key] = add(old / "inputs" / cfg[key])
cfg["weights"] = add(
    root
    / "artifacts/runs/shared/20260912T134236392057Z-8c08c1de/outputs"
    / "tabicl-classifier-v2-20260212.ckpt"
)
for view in cfg["views"]:
    view["files"] = [add(old / "inputs" / f) for f in view["files"]]
    view["name"] = "tabicl_" + view["name"]
cfg["model"] = {
    "package": "tabicl2.2.0",
    "n_estimators": 8,
    "batch_size": 8,
    "random_state": "42+outer",
    "device": "cuda",
    "allow_auto_download": False,
}
config = root / "configs/classification-tabicl/frozen-v1.json"
config.parent.mkdir(exist_ok=True)
config.write_text(json.dumps(cfg, indent=2))
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Sixty TabICL fits matching frozen TabPFN views with honest inner thresholds",
    "--config",
    str(config),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/tdi_tabicl_frozen.py"], cwd=root, check=True)
