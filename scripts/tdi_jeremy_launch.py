"""Record Jeremy600 matched eight-estimator classifier comparison."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
old = runs / "classification/20260912T115242271820Z-51023dac"
cfg = json.loads((root / "configs/classification-tabpfn/frozen-v1.json").read_text())
paths = []


def add(path):
    assert path.is_file(), path
    name = f"{len(paths):03d}-{path.name}"
    paths.append(path)
    return name


for key in ["labels", "assignments", "test_mapping"]:
    cfg[key] = add(old / "inputs" / cfg[key])
cfg["pfn_weights"] = add(old / "inputs" / cfg["weights"])
cfg["weights"] = add(
    runs / "shared/20260912T134236392057Z-8c08c1de/outputs" / "tabicl-classifier-v2-20260212.ckpt"
)
feature = add(
    runs / "shared/20260912T055819206227Z-c060d197/outputs" / "external-embeddings.parquet"
)
cfg["views"] = [
    dict(
        name=f"{backend}_jeremy_pca256",
        backend=backend,
        source_run="c060d197",
        prefix="emb_",
        files=[feature],
    )
    for backend in ["tabicl", "tabpfn"]
]
cfg["model"] = dict(
    n_estimators=8, seed="42+outer", scaler="StandardScaler", pca_components=256, pca_seed=42
)
config = root / "configs/classification-tabicl/jeremy-v1.json"
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
    "Forty Jeremy600 TabICL8 and TabPFN8 matched outer and honest inner classifiers",
    "--config",
    str(config),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/tdi_jeremy_readout.py"], cwd=root, check=True)
