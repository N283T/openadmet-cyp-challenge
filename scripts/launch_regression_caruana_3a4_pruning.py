"""Record and run an endpoint-specific internal Caruana compression experiment."""

import argparse
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--endpoint", choices=["1a2", "3a4"], default="3a4")
args = parser.parse_args()
endpoint = {
    "1a2": "CYP1A2_pIC50_direct_inhibition",
    "3a4": "CYP3A4_pIC50_direct_inhibition",
}[args.endpoint]
inputs = []


def add(relative: str) -> str:
    path = root / relative
    inputs.append(path)
    return f"{len(inputs) - 1:03d}-{path.name}"


original = "artifacts/runs/regression/20260912T054406239560Z-96fce32b"
shared = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
cfg = {
    "endpoint": endpoint,
    "seed": 42,
    "max_models": 20,
    "bootstrap_replicates": 1000,
    "bootstrap_seed": 20260913,
    "deployment_distinct_models": 20,
    "assignments": add(
        "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    "measurements": add(f"{shared}/measurements.parquet"),
    "compounds": add(f"{shared}/compounds.parquet"),
    "class_emb": add(
        "artifacts/runs/shared/20260912T051543210441Z-a3ea2ec3/outputs/external-embeddings.parquet"
    ),
    "mixed_emb": add(
        "artifacts/runs/shared/20260912T052853359792Z-2c0a0953/outputs/external-embeddings.parquet"
    ),
    "tree": add(
        "artifacts/runs/shared/20260912T053008533415Z-45e5a1eb/outputs/external-scores.parquet"
    ),
    "scorer": add("artifacts/metric-reference/custom_scoring_functions.py"),
}
component_runs = [
    ("nested_caruana", "20260912T054406239560Z-96fce32b"),
    ("jeremy_frozen_tabpfn", "20260912T060059705105Z-2f0fbfbd"),
    ("jeremy_scores_log2fc_tabpfn", "20260912T072303592078Z-07f2749a"),
    ("jeremy_finetuned_log2fc_lgbm", "20260912T072056445264Z-e6b6cbcb"),
    ("chemeleon_log2fc_tabicl", "20260912T071619253798Z-73a6fb7e"),
]
cfg["components"] = [
    {
        "name": name,
        "file": add(f"artifacts/runs/regression/{run}/outputs/external-oof.parquet"),
    }
    for name, run in component_runs
]
if args.endpoint == "1a2":
    cfg["ensemble_components"] = [
        "nested_caruana",
        "jeremy_frozen_tabpfn",
        "jeremy_scores_log2fc_tabpfn",
    ]
folds = []
for fold in range(5):
    directory = root / original / "outputs" / f"fold-{fold}" / cfg["endpoint"]
    selection_path = directory / "selected.json"
    selection = json.loads(selection_path.read_text())
    record = {
        "fold": fold,
        "selection": add(str(selection_path.relative_to(root))),
        "nested": add(
            "artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs/nested-log2fc.parquet"
        ),
        "models": {},
    }
    for item in selection["selected"]:
        model = directory / f"{item['name']}.txt"
        record["models"][item["name"]] = add(str(model.relative_to(root)))
    folds.append(record)
cfg["folds"] = folds

version = "v2" if args.endpoint == "1a2" else "v3"
config = root / f"configs/regression/caruana-{args.endpoint}-internal-pruning-{version}.json"
config.parent.mkdir(parents=True, exist_ok=True)
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
    (
        f"Evaluate {endpoint} Caruana internal compression using top-K weighted and "
        "equal averages, including outer-fold simulation of the full-data weight recipe"
    ),
    "--config",
    str(config),
]
for path in inputs:
    command.extend(["--input", str(path)])
subprocess.run(
    command + ["--", "python", "scripts/regression_caruana_3a4_pruning.py"],
    check=True,
    cwd=root,
)
