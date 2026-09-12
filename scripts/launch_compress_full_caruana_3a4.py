"""Record the fixed top-2 full-data CYP3A4 Caruana compression."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = []


def add(relative: str) -> str:
    path = root / relative
    inputs.append(path)
    return f"{len(inputs) - 1:03d}-{path.name}"


full = "artifacts/runs/regression/20260912T074545177210Z-0d561e17/outputs"
endpoint = "CYP3A4_pIC50_direct_inhibition"
selected_path = root / full / endpoint / "selected.json"
selected = json.loads(selected_path.read_text())
ranked = sorted(selected, key=lambda item: (-item["weight"], item["name"]))[:2]
cfg = {
    "endpoint": endpoint,
    "k": 2,
    "expected_candidates": ["class_quantile_0.65", "morgan_quantile_0.65"],
    "predictions": add(f"{full}/test-predictions.parquet"),
    "membership": add(f"{full}/membership.parquet"),
    "selection": add(f"{full}/{endpoint}/selected.json"),
    "compounds": add(
        "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"
    ),
    "test_mapping": add(
        "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet"
    ),
    "class_emb": add(
        "artifacts/runs/shared/20260912T051543210441Z-a3ea2ec3/outputs/external-embeddings.parquet"
    ),
    "models": {},
    "selection_evidence_run": "20260912T152832952478Z-3e4de7c8",
}
for item in ranked:
    cfg["models"][item["name"]] = add(f"{full}/{endpoint}/{item['name']}.txt")

config = root / "configs/regression/caruana-3a4-full-top2-v1.json"
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
    "Compress the full-data CYP3A4 Caruana component from 20 models to the fixed top two",
    "--config",
    str(config),
]
for path in inputs:
    command.extend(["--input", str(path)])
subprocess.run(
    command + ["--", "python", "scripts/compress_full_caruana_3a4.py"],
    check=True,
    cwd=root,
)
