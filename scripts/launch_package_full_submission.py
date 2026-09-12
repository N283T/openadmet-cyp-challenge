"""Package the completed full-data regression ensemble in a recorded run."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = []


def add(path):
    path = root / path
    assert path.is_file(), path
    name = f"{len(inputs):03d}-{path.name}"
    inputs.append(path)
    return name


base = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
cfg = dict(
    measurements=add(f"{base}/measurements.parquet"),
    aliases=add(f"{base}/aliases.parquet"),
    reference=add("data/raw/cyp-challenge-TEST-BLINDED.csv"),
    components=[],
)
for name, run in [
    ("caruana", "20260912T074545177210Z-0d561e17"),
    ("medium_tabpfn", "20260912T074409483212Z-f3024998"),
    ("scores_logfc_tabpfn", "20260912T074557646349Z-da857069"),
    ("ft_logfc_tuned_lgbm", "20260912T074624192288Z-bd58b762"),
    ("cheme_logfc_tabicl", "20260912T074558733332Z-98ab8ceb"),
]:
    path = f"artifacts/runs/regression/{run}"
    status = json.loads((root / path / "result.json").read_text())
    assert status["status"] == "completed"
    cfg["components"].append(
        dict(
            name=name,
            run=run,
            weight=0.2,
            predictions=add(f"{path}/outputs/test-predictions.parquet"),
            membership=add(f"{path}/outputs/membership.parquet"),
            summary=add(f"{path}/outputs/summary.json"),
            result=add(f"{path}/result.json"),
        )
    )
config = root / "configs/submission-full/final-five-mean-v1.json"
with config.open("x") as f:
    json.dump(cfg, f, indent=2)
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Verify full-data membership, 750 identities and five equal components; freeze submission",
    "--config",
    str(config),
]
for path in inputs:
    cmd += ["--input", str(path)]
subprocess.run(
    cmd + ["--", "python", "scripts/package_full_regression_submission.py"], cwd=root, check=True
)
