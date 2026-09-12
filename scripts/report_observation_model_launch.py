"""Launch recorded matched observation-MLP arms with immutable inputs."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
b = root / "artifacts/runs"
cfg = root / "configs/report_observation_model.json"
cfg.write_text(
    json.dumps(
        dict(
            arms=["independent", "hill_fixed", "hill_learned"],
            max_epochs=150,
            patience=20,
            split_run="20260911T104703655337Z-a3923255",
            architecture=[647, 256, 128, 16],
            seed=42,
            loss_weights=dict(direct=2, tdi=0.5, emax=0.25, logfc=0.5, hill_prior=0.01),
        ),
        indent=2,
    )
)
inputs = [
    b / "shared/20260912T093252891365Z-271e3020/outputs/observation-targets.parquet",
    b / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    b / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    b / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
    b / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    b / "shared/20260912T082554377951Z-d272774c/outputs/nested-log2fc.parquet",
    root / "artifacts/metric-reference/custom_scoring_functions.py",
]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Independent vs Hill observation MLP; held-out heads; inner epoch selection and refit",
    "--config",
    str(cfg),
]
for p in inputs:
    cmd.extend(["--input", str(p)])
subprocess.run(cmd + ["--", "python", "scripts/report_observation_model.py"], cwd=root, check=True)
