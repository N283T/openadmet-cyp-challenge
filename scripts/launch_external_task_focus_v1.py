"""Run intrinsic external-task comparisons without downstream challenge readouts."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "artifacts/runs/shared/20260916T082302554619Z-656c63e6/outputs"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
CONFIGS = [
    ROOT / "configs/external-training/task-focus-broad-continuous-v1.json",
    ROOT / "configs/external-training/task-focus-cyp-continuous-v1.json",
    ROOT / "configs/external-training/task-focus-cyp-continuous-ranking-v1.json",
]
MANIFEST = ROOT / "artifacts/launches/external_task_focus_v1_manifest.json"


def run(config: Path) -> str:
    purpose = json.loads(config.read_text())["purpose"]
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        purpose,
        "--config",
        str(config),
        "--input",
        f"external-compounds.parquet={CORPUS / 'external-compounds.parquet'}",
        "--input",
        f"observations.parquet={CORPUS / 'observations-endpoint-only.parquet'}",
        "--input",
        f"compounds.parquet={CHALLENGE / 'compounds.parquet'}",
        "--",
        "python",
        "scripts/train_external_medium.py",
    ]
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    print(result.stdout, flush=True)
    if result.stderr:
        print(result.stderr, flush=True)
    result.check_returncode()
    ids = re.findall(r"run_id=(\S+)", result.stdout)
    if len(ids) != 1:
        raise ValueError(ids)
    return ids[0]


manifest = []
for config in CONFIGS:
    run_id = run(config)
    manifest.append({"config": str(config.relative_to(ROOT)), "run_id": run_id})
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")

print(json.dumps(manifest, indent=2), flush=True)
