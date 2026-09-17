"""Train the nine-head CYP continuous model for intrinsic task evaluation."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "artifacts/runs/shared/20260916T082302554619Z-656c63e6/outputs"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
CONFIG = ROOT / "configs/external-training/task-focus-cyp-continuous-shared9-v1.json"
purpose = json.loads(CONFIG.read_text())["purpose"]
command = [
    "pixi", "run", "cyp", "run", "--track", "shared",
    "--purpose", purpose, "--config", str(CONFIG),
    "--input", f"external-compounds.parquet={CORPUS / 'external-compounds.parquet'}",
    "--input", f"observations.parquet={CORPUS / 'observations-endpoint-only.parquet'}",
    "--input", f"compounds.parquet={CHALLENGE / 'compounds.parquet'}",
    "--", "python", "scripts/train_external_medium.py",
]
result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
print(result.stdout, flush=True)
if result.stderr:
    print(result.stderr, flush=True)
result.check_returncode()
ids = re.findall(r"run_id=(\S+)", result.stdout)
if len(ids) != 1:
    raise ValueError(ids)
(ROOT / "artifacts/launches/external_shared_heads_v1_manifest.json").write_text(
    json.dumps({"config": str(CONFIG.relative_to(ROOT)), "run_id": ids[0]}, indent=2) + "\n"
)
