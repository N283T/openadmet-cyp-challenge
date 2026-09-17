"""Evaluate the focused 14-task encoder with the matched raw LightGBM readout."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
FEATURES = (
    ROOT
    / "artifacts/runs/shared/20260916T112240065253Z-6543c9c0/outputs"
    / "trial-05-w600_d3_ffn2/external-embeddings.parquet"
)
SCORER = ROOT / "artifacts/metric-reference/custom_scoring_functions.py"
READOUT_CONFIG = (
    ROOT / "configs/owned-encoder-v2/focused_cyp14_ffn2_raw600_lgbm_v1.json"
)
AUDIT_CONFIG = (
    ROOT / "configs/owned-encoder-v2/focused_cyp14_ffn2_raw600_lgbm_audit_v1.json"
)
MANIFEST = ROOT / "artifacts/launches/focused_cyp14_lgbm_v1_manifest.json"


def recorded(purpose: str, config: Path, inputs: dict[str, Path], script: str) -> str:
    command = [
        "pixi", "run", "cyp", "run", "--track", "regression",
        "--purpose", purpose, "--config", str(config),
    ]
    for alias, path in inputs.items():
        command.extend(["--input", f"{alias}={path}"])
    result = subprocess.run(
        command + ["--", "python", script], cwd=ROOT, text=True, capture_output=True
    )
    print(result.stdout, flush=True)
    if result.stderr:
        print(result.stderr, flush=True)
    result.check_returncode()
    ids = re.findall(r"run_id=(\S+)", result.stdout)
    if len(ids) != 1:
        raise ValueError(ids)
    return ids[0]


readout = recorded(
    "Matched raw LightGBM readout of the focused 14-task CYP continuous encoder",
    READOUT_CONFIG,
    {
        "external-embeddings.parquet": FEATURES,
        "measurements.parquet": CHALLENGE / "measurements.parquet",
        "parent-assignments.parquet": ASSIGNMENTS,
    },
    "scripts/external_embedding_readout.py",
)
audit = recorded(
    "Atomic audit of focused CYP14 raw LightGBM versus prior owned and Jeremy controls",
    AUDIT_CONFIG,
    {
        "custom_scoring_functions.py": SCORER,
        "focused.parquet": ROOT
        / f"artifacts/runs/regression/{readout}/outputs/external-oof.parquet",
        "previous-owned.parquet": ROOT
        / "artifacts/runs/regression/20260916T085659227651Z-a4bc5a2c/outputs"
        / "external-oof.parquet",
        "jeremy.parquet": ROOT
        / "artifacts/runs/regression/20260912T055849013223Z-f16a5dfd/outputs"
        / "external-oof.parquet",
    },
    "scripts/audit_atomic_readout_grid.py",
)
MANIFEST.write_text(json.dumps({"readout_run": readout, "audit_run": audit}, indent=2) + "\n")
print(json.dumps({"readout_run": readout, "audit_run": audit}, indent=2), flush=True)
