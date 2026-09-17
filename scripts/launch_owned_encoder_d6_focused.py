"""Train and evaluate a direct-CYP2D6-focused project-owned encoder."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "artifacts/runs/shared/20260916T073911873001Z-c23f09a3/outputs"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
)
WEIGHTS = Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt"
TRAIN_CONFIG = ROOT / "configs/external-training/owned-encoder-d6-focused-train-v1.json"
CONFIG_DIR = ROOT / "configs/owned-encoder-grid"
ENDPOINT = "CYP2D6_pIC50_direct_inhibition"


def recorded(command: list[str]) -> str:
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    print(result.stdout, flush=True)
    if result.stderr:
        print(result.stderr, flush=True)
    result.check_returncode()
    run_ids = re.findall(r"run_id=(\S+)", result.stdout)
    if len(run_ids) != 1:
        raise ValueError(f"Expected one run ID; got {run_ids}")
    return run_ids[0]


train_command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Train project-owned width-600 encoder on direct CYP2D6 external supervision",
    "--config",
    str(TRAIN_CONFIG),
    "--input",
    f"external-compounds.parquet={DATA / 'external-compounds.parquet'}",
    "--input",
    f"observations.parquet={DATA / 'observations.parquet'}",
    "--input",
    f"compounds.parquet={CHALLENGE / 'compounds.parquet'}",
    "--",
    "python",
    "scripts/train_external_medium.py",
]
encoder_run = recorded(train_command)
features = ROOT / f"artifacts/runs/shared/{encoder_run}/outputs/external-embeddings.parquet"

members = []
for backend in ["lgbm", "tabpfn"]:
    name = f"owned_external_d6focused_w600_raw600_{backend}_d6_v1"
    config = {
        "method": name,
        "backend": backend,
        "feature_file": "external-embeddings.parquet",
        "feature_prefix": "emb_",
        "pca": 0,
        "mask_observed_validation": False,
        "split_run": "20260911T104703655337Z-a3923255",
        "endpoints": [ENDPOINT],
        "encoder_run": encoder_run,
        "encoder_width": 600,
        "encoder_seed": 42,
        "data_variant": "direct_CYP2D6_only",
    }
    inputs = {
        "external-embeddings.parquet": features,
        "measurements.parquet": CHALLENGE / "measurements.parquet",
        "parent-assignments.parquet": ASSIGNMENTS,
    }
    if backend == "tabpfn":
        config["weights"] = "tabpfn-v2.6-regressor-v2.6_default.ckpt"
        inputs["tabpfn-v2.6-regressor-v2.6_default.ckpt"] = WEIGHTS
    config_path = CONFIG_DIR / f"{name}.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"CYP2D6 atomic {backend} readout of direct-D6-focused owned raw600 encoder",
        "--config",
        str(config_path),
    ]
    for alias, path in inputs.items():
        command.extend(["--input", f"{alias}={path}"])
    command.extend(["--", "python", "scripts/external_embedding_readout.py"])
    run_id = recorded(command)
    members.append(
        {
            "name": name,
            "backend": backend,
            "run_id": run_id,
            "oof": str(
                ROOT / f"artifacts/runs/regression/{run_id}/outputs/external-oof.parquet"
            ),
        }
    )

audit_config = {
    "endpoint": ENDPOINT,
    "scorer": "custom_scoring_functions.py",
    "protocol": "Independent focused readouts only; no blend or weight search.",
    "variants": [
        {"name": member["name"], "file": f"{member['name']}.parquet"} for member in members
    ],
}
audit_config_path = CONFIG_DIR / "owned_encoder_d6_focused_audit.json"
audit_config_path.write_text(json.dumps(audit_config, indent=2) + "\n")
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Matched CYP2D6 audit of direct-D6-focused project-owned encoder readouts",
    "--config",
    str(audit_config_path),
    "--input",
    "custom_scoring_functions.py="
    + str(ROOT / "artifacts/metric-reference/custom_scoring_functions.py"),
]
for member in members:
    command.extend(["--input", f"{member['name']}.parquet={member['oof']}"])
command.extend(["--", "python", "scripts/audit_atomic_readout_grid.py"])
audit_run = recorded(command)
manifest = {"encoder_run": encoder_run, "members": members, "audit_run": audit_run}
(ROOT / "artifacts/launches/owned_encoder_d6_focused_v1_manifest.json").write_text(
    json.dumps(manifest, indent=2) + "\n"
)
print(json.dumps(manifest, indent=2), flush=True)
