"""Run the v2 architecture screen and matched atomic audit."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRID_PATH = ROOT / "configs/external-training/owned-encoder-v2-architecture-screen.json"
GRID = json.loads(GRID_PATH.read_text())
CORPUS = ROOT / f"artifacts/runs/shared/{GRID['corpus_run']}/outputs"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
SCORER = ROOT / "artifacts/metric-reference/custom_scoring_functions.py"
CONFIG_DIR = ROOT / "configs/owned-encoder-v2"
CONFIG_DIR.mkdir(exist_ok=True)
MANIFEST = ROOT / "artifacts/launches/owned_encoder_v2_architecture_manifest.json"
ENDPOINTS = [
    "CYP1A2_pIC50_direct_inhibition",
    "CYP2C9_pIC50_direct_inhibition",
    "CYP2D6_pIC50_direct_inhibition",
    "CYP3A4_pIC50_direct_inhibition",
]


def recorded(track: str, purpose: str, config: Path, inputs: dict[str, Path], script: str):
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        track,
        "--purpose",
        purpose,
        "--config",
        str(config),
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


def tag(arch: dict) -> str:
    return f"w{arch['width']}_d{arch['depth']}_{arch['aggregation']}_f{arch['ffn_hidden_layers']}"


members = []
baseline_tag = "w600_d4_norm_f3"
members.append({"name": baseline_tag, **GRID["reuse"][baseline_tag]})
for arch in GRID["architectures"]:
    name = tag(arch)
    train = {
        key: value
        for key, value in GRID.items()
        if key not in {"corpus_run", "variant", "reuse", "architectures", "protocol"}
    }
    train.update(
        architecture="owned_width_grid",
        aggregation_norm=100.0,
        purpose=GRID["protocol"],
        **arch,
    )
    train_path = CONFIG_DIR / f"owned_v2_endpoint_{name}_train.json"
    train_path.write_text(json.dumps(train, indent=2) + "\n")
    encoder = recorded(
        "shared",
        f"Owned v2 endpoint-only architecture screen {name}",
        train_path,
        {
            "external-compounds.parquet": CORPUS / "external-compounds.parquet",
            "observations.parquet": CORPUS / "observations-endpoint-only.parquet",
            "compounds.parquet": CHALLENGE / "compounds.parquet",
        },
        "scripts/train_external_medium.py",
    )
    readout_config = {
        "method": f"owned_v2_endpoint_{name}_raw_lgbm",
        "backend": "lgbm",
        "feature_file": "external-embeddings.parquet",
        "feature_prefix": "emb_",
        "pca": 0,
        "mask_observed_validation": False,
        "split_run": "20260911T104703655337Z-a3923255",
        "endpoints": ENDPOINTS,
        "encoder_run": encoder,
        "architecture_tag": name,
    }
    readout_path = CONFIG_DIR / f"owned_v2_endpoint_{name}_raw_lgbm.json"
    readout_path.write_text(json.dumps(readout_config, indent=2) + "\n")
    readout = recorded(
        "regression",
        f"All-target raw LightGBM readout of owned v2 endpoint-only {name}",
        readout_path,
        {
            "external-embeddings.parquet": ROOT
            / f"artifacts/runs/shared/{encoder}/outputs/external-embeddings.parquet",
            "measurements.parquet": CHALLENGE / "measurements.parquet",
            "parent-assignments.parquet": ASSIGNMENTS,
        },
        "scripts/external_embedding_readout.py",
    )
    members.append({"name": name, "encoder_run": encoder, "readout_run": readout})
    MANIFEST.write_text(json.dumps(members, indent=2) + "\n")

controls = {
    "v2_compact_all_w600": "20260916T083209317276Z-d8372328",
    "v2_assay_w600": "20260916T082842598129Z-546bf741",
    "jeremy_raw600_lgbm": "20260912T055849013223Z-f16a5dfd",
}
variants = [{"name": member["name"], "file": f"{member['name']}.parquet"} for member in members] + [
    {"name": name, "file": f"{name}.parquet"} for name in controls
]
audit = {
    "endpoints": ENDPOINTS,
    "scorer": "custom_scoring_functions.py",
    "protocol": "Matched raw LightGBM architecture screen; no blending or bank insertion.",
    "variants": variants,
}
audit_path = CONFIG_DIR / "architecture-screen-audit.json"
audit_path.write_text(json.dumps(audit, indent=2) + "\n")
inputs = {"custom_scoring_functions.py": SCORER}
for member in members:
    inputs[f"{member['name']}.parquet"] = (
        ROOT / f"artifacts/runs/regression/{member['readout_run']}/outputs/external-oof.parquet"
    )
for name, run in controls.items():
    inputs[f"{name}.parquet"] = (
        ROOT / f"artifacts/runs/regression/{run}/outputs/external-oof.parquet"
    )
audit_run = recorded(
    "regression",
    "Matched all-target audit of owned encoder v2 architecture screen",
    audit_path,
    inputs,
    "scripts/audit_atomic_readout_grid.py",
)
result = {"members": members, "controls": controls, "audit_run": audit_run}
MANIFEST.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2), flush=True)
