"""Build and evaluate concatenated owned embedding views."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRID_PATH = ROOT / "configs/external-training/owned-encoder-v2-concat-screen.json"
GRID = json.loads(GRID_PATH.read_text())
CONFIG_DIR = ROOT / "configs/owned-encoder-v2"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
SCORER = ROOT / "artifacts/metric-reference/custom_scoring_functions.py"
MANIFEST = ROOT / "artifacts/launches/owned_encoder_v2_concat_manifest.json"
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


members = []
for recipe, names in GRID["recipes"].items():
    concat_config = {
        "recipe": recipe,
        "protocol": GRID["protocol"],
        "views": [
            {"name": name, "file": f"{name}.parquet", **GRID["views"][name]} for name in names
        ],
    }
    concat_path = CONFIG_DIR / f"concat_{recipe}.json"
    concat_path.write_text(json.dumps(concat_config, indent=2) + "\n")
    concat_inputs = {
        f"{name}.parquet": ROOT
        / f"artifacts/runs/shared/{GRID['views'][name]['run']}/outputs/external-embeddings.parquet"
        for name in names
    }
    concat_run = recorded(
        "shared",
        f"Concatenate owned raw embedding views: {recipe}",
        concat_path,
        concat_inputs,
        "scripts/concat_embedding_views.py",
    )
    readout_config = {
        "method": f"owned_v2_concat_{recipe}_lgbm",
        "backend": "lgbm",
        "feature_file": "concatenated-embeddings.parquet",
        "feature_prefix": "emb_",
        "pca": 0,
        "mask_observed_validation": False,
        "split_run": "20260911T104703655337Z-a3923255",
        "endpoints": ENDPOINTS,
        "concat_run": concat_run,
        "recipe": recipe,
    }
    readout_path = CONFIG_DIR / f"concat_{recipe}_lgbm.json"
    readout_path.write_text(json.dumps(readout_config, indent=2) + "\n")
    readout_run = recorded(
        "regression",
        f"All-target LightGBM on concatenated owned raw views: {recipe}",
        readout_path,
        {
            "concatenated-embeddings.parquet": ROOT
            / f"artifacts/runs/shared/{concat_run}/outputs/concatenated-embeddings.parquet",
            "measurements.parquet": CHALLENGE / "measurements.parquet",
            "parent-assignments.parquet": ASSIGNMENTS,
        },
        "scripts/external_embedding_readout.py",
    )
    members.append({"name": recipe, "concat_run": concat_run, "readout_run": readout_run})
    MANIFEST.write_text(json.dumps(members, indent=2) + "\n")

controls = {
    "owned_endpoint_d3": "20260916T085659227651Z-a4bc5a2c",
    "owned_endpoint_f2": "20260916T091257666544Z-b2ab5467",
    "owned_compact_all": "20260916T083209317276Z-d8372328",
    "jeremy_raw600_lgbm": "20260912T055849013223Z-f16a5dfd",
}
audit = {
    "endpoints": ENDPOINTS,
    "scorer": "custom_scoring_functions.py",
    "protocol": "Matched raw feature-level owned-view concatenation screen.",
    "variants": [{"name": item["name"], "file": f"{item['name']}.parquet"} for item in members]
    + [{"name": name, "file": f"{name}.parquet"} for name in controls],
}
audit_path = CONFIG_DIR / "concat-screen-audit.json"
audit_path.write_text(json.dumps(audit, indent=2) + "\n")
audit_inputs = {"custom_scoring_functions.py": SCORER}
for item in members:
    audit_inputs[f"{item['name']}.parquet"] = (
        ROOT / f"artifacts/runs/regression/{item['readout_run']}/outputs/external-oof.parquet"
    )
for name, run in controls.items():
    audit_inputs[f"{name}.parquet"] = (
        ROOT / f"artifacts/runs/regression/{run}/outputs/external-oof.parquet"
    )
audit_run = recorded(
    "regression",
    "Matched all-target audit of concatenated owned raw embedding views",
    audit_path,
    audit_inputs,
    "scripts/audit_atomic_readout_grid.py",
)
result = {"members": members, "controls": controls, "audit_run": audit_run}
MANIFEST.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2), flush=True)
