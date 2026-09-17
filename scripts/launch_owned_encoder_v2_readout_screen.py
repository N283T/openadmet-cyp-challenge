"""Screen LightGBM capacity and PCA on the strongest owned representations."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRID_PATH = ROOT / "configs/external-training/owned-encoder-v2-readout-screen.json"
GRID = json.loads(GRID_PATH.read_text())
CONFIG_DIR = ROOT / "configs/owned-encoder-v2"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
SCORER = ROOT / "artifacts/metric-reference/custom_scoring_functions.py"
MANIFEST = ROOT / "artifacts/launches/owned_encoder_v2_readout_manifest.json"
ENDPOINTS = [
    "CYP1A2_pIC50_direct_inhibition",
    "CYP2C9_pIC50_direct_inhibition",
    "CYP2D6_pIC50_direct_inhibition",
    "CYP3A4_pIC50_direct_inhibition",
]


def recorded(purpose: str, config: Path, inputs: dict[str, Path], script: str):
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
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
for feature_name, feature in GRID["features"].items():
    feature_path = ROOT / f"artifacts/runs/shared/{feature['run']}/outputs/{feature['file']}"
    arms = [(name, 0, overrides) for name, overrides in GRID["lgbm_grid"].items()]
    arms += [(f"pca{pca}", pca, {}) for pca in GRID["pca_dimensions"]]
    for arm, pca, overrides in arms:
        name = f"{feature_name}_{arm}"
        config = {
            "method": f"owned_v2_{name}_lgbm",
            "backend": "lgbm",
            "feature_file": "concatenated-embeddings.parquet",
            "feature_prefix": "emb_",
            "pca": pca,
            "mask_observed_validation": False,
            "split_run": "20260911T104703655337Z-a3923255",
            "endpoints": ENDPOINTS,
            "feature_run": feature["run"],
            "lgbm_overrides": overrides,
        }
        config_path = CONFIG_DIR / f"readout_{name}.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n")
        run = recorded(
            f"Owned v2 {feature_name} readout screen: {arm}",
            config_path,
            {
                "concatenated-embeddings.parquet": feature_path,
                "measurements.parquet": CHALLENGE / "measurements.parquet",
                "parent-assignments.parquet": ASSIGNMENTS,
            },
            "scripts/external_embedding_readout.py",
        )
        members.append({"name": name, "readout_run": run})
        MANIFEST.write_text(json.dumps(members, indent=2) + "\n")

controls = {
    "top_arch_default": "20260916T091505235124Z-4d54014a",
    "broad5_default": "20260916T091747327774Z-3cf2e8bb",
    "jeremy_raw600_lgbm": "20260912T055849013223Z-f16a5dfd",
}
audit = {
    "endpoints": ENDPOINTS,
    "scorer": "custom_scoring_functions.py",
    "protocol": "Target-specific readout screen; no blending or bank insertion.",
    "variants": [{"name": item["name"], "file": f"{item['name']}.parquet"} for item in members]
    + [{"name": name, "file": f"{name}.parquet"} for name in controls],
}
audit_path = CONFIG_DIR / "readout-screen-audit.json"
audit_path.write_text(json.dumps(audit, indent=2) + "\n")
inputs = {"custom_scoring_functions.py": SCORER}
for item in members:
    inputs[f"{item['name']}.parquet"] = (
        ROOT / f"artifacts/runs/regression/{item['readout_run']}/outputs/external-oof.parquet"
    )
for name, run in controls.items():
    inputs[f"{name}.parquet"] = (
        ROOT / f"artifacts/runs/regression/{run}/outputs/external-oof.parquet"
    )
audit_run = recorded(
    "Matched all-target audit of owned v2 readout screen",
    audit_path,
    inputs,
    "scripts/audit_atomic_readout_grid.py",
)
result = {"members": members, "controls": controls, "audit_run": audit_run}
MANIFEST.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2), flush=True)
