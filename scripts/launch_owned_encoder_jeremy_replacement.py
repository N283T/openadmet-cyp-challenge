"""Evaluate project-owned encoders against Jeremy frozen features on all targets."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
)
WEIGHTS = Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt"
SCORER = ROOT / "artifacts/metric-reference/custom_scoring_functions.py"
CONFIG_DIR = ROOT / "configs/owned-encoder-grid"
OUTPUT = ROOT / "artifacts/launches/owned_encoder_jeremy_replacement_v1_manifest.json"
ENDPOINTS = [
    "CYP1A2_pIC50_direct_inhibition",
    "CYP2C9_pIC50_direct_inhibition",
    "CYP2D6_pIC50_direct_inhibition",
    "CYP3A4_pIC50_direct_inhibition",
]
ENCODERS = {
    256: "20260916T073258258609Z-c05c3377",
    600: "20260916T073416372035Z-2dda20d9",
}
JEREMY = [
    {
        "name": "jeremy_raw600_lgbm",
        "file": ROOT
        / "artifacts/runs/regression/20260912T055849013223Z-f16a5dfd/outputs"
        / "external-oof.parquet",
    },
    {
        "name": "jeremy_raw600_tabpfn",
        "file": ROOT
        / "artifacts/runs/regression/20260912T060059705105Z-2f0fbfbd/outputs"
        / "external-oof.parquet",
    },
]


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


members = []
representations = [(256, 0), (600, 0), (600, 256)]
for width, pca in representations:
    encoder_run = ENCODERS[width]
    feature_path = (
        ROOT / f"artifacts/runs/shared/{encoder_run}/outputs/external-embeddings.parquet"
    )
    representation = f"pca{pca}" if pca else f"raw{width}"
    for backend in ["lgbm", "tabpfn"]:
        name = f"owned_w{width}_{representation}_{backend}_alltargets_v1"
        config = {
            "method": name,
            "backend": backend,
            "feature_file": "external-embeddings.parquet",
            "feature_prefix": "emb_",
            "pca": pca,
            "mask_observed_validation": False,
            "split_run": "20260911T104703655337Z-a3923255",
            "endpoints": ENDPOINTS,
            "encoder_run": encoder_run,
            "encoder_width": width,
            "encoder_seed": 42,
            "comparison_role": "Jeremy frozen encoder replacement",
        }
        inputs = {
            "external-embeddings.parquet": feature_path,
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
            f"All-target {backend} Jeremy-replacement test of owned {representation} encoder",
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
                "width": width,
                "pca": pca,
                "backend": backend,
                "run_id": run_id,
                "oof": str(
                    ROOT / f"artifacts/runs/regression/{run_id}/outputs/external-oof.parquet"
                ),
            }
        )
        OUTPUT.write_text(json.dumps(members, indent=2) + "\n")

audit_variants = [
    {"name": member["name"], "file": f"{member['name']}.parquet"} for member in members
] + [{"name": item["name"], "file": f"{item['name']}.parquet"} for item in JEREMY]
audit_config = {
    "endpoints": ENDPOINTS,
    "scorer": "custom_scoring_functions.py",
    "protocol": (
        "Matched atomic comparison against Jeremy frozen raw600. No blending, "
        "bank insertion, or weight search."
    ),
    "variants": audit_variants,
}
audit_path = CONFIG_DIR / "owned_encoder_jeremy_replacement_audit.json"
audit_path.write_text(json.dumps(audit_config, indent=2) + "\n")
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Matched all-target atomic audit of project-owned encoders versus Jeremy frozen raw600",
    "--config",
    str(audit_path),
    "--input",
    f"custom_scoring_functions.py={SCORER}",
]
for member in members:
    command.extend(["--input", f"{member['name']}.parquet={member['oof']}"])
for item in JEREMY:
    command.extend(["--input", f"{item['name']}.parquet={item['file']}"])
command.extend(["--", "python", "scripts/audit_atomic_readout_grid.py"])
audit_run = recorded(command)
result = {"members": members, "jeremy_controls": JEREMY, "audit_run": audit_run}
OUTPUT.write_text(json.dumps(result, indent=2, default=str) + "\n")
print(json.dumps(result, indent=2, default=str), flush=True)
