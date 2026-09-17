"""Run selected TabPFN readouts for project-owned external encoders."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_MANIFEST = ROOT / "artifacts/launches/owned_encoder_grid_v1_manifest.json"
OUTPUT = ROOT / "artifacts/launches/owned_encoder_tabpfn_d6_v1_manifest.json"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
)
WEIGHTS = Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt"
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


base = json.loads(BASE_MANIFEST.read_text())
encoder_runs = {member["width"]: member["encoder_run"] for member in base}
screen = {256: [0], 600: [0, 128, 256]}
members = []
for width, dimensions in screen.items():
    encoder_run = encoder_runs[width]
    feature_path = (
        ROOT / f"artifacts/runs/shared/{encoder_run}/outputs/external-embeddings.parquet"
    )
    for pca in dimensions:
        representation = f"pca{pca}" if pca else f"raw{width}"
        name = f"owned_external_w{width}_s42_{representation}_tabpfn_d6_v1"
        config = {
            "method": name,
            "backend": "tabpfn",
            "feature_file": "external-embeddings.parquet",
            "feature_prefix": "emb_",
            "pca": pca,
            "weights": "tabpfn-v2.6-regressor-v2.6_default.ckpt",
            "mask_observed_validation": False,
            "split_run": "20260911T104703655337Z-a3923255",
            "endpoints": [ENDPOINT],
            "encoder_run": encoder_run,
            "encoder_width": width,
            "encoder_seed": 42,
        }
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
            f"CYP2D6 atomic TabPFN readout of owned width {width}, {representation}",
            "--config",
            str(config_path),
            "--input",
            f"external-embeddings.parquet={feature_path}",
            "--input",
            f"measurements.parquet={CHALLENGE / 'measurements.parquet'}",
            "--input",
            f"parent-assignments.parquet={ASSIGNMENTS}",
            "--input",
            f"tabpfn-v2.6-regressor-v2.6_default.ckpt={WEIGHTS}",
            "--",
            "python",
            "scripts/external_embedding_readout.py",
        ]
        run_id = recorded(command)
        members.append(
            {
                "name": name,
                "width": width,
                "pca": pca,
                "backend": "tabpfn",
                "encoder_run": encoder_run,
                "readout_run": run_id,
                "oof": str(
                    ROOT
                    / f"artifacts/runs/regression/{run_id}/outputs/external-oof.parquet"
                ),
            }
        )
        OUTPUT.write_text(json.dumps(members, indent=2) + "\n")

audit_config = {
    "endpoint": ENDPOINT,
    "scorer": "custom_scoring_functions.py",
    "protocol": "Selected independent TabPFN readouts; no blend or weight search.",
    "variants": [
        {"name": member["name"], "file": f"{member['name']}.parquet"} for member in members
    ],
}
audit_config_path = CONFIG_DIR / "owned_encoder_tabpfn_d6_audit.json"
audit_config_path.write_text(json.dumps(audit_config, indent=2) + "\n")
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Matched CYP2D6 audit of selected project-owned encoder TabPFN readouts",
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
print(json.dumps({"members": members, "audit_run": audit_run}, indent=2), flush=True)
