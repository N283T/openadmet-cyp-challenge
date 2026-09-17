"""Run and audit CYP2D6 raw/PCA atomic readouts of ATTNSOM embeddings."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRID_PATH = ROOT / "configs/regression/attnsom-atomic-grid-d6-v1.json"
GRID = json.loads(GRID_PATH.read_text())
FEATURES = (
    ROOT
    / f"artifacts/runs/classification/{GRID['embedding_run']}/outputs/attnsom-embeddings.parquet"
)
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
)
SCORER = ROOT / "artifacts/metric-reference/custom_scoring_functions.py"
CONFIG_DIR = ROOT / "configs/attnsom-atomic-grid"
CONFIG_DIR.mkdir(exist_ok=True)


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


variants = []
for backend in GRID["backends"]:
    for dimensions in GRID["dimensions"]:
        dimension_name = "raw512" if dimensions == 0 else f"pca{dimensions}"
        name = f"attnsom_{dimension_name}_{backend}_d6"
        config = {
            "method": name,
            "backend": backend,
            "feature_file": "attnsom-embeddings.parquet",
            "feature_prefix": "attnsom_emb_",
            "pca": dimensions,
            "mask_observed_validation": False,
            "split_run": "20260911T104703655337Z-a3923255",
            "endpoints": [GRID["endpoint"]],
            "source_training_run": GRID["training_run"],
            "source_embedding_run": GRID["embedding_run"],
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
            f"CYP2D6 atomic {backend} readout of ATTNSOM {dimension_name} embedding",
            "--config",
            str(config_path),
            "--input",
            f"attnsom-embeddings.parquet={FEATURES}",
            "--input",
            f"measurements.parquet={CHALLENGE / 'measurements.parquet'}",
            "--input",
            f"parent-assignments.parquet={ASSIGNMENTS}",
            "--",
            "python",
            "scripts/external_embedding_readout.py",
        ]
        run_id = recorded(command)
        variants.append(
            {
                "name": name,
                "run_id": run_id,
                "path": ROOT
                / f"artifacts/runs/regression/{run_id}/outputs/external-oof.parquet",
            }
        )

audit_config = {
    **GRID,
    "scorer": "custom_scoring_functions.py",
    "variants": [
        {"name": variant["name"], "file": f"{variant['name']}.parquet"}
        for variant in variants
    ],
}
audit_config_path = CONFIG_DIR / "attnsom_atomic_grid_d6_audit.json"
audit_config_path.write_text(json.dumps(audit_config, indent=2) + "\n")
audit_command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Matched CYP2D6 audit of independent ATTNSOM raw and PCA atomic readouts",
    "--config",
    str(audit_config_path),
    "--input",
    f"custom_scoring_functions.py={SCORER}",
]
for variant in variants:
    audit_command.extend(["--input", f"{variant['name']}.parquet={variant['path']}"])
audit_command.extend(["--", "python", "scripts/audit_atomic_readout_grid.py"])
audit_run = recorded(audit_command)
manifest = {
    "readouts": [{**variant, "path": str(variant["path"])} for variant in variants],
    "audit_run": audit_run,
}
(ROOT / "artifacts/launches/attnsom_atomic_grid_d6_v1_manifest.json").write_text(
    json.dumps(manifest, indent=2) + "\n"
)
print(json.dumps(manifest, indent=2), flush=True)
