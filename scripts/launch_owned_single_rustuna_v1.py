"""Run single-checkpoint Rustuna search, full OOF readout and Jeremy audit."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
corpus = root / "artifacts/runs/shared/20260916T082302554619Z-656c63e6/outputs"
challenge = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
assignments = (
    root
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
jeremy = (
    root
    / "artifacts/runs/regression/20260912T055849013223Z-f16a5dfd/outputs"
    / "external-oof.parquet"
)
control = (
    root
    / "artifacts/runs/regression/20260916T085659227651Z-a4bc5a2c/outputs"
    / "external-oof.parquet"
)
scorer = root / "artifacts/metric-reference/custom_scoring_functions.py"
hpo_config = root / "configs/external-training/owned-single-rustuna-v1.json"
readout_config = root / "configs/owned-encoder-v2/owned_single_rustuna_v1_raw_lgbm.json"
audit_config = root / "configs/owned-encoder-v2/owned_single_rustuna_v1_audit.json"
manifest_path = root / "artifacts/launches/owned_single_rustuna_v1_manifest.json"


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
        command + ["--", "python", script], cwd=root, capture_output=True, text=True
    )
    print(result.stdout, flush=True)
    if result.stderr:
        print(result.stderr, flush=True)
    result.check_returncode()
    run_ids = re.findall(r"run_id=(\S+)", result.stdout)
    if len(run_ids) != 1:
        raise ValueError(run_ids)
    return run_ids[0]


hpo_run = recorded(
    "regression",
    "Rustuna search for one owned encoder checkpoint; fold 3 tunes and fold 4 confirms",
    hpo_config,
    {
        "external-compounds.parquet": corpus / "external-compounds.parquet",
        "observations.parquet": corpus / "observations-endpoint-only.parquet",
        "compounds.parquet": challenge / "compounds.parquet",
        "measurements.parquet": challenge / "measurements.parquet",
        "parent-assignments.parquet": assignments,
        "custom_scoring_functions.py": scorer,
        "jeremy-oof.parquet": jeremy,
    },
    "scripts/tune_owned_single_encoder.py",
)
hpo_output = root / f"artifacts/runs/regression/{hpo_run}/outputs"
readout_run = recorded(
    "regression",
    "All-fold matched raw LightGBM for the Rustuna-selected single owned encoder",
    readout_config,
    {
        "selected-embeddings.parquet": hpo_output / "selected-embeddings.parquet",
        "measurements.parquet": challenge / "measurements.parquet",
        "parent-assignments.parquet": assignments,
    },
    "scripts/external_embedding_readout.py",
)
readout_output = root / f"artifacts/runs/regression/{readout_run}/outputs"
audit_run = recorded(
    "regression",
    "Matched audit of Rustuna-selected single owned encoder versus Jeremy and d3 control",
    audit_config,
    {
        "custom_scoring_functions.py": scorer,
        "owned.parquet": readout_output / "external-oof.parquet",
        "jeremy.parquet": jeremy,
        "owned-d3-control.parquet": control,
    },
    "scripts/audit_atomic_readout_grid.py",
)
manifest = {"hpo_run": hpo_run, "readout_run": readout_run, "audit_run": audit_run}
manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps(manifest, indent=2), flush=True)
