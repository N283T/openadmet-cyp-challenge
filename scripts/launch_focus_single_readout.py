"""Fixed LGBM+safe nestedlogfc readouts and25/50%ensemble controls."""

import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
split = "20260911T104703655337Z-a3923255"
data = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
manifest = {}


def launch(name, cfg, inputs, script):
    path = root / f"configs/focus-single-target/{name}.json"
    with path.open("x") as stream:
        json.dump(cfg, stream, indent=2)
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        name,
        "--config",
        str(path),
    ]
    for item in inputs:
        command += ["--input", str(item)]
    result = subprocess.run(
        command + ["--", "python", script], cwd=root, capture_output=True, text=True
    )
    (root / f"artifacts/launches/{name}.log").write_text(result.stdout + result.stderr)
    print(result.stdout, flush=True)
    result.check_returncode()
    runs = re.findall(r"run_id=(\S+)", result.stdout)
    assert len(runs) == 1
    return runs[0]


for target, fit, n in [
    ("CYP1A2", "20260912T082204145061Z-26558f40", 1372),
    ("CYP2D6", "20260912T082308319207Z-5c3b657c", 1457),
]:
    source = root / f"artifacts/runs/regression/{fit}/outputs"
    inputs = [source / f"fold{fold}/external-embeddings.parquet" for fold in range(5)] + [
        data / "measurements.parquet",
        root / f"artifacts/runs/regression/{split}/outputs/parent-assignments.parquet",
        root
        / "artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs/nested-log2fc.parquet",
    ]
    cfg = dict(
        endpoint=target + "_pIC50_direct_inhibition",
        method="focus_single_" + target + "_lgbm_logfc",
        backend="lgbm",
        split_run=split,
        expected_oof_rows=n,
        feature_files={str(f): f"{f:03d}-external-embeddings.parquet" for f in range(5)},
        nested_file="nested-log2fc.parquet",
        lgbm_overrides=dict(n_jobs=2),
    )
    run = launch(
        "focus_single_" + target + "_lgbm_logfc_v1", cfg, inputs, "scripts/focus_single_readout.py"
    )
    audit_inputs = [
        root / f"artifacts/runs/regression/{run}/outputs/external-oof.parquet",
        root
        / "artifacts/runs/regression/20260912T072935321335Z-1e8d2643/outputs/external-oof.parquet",
        root / "artifacts/metric-reference/custom_scoring_functions.py",
    ] + [source / f"fold{fold}/membership.parquet" for fold in range(5)]
    audit = launch(
        "focus_single_" + target + "_lgbm_logfc_audit_v1",
        dict(endpoint=target, n=n, prediction_file="000-external-oof.parquet"),
        audit_inputs,
        "scripts/audit_focus_single_target.py",
    )
    manifest[target] = dict(fit=fit, readout=run, audit=audit)
    (root / "artifacts/launches/focus_single_readout_manifest.json").write_text(
        json.dumps(manifest, indent=2)
    )
