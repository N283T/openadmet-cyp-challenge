"""Launch fixed paired readouts of old and alternative primary features."""

import argparse
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
base = root / "artifacts/runs"
parser = argparse.ArgumentParser()
parser.add_argument("--mode", choices=["old", "new", "both"], default="both")
args = parser.parse_args()
for mode, run in [
    ("old", "20260912T052549384296Z-bac5fb3a"),
    ("new", "20260912T082554377951Z-d272774c"),
]:
    if args.mode not in [mode, "both"]:
        continue
    assert (base / f"shared/{run}/result.json").exists(), "Wait for primary run completion"
    cfg = dict(
        method=f"focus_primary_{mode}_lgbm_v1",
        backend="lgbm",
        feature_file="external-embeddings.parquet",
        feature_prefix="emb_",
        nested_file="nested-log2fc.parquet",
        split_run="20260911T104703655337Z-a3923255",
        lgbm_overrides={"n_jobs": 3},
    )
    p = root / f"configs/focus-primary/readout-{mode}-v1.json"
    p.write_text(json.dumps(cfg, indent=2))
    inputs = [
        base / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
        base / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
        base / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
        base / f"shared/{run}/outputs/nested-log2fc.parquet",
    ]
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Fixed600emb pIC50LGBM paired control with {mode} nested primary prediction",
        "--config",
        str(p),
    ]
    for f in inputs:
        command.extend(["--input", str(f)])
    r = subprocess.run(
        command + ["--", "python", "scripts/external_embedding_readout.py"],
        cwd=root,
        text=True,
        capture_output=True,
    )
    print(r.stdout, r.stderr, flush=True)
    r.check_returncode()
    rid = re.search(r"run_id=(\S+)", r.stdout)[1]
    acfg = {
        "model_run": rid,
        "database": {
            "host": str(root / ".runtime/socket"),
            "port": 5434,
            "dbname": "cyp_challenge",
        },
    }
    ap = root / f"configs/focus-primary/audit-{mode}-v1.json"
    ap.write_text(json.dumps(acfg, indent=2))
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Audit focus primary {mode} readout",
        "--config",
        str(ap),
    ]
    for f in [
        base / f"regression/{rid}/outputs/external-oof.parquet",
        base / "regression/20260911T153808228197Z-dec24994/outputs/aligned-oof.parquet",
        root / "artifacts/metric-reference/custom_scoring_functions.py",
    ]:
        command.extend(["--input", str(f)])
    r = subprocess.run(
        command + ["--", "python", "scripts/audit_external_readout.py"],
        cwd=root,
        text=True,
        capture_output=True,
    )
    print(r.stdout, r.stderr, flush=True)
    r.check_returncode()
