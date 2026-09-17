"""Test whether the CYP2C9-specialist embedding complements owned broad5."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_CONFIG = (
    ROOT / "configs/external-training/owned-encoder-v2-2c9-specialist-concat.json"
)
CONFIG_DIR = ROOT / "configs/owned-encoder-v2"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
SCORER = ROOT / "artifacts/metric-reference/custom_scoring_functions.py"
ENDPOINT = "CYP2C9_pIC50_direct_inhibition"


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


source = json.loads(SOURCE_CONFIG.read_text())
concat_run = recorded(
    "shared",
    "Append owned CYP2C9-specialist raw embedding to owned broad5",
    SOURCE_CONFIG,
    {
        "broad5.parquet": ROOT
        / "artifacts/runs/shared/20260916T091746252797Z-e571bd53/outputs"
        / "concatenated-embeddings.parquet",
        "2c9-specialist.parquet": ROOT
        / "artifacts/runs/shared/20260916T094226751324Z-20c2d3d6/outputs"
        / "external-embeddings.parquet",
    },
    "scripts/concat_embedding_views.py",
)

variants = {
    "broad5_plus_2c9_specialist_default": {},
    "broad5_plus_2c9_specialist_leaves31": {
        "n_estimators": 500,
        "num_leaves": 31,
        "min_child_samples": 20,
        "learning_rate": 0.03,
    },
}
readouts = {}
for name, overrides in variants.items():
    config = {
        "method": name,
        "backend": "lgbm",
        "feature_file": "concatenated-embeddings.parquet",
        "feature_prefix": "emb_",
        "pca": 0,
        "mask_observed_validation": False,
        "split_run": "20260911T104703655337Z-a3923255",
        "endpoints": [ENDPOINT],
        "concat_run": concat_run,
        "lgbm_overrides": overrides,
    }
    path = CONFIG_DIR / f"{name}.json"
    path.write_text(json.dumps(config, indent=2) + "\n")
    readouts[name] = recorded(
        "regression",
        f"CYP2C9-only readout for {name}",
        path,
        {
            "concatenated-embeddings.parquet": ROOT
            / f"artifacts/runs/shared/{concat_run}/outputs/concatenated-embeddings.parquet",
            "measurements.parquet": CHALLENGE / "measurements.parquet",
            "parent-assignments.parquet": ASSIGNMENTS,
        },
        "scripts/external_embedding_readout.py",
    )

controls = {
    "broad5_default": "20260916T091747327774Z-3cf2e8bb",
    "broad5_leaves31": "20260916T093341774576Z-6dedf75d",
    "jeremy_raw600_lgbm": "20260912T055849013223Z-f16a5dfd",
    "2c9_specialist_alone": "20260916T094416542425Z-8d87951f",
}
audit = {
    "endpoint": ENDPOINT,
    "scorer": "custom_scoring_functions.py",
    "protocol": source["protocol"],
    "variants": [
        {"name": name, "file": f"{name}.parquet"}
        for name in [*readouts, *controls]
    ],
}
audit_path = CONFIG_DIR / "2c9-specialist-concat-audit.json"
audit_path.write_text(json.dumps(audit, indent=2) + "\n")
audit_inputs = {"custom_scoring_functions.py": SCORER}
for name, run_id in {**readouts, **controls}.items():
    audit_inputs[f"{name}.parquet"] = (
        ROOT / f"artifacts/runs/regression/{run_id}/outputs/external-oof.parquet"
    )
audit_run = recorded(
    "regression",
    "Matched CYP2C9 audit of specialist embedding appended to broad5",
    audit_path,
    audit_inputs,
    "scripts/audit_atomic_readout_grid.py",
)
print(
    json.dumps(
        {"concat_run": concat_run, "readouts": readouts, "audit_run": audit_run},
        indent=2,
    ),
    flush=True,
)
