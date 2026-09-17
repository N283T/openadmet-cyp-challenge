"""Train corpus variants and evaluate frozen raw embeddings on all regression targets."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRID_PATH = ROOT / "configs/external-training/owned-encoder-v2-data-screen.json"
GRID = json.loads(GRID_PATH.read_text())
CORPUS = ROOT / f"artifacts/runs/shared/{GRID['corpus_run']}/outputs"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
CONFIG_DIR = ROOT / "configs/owned-encoder-v2"
CONFIG_DIR.mkdir(exist_ok=True)
MANIFEST = ROOT / "artifacts/launches/owned_encoder_v2_data_screen_manifest.json"


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


manifest = []
for variant in GRID["variants"]:
    tag = f"owned_v2_{variant}_w600_d4_norm"
    train_config = {
        key: value
        for key, value in GRID.items()
        if key not in {"corpus_run", "variants", "readout", "protocol"}
    }
    train_config["purpose"] = GRID["protocol"]
    train_path = CONFIG_DIR / f"{tag}_train.json"
    train_path.write_text(json.dumps(train_config, indent=2) + "\n")
    encoder = recorded(
        "shared",
        f"Owned encoder v2 data screen: {variant}, width600 depth4 norm",
        train_path,
        {
            "external-compounds.parquet": CORPUS / "external-compounds.parquet",
            "observations.parquet": CORPUS / f"observations-{variant.replace('_', '-')}.parquet",
            "compounds.parquet": CHALLENGE / "compounds.parquet",
        },
        "scripts/train_external_medium.py",
    )
    readout_config = {
        "method": f"{tag}_raw_lgbm",
        "backend": "lgbm",
        "feature_file": "external-embeddings.parquet",
        "feature_prefix": "emb_",
        "pca": 0,
        "mask_observed_validation": False,
        "split_run": "20260911T104703655337Z-a3923255",
        "endpoints": [
            "CYP1A2_pIC50_direct_inhibition",
            "CYP2C9_pIC50_direct_inhibition",
            "CYP2D6_pIC50_direct_inhibition",
            "CYP3A4_pIC50_direct_inhibition",
        ],
        "encoder_run": encoder,
        "corpus_variant": variant,
    }
    readout_path = CONFIG_DIR / f"{tag}_raw_lgbm.json"
    readout_path.write_text(json.dumps(readout_config, indent=2) + "\n")
    readout = recorded(
        "regression",
        f"All-target raw LightGBM readout of owned encoder v2 {variant}",
        readout_path,
        {
            "external-embeddings.parquet": ROOT
            / f"artifacts/runs/shared/{encoder}/outputs/external-embeddings.parquet",
            "measurements.parquet": CHALLENGE / "measurements.parquet",
            "parent-assignments.parquet": ASSIGNMENTS,
        },
        "scripts/external_embedding_readout.py",
    )
    manifest.append({"variant": variant, "encoder_run": encoder, "readout_run": readout})
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")

print(json.dumps(manifest, indent=2), flush=True)
