"""Build and train the v2 CYP2C9-specialist owned encoder."""

import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "artifacts/runs/shared/20260916T082302554619Z-656c63e6/outputs"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
SUBSET_CONFIG = ROOT / "configs/external-training/owned-encoder-2c9-only-v2.json"
TRAIN_CONFIG = ROOT / "configs/external-training/owned-encoder-2c9-train-v2.json"


def run(command: list[str]) -> str:
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    print(result.stdout, flush=True)
    if result.stderr:
        print(result.stderr, flush=True)
    result.check_returncode()
    ids = re.findall(r"run_id=(\S+)", result.stdout)
    if len(ids) != 1:
        raise ValueError(ids)
    return ids[0]


subset = run(
    [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        "Create consolidated CYP2C9-only external supervision subset",
        "--config",
        str(SUBSET_CONFIG),
        "--input",
        f"external-compounds.parquet={CORPUS / 'external-compounds.parquet'}",
        "--input",
        f"observations.parquet={CORPUS / 'observations-compact-all.parquet'}",
        "--",
        "python",
        "scripts/external_observation_subset.py",
    ]
)
subset_out = ROOT / f"artifacts/runs/shared/{subset}/outputs"
encoder = run(
    [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        "Train v2 consolidated CYP2C9-specialist width600 depth3 encoder",
        "--config",
        str(TRAIN_CONFIG),
        "--input",
        f"external-compounds.parquet={subset_out / 'external-compounds.parquet'}",
        "--input",
        f"observations.parquet={subset_out / 'observations.parquet'}",
        "--input",
        f"compounds.parquet={CHALLENGE / 'compounds.parquet'}",
        "--",
        "python",
        "scripts/train_external_medium.py",
    ]
)
encoder_out = ROOT / f"artifacts/runs/shared/{encoder}/outputs"
readout_config = ROOT / "configs/owned-encoder-v2/2c9-specialist-raw-lgbm.json"
readout_config.write_text(
    "{\n"
    '  "method": "owned_v2_2c9_specialist_raw_lgbm",\n'
    '  "backend": "lgbm",\n'
    '  "feature_file": "external-embeddings.parquet",\n'
    '  "feature_prefix": "emb_",\n'
    '  "pca": 0,\n'
    '  "mask_observed_validation": false,\n'
    '  "split_run": "20260911T104703655337Z-a3923255",\n'
    '  "endpoints": ["CYP1A2_pIC50_direct_inhibition", '
    '"CYP2C9_pIC50_direct_inhibition", "CYP2D6_pIC50_direct_inhibition", '
    '"CYP3A4_pIC50_direct_inhibition"],\n'
    f'  "encoder_run": "{encoder}"\n'
    "}\n"
)
readout = run(
    [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "All-target raw LightGBM diagnostic of v2 CYP2C9-specialist encoder",
        "--config",
        str(readout_config),
        "--input",
        f"external-embeddings.parquet={encoder_out / 'external-embeddings.parquet'}",
        "--input",
        f"measurements.parquet={CHALLENGE / 'measurements.parquet'}",
        "--input",
        f"parent-assignments.parquet={ASSIGNMENTS}",
        "--",
        "python",
        "scripts/external_embedding_readout.py",
    ]
)
print({"subset_run": subset, "encoder_run": encoder, "readout_run": readout})
