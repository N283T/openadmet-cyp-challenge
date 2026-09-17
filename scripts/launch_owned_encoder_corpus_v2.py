"""Launch reproducible curation of TDC/Novartis encoder corpus variants."""

import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "configs/external-training/owned-encoder-corpus-v2.json"
config = json.loads(CFG.read_text())
base = ROOT / f"artifacts/runs/shared/{config['base_run']}/outputs"
tdc = ROOT / f"artifacts/runs/shared/{config['tdc_run']}/outputs"
nov = ROOT / f"artifacts/runs/shared/{config['novartis_run']}/outputs"
challenge_connectivity = (
    ROOT
    / "artifacts/runs/shared/20260911T052745709989Z-8792f422/outputs"
    / "challenge-connectivity.parquet"
)

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Curate leak-excluded endpoint variants from the full TDC ADME and Polaris Novartis snapshots",
    "--config",
    str(CFG),
    "--input",
    f"external-compounds.parquet={base / 'external-compounds.parquet'}",
    "--input",
    f"observations.parquet={base / 'observations.parquet'}",
    "--input",
    f"challenge-connectivity.parquet={challenge_connectivity}",
    "--input",
    f"novartis-train-test.parquet={nov / 'train_test.parquet'}",
]
for name in config["tdc_datasets"]:
    command.extend(["--input", f"{name}.tab={tdc / (name + '.tab')}"])
command.extend(["--", "python", "scripts/curate_owned_encoder_corpus.py"])
subprocess.run(command, cwd=ROOT, check=True)
