"""Launch a recorded CYP/HEME/ligand trunk-input preparation run."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

REFERENCE_RUNS = {
    "cyp1a2": "20260928T111653134056Z-84ffca6a",
    "cyp2c9": "20260928T111758347957Z-0357cb7c",
    "cyp2d6": "20260928T111758358079Z-0d42000d",
    "cyp3a4": "20260928T111713518371Z-c156bd1a",
}
ENSEMBLE_RUNS = {
    "cyp1a2": "20260928T113128124677Z-66cecb8c",
    "cyp2c9": "20260928T113128145767Z-366f1503",
    "cyp2d6": "20260928T113128150604Z-acf3b0b7",
    "cyp3a4": "20260928T113128171818Z-ca17bfa3",
}
ACCESSIONS = {
    "cyp1a2": "P05177",
    "cyp2c9": "P11712",
    "cyp2d6": "P10635",
    "cyp3a4": "P08684",
}
COMPOUND_RUN = "20260909T132734119233Z-a7063ccd"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("target", choices=sorted(REFERENCE_RUNS))
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs/shared"
    reference = runs / REFERENCE_RUNS[args.target] / "outputs"
    ensemble = runs / ENSEMBLE_RUNS[args.target] / "outputs"
    compounds = runs / COMPOUND_RUN / "outputs/compounds.parquet"
    aliases = runs / COMPOUND_RUN / "outputs/aliases.parquet"
    measurements = runs / COMPOUND_RUN / "outputs/measurements.parquet"
    measurement_sources = runs / COMPOUND_RUN / "outputs/measurement_sources.parquet"
    accession = ACCESSIONS[args.target]
    config = root / "configs/boltz-trunk" / (
        "prepare-smoke-v1.json" if args.smoke else "prepare-v1.json"
    )
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        f"Prepare template-free CYP/HEME/ligand trunk inputs for {args.target}",
        "--config",
        str(config),
        "--input",
        f"target.fasta={reference / f'{accession}.fasta'}",
        "--input",
        f"target.a3m={reference / f'{accession}.a3m'}",
        "--input",
        f"pocket_residues.json={ensemble / 'pocket_residues.json'}",
        "--input",
        f"compounds.parquet={compounds}",
        "--input",
        f"aliases.parquet={aliases}",
        "--input",
        f"measurements.parquet={measurements}",
        "--input",
        f"measurement_sources.parquet={measurement_sources}",
        "--",
        "python",
        "scripts/prepare_boltz_trunk_inputs.py",
    ]
    subprocess.run(command, cwd=root, check=True)


if __name__ == "__main__":
    main()
