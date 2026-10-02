"""Launch a recorded protein-only ensemble and pocket-definition build."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import polars as pl


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("target", choices=["cyp1a2", "cyp2c9", "cyp2d6", "cyp3a4"])
    parser.add_argument("reference_run")
    parser.add_argument("structure_run")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs/shared"
    reference = runs / args.reference_run / "outputs"
    structures = runs / args.structure_run / "outputs"
    target_config = root / f"configs/structure/targets/{args.target}.json"
    config = root / "configs/structure/ensemble-complex.json"
    accession = {
        "cyp1a2": "P05177",
        "cyp2c9": "P11712",
        "cyp2d6": "P10635",
        "cyp3a4": "P08684",
    }[args.target]
    manifest = pl.read_csv(structures / "source_manifest.csv")
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        f"Build AFDB-aligned CYP/HEME/ligand ensemble and contact pocket for {args.target}",
        "--config",
        str(config),
        "--input",
        f"target_config.json={target_config}",
        "--input",
        f"target.fasta={reference / f'{accession}.fasta'}",
        "--input",
        f"afdb_model.cif={reference / 'afdb_model.cif'}",
        "--input",
        f"source_manifest.csv={structures / 'source_manifest.csv'}",
    ]
    for pdb_id in manifest["pdb_id"]:
        command.extend(
            ["--input", f"pdb_{pdb_id}.cif={structures / 'pdb' / f'{pdb_id}.cif'}"]
        )
    command.extend(["--", "python", "scripts/prepare_cyp_structure_ensemble.py"])
    subprocess.run(command, cwd=root, check=True)


if __name__ == "__main__":
    main()
