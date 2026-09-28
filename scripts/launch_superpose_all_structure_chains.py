"""Declare every original CYP3A4 complex for independent chain superposition."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import polars as pl


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    config_path = root / "configs/structure/all-chains-to-afdb.json"
    config = json.loads(config_path.read_text())
    runs = root / "artifacts/runs/shared"
    source_outputs = runs / config["source_run_id"] / "outputs"
    afdb_outputs = runs / config["afdb_run_id"] / "outputs"
    manifest_path = source_outputs / "source_manifest.csv"
    manifest = pl.read_csv(manifest_path)
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        "Superpose every public CYP3A4 chain with its HEM and ligands onto AFDB",
        "--config",
        str(config_path),
        "--input",
        f"source_manifest.csv={manifest_path}",
        "--input",
        f"afdb_model.cif={afdb_outputs / 'afdb_model.cif'}",
    ]
    for row in manifest.iter_rows(named=True):
        pdb_id = row["pdb_id"]
        command.extend(["--input", f"pdb_{pdb_id}.cif={source_outputs / row['file']}"])
    command.extend(["--", "python", "scripts/superpose_all_structure_chains.py"])
    subprocess.run(command, check=True, cwd=root)


if __name__ == "__main__":
    main()
