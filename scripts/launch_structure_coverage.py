"""Declare immutable inputs for the public CYP3A4 sequence-coverage analysis."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import polars as pl


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    config_path = root / "configs/structure/public-coverage.json"
    config = json.loads(config_path.read_text())
    run_root = root / "artifacts/runs/shared"
    source_outputs = run_root / config["source_run_id"] / "outputs"
    afdb_outputs = run_root / config["afdb_run_id"] / "outputs"
    alignment_outputs = run_root / config["alignment_run_id"] / "outputs"
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
        "Map resolved CYP3A4 PDB residues to the 503-residue AFDB reference sequence",
        "--config",
        str(config_path),
        "--input",
        f"source_manifest.csv={manifest_path}",
        "--input",
        f"alignment_manifest.csv={alignment_outputs / 'alignment_manifest.csv'}",
        "--input",
        f"afdb_model.cif={afdb_outputs / 'afdb_model.cif'}",
    ]
    for row in manifest.iter_rows(named=True):
        pdb_id = row["pdb_id"]
        command.extend(["--input", f"pdb_{pdb_id}.cif={source_outputs / row['file']}"])
    command.extend(["--", "python", "scripts/summarize_structure_coverage.py"])
    subprocess.run(command, check=True, cwd=root)


if __name__ == "__main__":
    main()
