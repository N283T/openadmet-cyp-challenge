"""Declare every immutable structure input for the AFDB alignment run."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import polars as pl


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    config_path = root / "configs/structure/align-public-to-afdb.json"
    config = json.loads(config_path.read_text())
    source_outputs = root / "artifacts/runs/shared" / config["source_run_id"] / "outputs"
    afdb_outputs = root / "artifacts/runs/shared" / config["afdb_run_id"] / "outputs"
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
        "Superpose public CYP3A4 heme structures onto AFDB P08684 coordinates",
        "--config",
        str(config_path),
        "--input",
        f"source_manifest.csv={manifest_path}",
        "--input",
        f"afdb_model.cif={afdb_outputs / 'afdb_model.cif'}",
    ]
    for row in manifest.iter_rows(named=True):
        pdb_id = row["pdb_id"]
        command.extend(
            ["--input", f"pdb_{pdb_id}.cif={source_outputs / row['file']}"]
        )
    command.extend(["--", "python", "scripts/align_structure_sources.py"])
    subprocess.run(command, check=True, cwd=root)


if __name__ == "__main__":
    main()
