"""Declare all aligned public CIF inputs for the combined overlay run."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import polars as pl


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    config_path = root / "configs/structure/combined-public-overlay.json"
    config = json.loads(config_path.read_text())
    outputs = root / "artifacts/runs/shared" / config["alignment_run_id"] / "outputs"
    manifest_path = outputs / "alignment_manifest.csv"
    manifest = pl.read_csv(manifest_path)
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        "Combine 122 AFDB-superposed public CYP3A4/HEM complexes into one mmCIF",
        "--config",
        str(config_path),
        "--input",
        f"alignment_manifest.csv={manifest_path}",
        "--input",
        f"afdb_reference.cif={outputs / 'afdb_reference.cif'}",
    ]
    for row in manifest.iter_rows(named=True):
        pdb_id = row["pdb_id"]
        command.extend(
            ["--input", f"aligned_{pdb_id}.cif={outputs / row['aligned_file']}"]
        )
    command.extend(["--", "python", "scripts/build_combined_structure_overlay.py"])
    subprocess.run(command, check=True, cwd=root)


if __name__ == "__main__":
    main()
