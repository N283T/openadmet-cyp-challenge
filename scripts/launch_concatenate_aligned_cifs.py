"""Declare AFDB and all aligned CIF inputs for the plain concatenation run."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import polars as pl


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    config_path = root / "configs/structure/concatenate-aligned-cifs.json"
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
        "Concatenate AFDB reference and 122 aligned CYP3A4 mmCIF data blocks unchanged",
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
    command.extend(["--", "python", "scripts/concatenate_aligned_cifs.py"])
    subprocess.run(command, check=True, cwd=root)


if __name__ == "__main__":
    main()
