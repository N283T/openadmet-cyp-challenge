"""Record immutable inputs for one CIF with one fitted chain per data block."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    config_path = root / "configs/structure/multiblock-superposed.json"
    config = json.loads(config_path.read_text())
    outputs = root / "artifacts/runs/shared" / config["all_chain_run_id"] / "outputs"
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        "Store AFDB and 142 fitted CYP3A4 chains in separate data blocks of one CIF",
        "--config",
        str(config_path),
        "--input",
        f"all_chains.cif={outputs / 'cyp3a4_all_142_chains_on_afdb.cif'}",
        "--input",
        f"all_chains_summary.json={outputs / 'summary.json'}",
        "--input",
        f"chain_alignment.csv={outputs / 'chain_alignment.csv'}",
        "--",
        "python",
        "scripts/build_multiblock_superposed_cif.py",
    ]
    subprocess.run(command, check=True, cwd=root)


if __name__ == "__main__":
    main()
