"""Record all raw tensor inputs and pool each target without duplicating 840 GB."""

import argparse
import json
import subprocess
from pathlib import Path

import polars as pl


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("target", choices=["cyp1a2", "cyp2c9", "cyp2d6", "cyp3a4"])
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    target = args.target
    canonical = root / "artifacts/boltz-trunk/post-r1"
    manifest = canonical / target / "manifest.parquet"
    definition = root / "configs/structure/targets" / f"{target}.json"
    length = json.loads(definition.read_text())["expected_length"]
    config = root / "configs/boltz-trunk" / f"pool-{target}-v1.json"
    if not config.exists():
        config.write_text(
            json.dumps(
                {
                    "target": target,
                    "protein_length": length,
                    "workers": 8,
                    "representation": "post-r1",
                    "input_storage": "hardlinks to retained frozen canonical tensors",
                },
                indent=2,
            )
            + "\n"
        )
    boltz = Path("/home/nagaet/ghq/github.com/N283T/boltz-trunk")
    files = {
        "manifest.parquet": manifest,
        "pocket.json": root
        / "artifacts/tracked/cyp-structure-prep-20260928"
        / target
        / "pocket_residues.json",
        "HEM.pkl": Path.home() / ".boltz/mols/HEM.pkl",
        "boltz_schema.py": boltz / "src/boltz/data/parse/schema.py",
        "boltz_tokenizer.py": boltz / "src/boltz/data/tokenize/boltz2.py",
    }
    for row in pl.read_parquet(manifest).iter_rows(named=True):
        files[f"record_{row['record_index']:06d}.npz"] = canonical / row["path"]
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        f"Pool all {target} post-r1 CYP/HEM tensors for PXR-style region and HEM ablations",
        "--config",
        str(config),
        "--input-mode",
        "hardlink",
    ]
    for alias, path in files.items():
        command += ["--input", f"{alias}={path}"]
    subprocess.run(
        command + ["--", "python", "scripts/pool_cyp_boltz_trunk.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
