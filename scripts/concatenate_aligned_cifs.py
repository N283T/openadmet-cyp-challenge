"""Concatenate aligned mmCIF data blocks without changing their contents."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import UTC, datetime
from pathlib import Path

import gemmi
import polars as pl


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    manifest = pl.read_csv(input_dir / "alignment_manifest.csv").sort("pdb_id")
    if manifest.height != config["expected_source_count"]:
        raise ValueError("Unexpected number of aligned CIFs")

    output_path = output_dir / "cyp3a4_aligned_122.cif"
    rows = []
    block_names = set()
    with output_path.open("wb") as combined:
        for source in manifest.iter_rows(named=True):
            pdb_id = source["pdb_id"]
            data = (input_dir / f"aligned_{pdb_id}.cif").read_bytes()
            if sha256(data) != source["aligned_sha256"]:
                raise ValueError(f"Aligned CIF hash mismatch for {pdb_id}")
            document = gemmi.cif.read_string(data.decode("utf-8"))
            if len(document) != 1:
                raise ValueError(f"Expected one CIF data block in {pdb_id}")
            block = document[0]
            if block.name.upper() != pdb_id.upper() or block.name in block_names:
                raise ValueError(f"Unexpected or duplicate data block: {block.name}")
            block_names.add(block.name)
            atoms = len(block.find_values("_atom_site.id"))
            if atoms != source["all_atom_sites"]:
                raise ValueError(f"Atom count mismatch for {pdb_id}")
            offset = combined.tell()
            combined.write(data)
            if not data.endswith(b"\n"):
                combined.write(b"\n")
            combined.write(b"\n")
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "data_block": block.name,
                    "source_sha256": source["aligned_sha256"],
                    "source_bytes": len(data),
                    "start_byte": offset,
                    "atoms": atoms,
                }
            )

    combined_document = gemmi.cif.read_file(str(output_path))
    if len(combined_document) != manifest.height:
        raise ValueError("Combined file has an unexpected number of data blocks")
    if [block.name for block in combined_document] != [row["data_block"] for row in rows]:
        raise ValueError("Combined file changed the data block order or identifiers")
    for block, row in zip(combined_document, rows, strict=True):
        if len(block.find_values("_atom_site.id")) != row["atoms"]:
            raise ValueError(f"Atom count changed in combined block {row['pdb_id']}")

    pl.DataFrame(rows).write_csv(output_dir / "block_map.csv")
    summary = {
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "alignment_run_id": config["alignment_run_id"],
        "data_blocks": len(combined_document),
        "atoms_across_blocks": sum(row["atoms"] for row in rows),
        "file": output_path.name,
        "bytes": output_path.stat().st_size,
        "sha256": sha256(output_path.read_bytes()),
        "method": "Byte-for-byte concatenation of source blocks, with newline separators only",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
