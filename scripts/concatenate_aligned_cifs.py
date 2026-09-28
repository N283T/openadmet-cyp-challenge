"""Concatenate AFDB and aligned mmCIF data blocks without changing contents."""

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


def append_block(
    combined, data: bytes, source_id: str, expected_sha256: str, expected_atoms: int,
    block_names: set[str],
) -> dict:
    if sha256(data) != expected_sha256:
        raise ValueError(f"CIF hash mismatch for {source_id}")
    document = gemmi.cif.read_string(data.decode("utf-8"))
    if len(document) != 1:
        raise ValueError(f"Expected one CIF data block in {source_id}")
    block = document[0]
    if block.name in block_names:
        raise ValueError(f"Duplicate data block: {block.name}")
    block_names.add(block.name)
    atoms = len(block.find_values("_atom_site.id"))
    if atoms != expected_atoms:
        raise ValueError(f"Atom count mismatch for {source_id}")
    offset = combined.tell()
    combined.write(data)
    if not data.endswith(b"\n"):
        combined.write(b"\n")
    combined.write(b"\n")
    return {
        "source_id": source_id,
        "data_block": block.name,
        "source_sha256": expected_sha256,
        "source_bytes": len(data),
        "start_byte": offset,
        "atoms": atoms,
    }


def main() -> None:
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    manifest = pl.read_csv(input_dir / "alignment_manifest.csv").sort("pdb_id")
    if manifest.height != config["expected_source_count"]:
        raise ValueError("Unexpected number of aligned CIFs")

    output_path = output_dir / "cyp3a4_afdb_aligned_123.cif"
    rows = []
    block_names = set()
    with output_path.open("wb") as combined:
        reference_data = (input_dir / "afdb_reference.cif").read_bytes()
        reference_atoms = len(
            gemmi.cif.read_string(reference_data.decode("utf-8"))[0].find_values("_atom_site.id")
        )
        if reference_atoms != config["reference_atoms"]:
            raise ValueError("Unexpected AFDB reference atom count")
        rows.append(
            append_block(
                combined, reference_data, "AFDB", config["reference_sha256"],
                reference_atoms, block_names,
            )
        )
        if rows[0]["data_block"] != "AF-P08684-F1":
            raise ValueError("Unexpected AFDB data block identifier")
        for source in manifest.iter_rows(named=True):
            pdb_id = source["pdb_id"]
            data = (input_dir / f"aligned_{pdb_id}.cif").read_bytes()
            row = append_block(
                combined, data, pdb_id, source["aligned_sha256"],
                source["all_atom_sites"], block_names,
            )
            if row["data_block"].upper() != pdb_id.upper():
                raise ValueError(f"Unexpected data block identifier for {pdb_id}")
            rows.append(row)

    combined_document = gemmi.cif.read_file(str(output_path))
    if len(combined_document) != manifest.height + 1:
        raise ValueError("Combined file has an unexpected number of data blocks")
    if [block.name for block in combined_document] != [row["data_block"] for row in rows]:
        raise ValueError("Combined file changed the data block order or identifiers")
    for block, row in zip(combined_document, rows, strict=True):
        if len(block.find_values("_atom_site.id")) != row["atoms"]:
            raise ValueError(f"Atom count changed in combined block {row['source_id']}")

    pl.DataFrame(rows).write_csv(output_dir / "block_map.csv")
    summary = {
        "completed_at_utc": datetime.now(UTC).isoformat(),
        "alignment_run_id": config["alignment_run_id"],
        "afdb_reference_block": rows[0]["data_block"],
        "afdb_reference_atoms": rows[0]["atoms"],
        "aligned_source_blocks": manifest.height,
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
