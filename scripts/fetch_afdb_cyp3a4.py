"""Snapshot the official AFDB CYP3A4 prediction and its metadata."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import UTC, datetime
from pathlib import Path

import gemmi
import httpx


def fetch(url: str) -> tuple[bytes, str]:
    response = httpx.get(url, timeout=45, follow_redirects=True)
    response.raise_for_status()
    return response.content, str(response.url)


def main() -> None:
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    fasta_lines = (input_dir / "cyp3a4_sequence.fasta").read_text().splitlines()
    sequence = "".join(fasta_lines[1:])
    if config["accession"] not in fasta_lines[0]:
        raise ValueError("FASTA accession does not match the AFDB request")

    metadata_bytes, metadata_url = fetch(config["api_url"])
    predictions = json.loads(metadata_bytes)
    if len(predictions) != 1:
        raise ValueError(f"Expected one AFDB prediction, got {len(predictions)}")
    prediction = predictions[0]
    if (
        prediction.get("entryId") != config["expected_entry_id"]
        or prediction.get("uniprotAccession") != config["accession"]
        or prediction.get("sequence") != sequence
    ):
        raise ValueError("AFDB prediction does not match the official CYP3A4 sequence")

    files = {}
    for key, filename in (("cifUrl", "afdb_model.cif"), ("pdbUrl", "afdb_model.pdb")):
        data, resolved_url = fetch(prediction[key])
        if len(data) < 100_000:
            raise ValueError(f"Unexpectedly small AFDB file: {filename}")
        path = output_dir / filename
        path.write_bytes(data)
        files[filename] = {
            "url": resolved_url,
            "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
        }
    cif_structure = gemmi.read_structure(str(output_dir / "afdb_model.cif"))
    pdb_structure = gemmi.read_structure(str(output_dir / "afdb_model.pdb"))
    cif_ca = sum(
        atom.name == "CA" for chain in cif_structure[0] for residue in chain for atom in residue
    )
    pdb_ca = sum(
        atom.name == "CA" for chain in pdb_structure[0] for residue in chain for atom in residue
    )
    if cif_ca != len(sequence) or pdb_ca != len(sequence):
        raise ValueError(f"AFDB model CA count mismatch: CIF={cif_ca}, PDB={pdb_ca}")

    (output_dir / "afdb_prediction.json").write_bytes(metadata_bytes)
    summary = {
        "collected_at_utc": datetime.now(UTC).isoformat(),
        "entry_id": prediction["entryId"],
        "accession": prediction["uniprotAccession"],
        "model_created_date": prediction.get("modelCreatedDate"),
        "api_url": metadata_url,
        "metadata_sha256": hashlib.sha256(metadata_bytes).hexdigest(),
        "sequence_length": len(sequence),
        "ca_count": cif_ca,
        "files": files,
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
