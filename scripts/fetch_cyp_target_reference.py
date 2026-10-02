"""Snapshot one CYP AFDB reference and its pinned LiteFold MSA."""

from __future__ import annotations

import hashlib
import io
import json
import os
from datetime import UTC, datetime
from pathlib import Path

import gemmi
import httpx
import polars as pl


def fetch(url: str) -> tuple[bytes, str]:
    response = httpx.get(url, timeout=120, follow_redirects=True)
    response.raise_for_status()
    return response.content, str(response.url)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    base = (
        f"https://huggingface.co/datasets/{config['msa_dataset']}"
        f"/resolve/{config['msa_revision']}"
    )
    manifest_data, manifest_url = fetch(f"{base}/manifest.tsv")
    manifest = pl.read_csv(io.BytesIO(manifest_data), separator="\t")
    matches = manifest.filter(pl.col("acc") == config["accession"])
    if matches.height != 1:
        raise ValueError(f"Expected one LiteFold row for {config['accession']}")
    row = matches.row(0, named=True)
    sequence = row["sequence"]
    if (
        row["length"] != config["expected_length"]
        or len(sequence) != config["expected_length"]
        or row["a3m_path"] != config["msa_path"]
        or row["a3m_bytes"] != config["msa_bytes"]
    ):
        raise ValueError("Pinned LiteFold manifest row differs from target config")

    a3m_data, a3m_url = fetch(f"{base}/{config['msa_path']}")
    if len(a3m_data) != config["msa_bytes"]:
        raise ValueError("MSA byte count differs from pinned manifest")
    lines = a3m_data.decode().splitlines()
    header_indices = [index for index, line in enumerate(lines) if line.startswith(">")]
    if not header_indices:
        raise ValueError("MSA contains no FASTA records")
    query_end = header_indices[1] if len(header_indices) > 1 else len(lines)
    query = "".join(lines[header_indices[0] + 1 : query_end]).replace("-", "").upper()
    if query != sequence:
        raise ValueError("MSA query does not match the LiteFold target sequence")
    core_lengths = {
        sum(character.isupper() or character == "-" for character in line)
        for line in lines
        if line and not line.startswith(">")
    }
    if core_lengths != {len(sequence)}:
        raise ValueError("MSA rows do not span the full query alignment")

    metadata_data, metadata_url = fetch(config["afdb_api_url"])
    predictions = json.loads(metadata_data)
    canonical = [
        item
        for item in predictions
        if item.get("entryId") == config["afdb_entry_id"]
        and item.get("uniprotAccession") == config["accession"]
    ]
    if len(canonical) != 1:
        raise ValueError("Expected exactly one matching canonical AFDB prediction")
    prediction = canonical[0]
    if prediction.get("sequence") != sequence:
        raise ValueError("AFDB metadata differs from the pinned LiteFold sequence")
    cif_data, cif_url = fetch(prediction["cifUrl"])
    cif_path = output_dir / "afdb_model.cif"
    cif_path.write_bytes(cif_data)
    structure = gemmi.read_structure(str(cif_path))
    if len(structure) != 1 or len(structure[0]) != 1:
        raise ValueError("AFDB reference is not a single-chain model")
    polymer = structure[0][0].get_polymer()
    observed = gemmi.one_letter_code([residue.name for residue in polymer])
    if len(polymer) != len(sequence) or observed != sequence:
        raise ValueError("AFDB coordinates do not contain the full target sequence")

    fasta_name = f"{config['accession']}.fasta"
    msa_name = f"{config['accession']}.a3m"
    (output_dir / fasta_name).write_text(
        f">sp|{config['accession']}|{config['target']} canonical\n{sequence}\n"
    )
    (output_dir / msa_name).write_bytes(a3m_data)
    (output_dir / "afdb_prediction.json").write_bytes(metadata_data)
    summary = {
        "collected_at_utc": datetime.now(UTC).isoformat(),
        "target": config["target"],
        "accession": config["accession"],
        "sequence_length": len(sequence),
        "afdb_entry_id": prediction["entryId"],
        "afdb_model_created_date": prediction.get("modelCreatedDate"),
        "afdb_metadata_url": metadata_url,
        "afdb_metadata_sha256": sha256(metadata_data),
        "afdb_cif_url": cif_url,
        "afdb_cif_bytes": len(cif_data),
        "afdb_cif_sha256": sha256(cif_data),
        "msa_dataset": config["msa_dataset"],
        "msa_revision": config["msa_revision"],
        "msa_manifest_url": manifest_url,
        "msa_manifest_sha256": sha256(manifest_data),
        "msa_url": a3m_url,
        "msa_bytes": len(a3m_data),
        "msa_sha256": sha256(a3m_data),
        "msa_sequences": len(header_indices),
        "msa_aligned_columns": len(sequence),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
