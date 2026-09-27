"""Snapshot and validate a full-length CYP3A4 MSA from LiteFold."""

from __future__ import annotations

import hashlib
import io
import json
import os
from datetime import UTC, datetime
from pathlib import Path

import httpx
import polars as pl


def fetch(url: str) -> bytes:
    response = httpx.get(url, timeout=120, follow_redirects=True)
    response.raise_for_status()
    return response.content


def main() -> None:
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    fasta_lines = (input_dir / "cyp3a4_sequence.fasta").read_text().splitlines()
    sequence = "".join(fasta_lines[1:])
    if config["accession"] not in fasta_lines[0] or len(sequence) != 503:
        raise ValueError("Expected the official full-length P08684 sequence")

    base = f"https://huggingface.co/datasets/{config['dataset']}/resolve/{config['revision']}"
    manifest_data = fetch(f"{base}/manifest.tsv")
    manifest = pl.read_csv(io.BytesIO(manifest_data), separator="\t")
    matches = manifest.filter(pl.col("acc") == config["accession"])
    if matches.height != 1:
        raise ValueError(f"Expected one manifest row for {config['accession']}")
    row = matches.row(0, named=True)
    if (
        row["sequence"] != sequence
        or row["length"] != len(sequence)
        or row["a3m_path"] != config["expected_a3m_path"]
        or row["a3m_bytes"] != config["expected_a3m_bytes"]
    ):
        raise ValueError("Manifest P08684 row disagrees with the official sequence/config")

    a3m_data = fetch(f"{base}/{row['a3m_path']}")
    if len(a3m_data) != row["a3m_bytes"]:
        raise ValueError("A3M size differs from the pinned manifest")
    lines = a3m_data.decode("utf-8").splitlines()
    first_header = next((i for i, line in enumerate(lines) if line.startswith(">")), None)
    if first_header is None:
        raise ValueError("A3M has no FASTA headers")
    next_header = next(
        (i for i in range(first_header + 1, len(lines)) if lines[i].startswith(">")),
        len(lines),
    )
    query = "".join(lines[first_header + 1 : next_header]).replace("-", "").upper()
    if query != sequence:
        raise ValueError("A3M query differs from full-length CYP3A4")
    n_sequences = sum(line.startswith(">") for line in lines)
    if n_sequences < 2:
        raise ValueError("A3M contains no homologous sequences")
    core_lengths = [
        sum(char.isupper() or char == "-" for char in line)
        for line in lines
        if line and not line.startswith(">")
    ]
    if len(core_lengths) != n_sequences or set(core_lengths) != {len(sequence)}:
        raise ValueError("A3M aligned rows do not all span the 503-residue query")

    filename = "P08684.a3m"
    (output_dir / filename).write_bytes(a3m_data)
    manifest_row = {key: row[key] for key in row}
    (output_dir / "manifest_row.json").write_text(json.dumps(manifest_row, indent=2) + "\n")
    summary = {
        "collected_at_utc": datetime.now(UTC).isoformat(),
        "dataset": config["dataset"],
        "revision": config["revision"],
        "generation": "ColabFold/MMseqs2 GPU; UniRef30 2023-02 and ColabFold environmental 2021-08",
        "license": "cc-by-4.0",
        "accession": config["accession"],
        "sequence_length": len(sequence),
        "n_sequences": n_sequences,
        "aligned_columns": len(sequence),
        "query_header": lines[first_header],
        "manifest_url": f"{base}/manifest.tsv",
        "manifest_sha256": hashlib.sha256(manifest_data).hexdigest(),
        "a3m_url": f"{base}/{row['a3m_path']}",
        "a3m_bytes": len(a3m_data),
        "a3m_sha256": hashlib.sha256(a3m_data).hexdigest(),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
