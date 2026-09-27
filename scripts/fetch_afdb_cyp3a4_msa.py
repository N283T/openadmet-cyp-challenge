"""Fetch the AFDB CYP3A4 MSA, or record the remote access failure."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import UTC, datetime
from pathlib import Path

import httpx


def main() -> None:
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    metadata = json.loads((input_dir / "afdb_prediction.json").read_text())[0]
    fasta_lines = (input_dir / "cyp3a4_sequence.fasta").read_text().splitlines()
    sequence = "".join(fasta_lines[1:])
    if (
        config["accession"] != metadata["uniprotAccession"]
        or config["entry_id"] != metadata["entryId"]
        or config["model_version"] != metadata["latestVersion"]
        or config["msa_url"] != metadata["msaUrl"]
        or sequence != metadata["sequence"]
    ):
        raise ValueError("AFDB MSA request does not match the saved full-length sequence")

    summary = {
        "checked_at_utc": datetime.now(UTC).isoformat(),
        "accession": config["accession"],
        "entry_id": config["entry_id"],
        "model_version": config["model_version"],
        "sequence_length": len(sequence),
        "requested_url": config["msa_url"],
    }
    try:
        with httpx.stream("GET", config["msa_url"], timeout=60, follow_redirects=True) as response:
            summary["resolved_url"] = str(response.url)
            summary["http_status"] = response.status_code
            if response.status_code != 200:
                summary["status"] = "unavailable"
                summary["response_excerpt"] = response.read()[:500].decode(
                    "utf-8", errors="replace"
                )
            else:
                data = response.read()
                records = data.decode("utf-8").splitlines()
                if len(records) < 2 or not records[0].startswith(">"):
                    raise ValueError("Downloaded MSA has no leading A3M query record")
                next_header = next(
                    (index for index in range(1, len(records)) if records[index].startswith(">")),
                    len(records),
                )
                query = "".join(records[1:next_header])
                if query.replace("-", "").upper() != sequence:
                    raise ValueError("Downloaded MSA query differs from full-length CYP3A4")
                n_sequences = sum(line.startswith(">") for line in records)
                if n_sequences < 2:
                    raise ValueError("Downloaded MSA has no homologous sequences")
                filename = "AF-P08684-F1-msa_v6.a3m"
                (output_dir / filename).write_bytes(data)
                summary.update(
                    status="available",
                    filename=filename,
                    bytes=len(data),
                    sha256=hashlib.sha256(data).hexdigest(),
                    n_sequences=n_sequences,
                )
    except httpx.HTTPError as exc:
        summary.update(status="unavailable", error=str(exc))
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
