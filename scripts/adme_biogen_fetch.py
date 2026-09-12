"""Acquire the primary-author Biogen public ADME CSV at a resolved Git commit."""

import hashlib
import json
import os
from pathlib import Path

import httpx

out = Path(os.environ["CYP_OUTPUT_DIR"])
with httpx.Client(timeout=60, follow_redirects=True) as client:
    url = "https://api.github.com/repos/molecularinformatics/Computational-ADME/commits/main"
    response = client.get(url)
    response.raise_for_status()
    (out / "github-commit.json").write_bytes(response.content)
    commit = response.json()["sha"]
    records = []
    for name in ["ADME_public_set_3521.csv", "README.md", "LICENSE"]:
        url = f"https://raw.githubusercontent.com/molecularinformatics/Computational-ADME/{commit}/{name}"
        r = client.get(url)
        r.raise_for_status()
        (out / name).write_bytes(r.content)
        records.append(
            dict(
                file=name,
                url=url,
                sha256=hashlib.sha256(r.content).hexdigest(),
                bytes=len(r.content),
            )
        )
(out / "source.json").write_text(
    json.dumps(
        dict(repository="molecularinformatics/Computational-ADME", commit=commit, files=records),
        indent=2,
    )
)
print(json.dumps(dict(commit=commit, files=records), indent=2))
