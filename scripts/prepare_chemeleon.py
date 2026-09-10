"""Verify existing PXR checkpoint against Zenodo and snapshot official reference code."""

import hashlib
import json
import os
import shutil
from pathlib import Path

import httpx

from cyp_challenge.storage import digest, now, write_json


def main():
    inputs, output = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    checkpoint = next(inputs.glob("*-chemeleon_mp.pt"))
    with httpx.Client(timeout=90, follow_redirects=True) as client:
        verification = {}
        try:
            if config.get("reuse_local_without_remote_verification", False):
                raise httpx.ConnectError(
                    "Remote verification deferred after recorded Zenodo outages"
                )
            response = client.get(config["record_api"])
            (output / "zenodo-api-response.txt").write_bytes(response.content)
            response.raise_for_status()
            metadata = response.json()
            write_json(output / "zenodo-record.json", metadata)
            entry = next(f for f in metadata["files"] if f["key"] == "chemeleon_mp.pt")
            algorithm, expected = entry["checksum"].split(":", 1)
            with checkpoint.open("rb") as stream:
                actual = hashlib.file_digest(stream, algorithm).hexdigest()
            if actual != expected or checkpoint.stat().st_size != entry["size"]:
                raise ValueError("Existing checkpoint does not match official Zenodo file")
            verification = {"method": "record checksum", "checksum": entry["checksum"]}
        except httpx.HTTPError as exc:
            write_json(output / "api-failure.json", {"error": str(exc)})
            if config.get("reuse_local_without_remote_verification", False):
                verification = {
                    "method": "local SHA256 pin only",
                    "official_checksum_verified": False,
                    "reason": (
                        "Zenodo API and direct download returned HTTP 504 in prior recorded runs"
                    ),
                }
            else:
                downloaded = output / "official-download.pt"
                with client.stream("GET", config["checkpoint_url"] + "?download=1") as response:
                    response.raise_for_status()
                    with downloaded.open("xb") as stream:
                        for chunk in response.iter_bytes():
                            stream.write(chunk)
                if digest(downloaded) != digest(checkpoint):
                    raise ValueError(
                        "Existing checkpoint differs from official file download"
                    ) from exc
                verification = {"method": "SHA256 equality with fresh official file download"}
        revision = config["upstream_revision"]
        for filename in ["chemeleon_fingerprint.py", "LICENSE", "README.md"]:
            url = f"https://raw.githubusercontent.com/JacksonBurns/chemeleon/{revision}/{filename}"
            reference = client.get(url)
            reference.raise_for_status()
            (output / filename).write_bytes(reference.content)
    shutil.copy2(checkpoint, output / "chemeleon_mp.pt")
    write_json(
        output / "provenance.json",
        {
            "retrieved_at": now(),
            "config": config,
            "verification": verification,
            "checkpoint_sha256": digest(checkpoint),
            "bytes": checkpoint.stat().st_size,
            "origin": "Existing PXR checkpoint; see verification for actual verification status",
        },
    )
    print("Prepared", digest(checkpoint), flush=True)


if __name__ == "__main__":
    main()
