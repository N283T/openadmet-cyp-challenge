"""Snapshot declared public GET sources, including unsuccessful response attempts."""

import json
import os
import time
from pathlib import Path

import httpx

from cyp_challenge.storage import digest, now, write_json


def main():
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    records = []
    with httpx.Client(timeout=180, follow_redirects=True) as client:
        for item in cfg["requests"]:
            for attempt in range(3):
                path = out / f"{item['name']}.attempt{attempt}"
                started = now()
                response = client.get(item["url"])
                path.write_bytes(response.content)
                record = dict(
                    name=item["name"],
                    url=item["url"],
                    final_url=str(response.url),
                    retrieved_at=started,
                    status=response.status_code,
                    sha256=digest(path),
                    bytes=path.stat().st_size,
                    headers=dict(response.headers),
                    attempt=attempt,
                )
                records.append(record)
                write_json(out / f"{item['name']}.attempt{attempt}.metadata.json", record)
                print(
                    json.dumps({k: record[k] for k in ["name", "status", "bytes", "attempt"]}),
                    flush=True,
                )
                if response.is_success:
                    (out / item["name"]).write_bytes(response.content)
                    break
                if attempt == 2 or response.status_code not in [429, 500, 502, 503, 504]:
                    response.raise_for_status()
                time.sleep(3 * (attempt + 1))
    write_json(out / "requests.json", records)


if __name__ == "__main__":
    main()
