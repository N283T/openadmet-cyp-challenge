"""Record HTTP source acquisition without mutating source files."""

import hashlib
import json
import os
from pathlib import Path

import httpx

out = Path(os.environ["CYP_OUTPUT_DIR"])
config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
records = []
for item in config["downloads"]:
    try:
        with httpx.Client(timeout=120, follow_redirects=True) as client:
            response = client.get(item["url"])
        path = out / item["filename"]
        path.write_bytes(response.content)
        records.append(
            {
                **item,
                "status": response.status_code,
                "bytes": len(response.content),
                "sha256": hashlib.sha256(response.content).hexdigest(),
                "resolved_url": str(response.url),
            }
        )
    except Exception as exc:
        records.append({**item, "error": str(exc)})
    print(records[-1], flush=True)
(out / "manifest.json").write_text(json.dumps(records, indent=2))
