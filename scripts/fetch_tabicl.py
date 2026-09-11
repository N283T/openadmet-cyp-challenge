"""Snapshot the fixed upstream TabICL regression checkpoint."""

import os
import shutil
from pathlib import Path

from huggingface_hub import hf_hub_download

from cyp_challenge.storage import digest, write_json

out = Path(os.environ["CYP_OUTPUT_DIR"])
name = "tabicl-regressor-v2-20260212.ckpt"
p = hf_hub_download("jingang/TabICL", name)
shutil.copy2(p, out / name)
write_json(
    out / "summary.json",
    {
        "repo": "jingang/TabICL",
        "filename": name,
        "resolved_cache_path": p,
        "sha256": digest(out / name),
    },
)
