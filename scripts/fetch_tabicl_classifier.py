"""Pin the classifier checkpoint separately from previously used regression weights."""
import json
import os
import shutil
from pathlib import Path

from huggingface_hub import HfApi, hf_hub_download

from cyp_challenge.storage import digest

out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
revision = HfApi().model_info(cfg["repo"]).sha
path = hf_hub_download(cfg["repo"], cfg["filename"], revision=revision)
shutil.copy2(path, out / cfg["filename"])
(out / "summary.json").write_text(json.dumps(dict(
    **cfg, revision=revision, sha256=digest(out / cfg["filename"])), indent=2))
print(revision, cfg["filename"], flush=True)
