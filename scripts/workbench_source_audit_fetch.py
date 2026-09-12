"""Pin and snapshot public Workbench CYP recipes and Chemprop implementation."""

import json
import os
from pathlib import Path

import httpx

out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
client = httpx.Client(timeout=60, follow_redirects=True)


def get(url):
    r = client.get(url)
    r.raise_for_status()
    return r


base = f"{cfg['api']}/repos/{cfg['repository']}"
commit = get(f"{base}/commits/{cfg['ref']}").json()
sha = commit["sha"]
(out / "commit.json").write_text(json.dumps(commit, indent=2))
tree = get(f"{base}/git/trees/{sha}?recursive=1").json()
(out / "tree.json").write_text(json.dumps(tree, indent=2))
assert not tree.get("truncated")
paths = [
    r["path"]
    for r in tree["tree"]
    if r["type"] == "blob"
    and (
        (
            r["path"].startswith("ml_pipelines/OpenADMET/cyp/")
            and r["path"].endswith((".py", ".md", ".json", ".yaml"))
        )
        or ("chemprop" in r["path"].lower() and r["path"].endswith((".py", ".json", ".yaml")))
        or r["path"] in ["pyproject.toml", "requirements.txt"]
    )
]
for p in paths:
    dest = out / "source" / p
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(get(f"{cfg['raw']}/{cfg['repository']}/{sha}/{p}").content)
for p in cfg["history_paths"]:
    h = get(f"{base}/commits?path={p}&per_page=30").json()
    dest = out / "history" / f"{Path(p).stem}.json"
    dest.parent.mkdir(exist_ok=True)
    dest.write_text(json.dumps(h, indent=2))
    # Every recent recipe version is small; capture contents to distinguish 18 vs26head recipe.
    for entry in h:
        oldsha = entry["sha"]
        r = client.get(f"{cfg['raw']}/{cfg['repository']}/{oldsha}/{p}")
        if r.status_code == 200:
            (out / "history" / f"{Path(p).stem}-{oldsha}.py").write_bytes(r.content)
(out / "summary.json").write_text(
    json.dumps(dict(repository=cfg["repository"], commit=sha, paths=paths), indent=2)
)
print(sha, flush=True)
