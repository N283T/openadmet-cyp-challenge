"""Bounded public author repository/history inventory with immutable HTTP snapshots."""

import hashlib
import json
import os
import re
import urllib.error
import urllib.request
from pathlib import Path

out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
records = []


def fetch(url, label):
    request = urllib.request.Request(url, headers={"User-Agent": "CYP-public-source-inventory"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            raw, status = response.read(5_000_001), response.status
        assert len(raw) <= 5_000_000
    except urllib.error.HTTPError as error:
        raw, status = error.read(), error.code
    except (TimeoutError, urllib.error.URLError) as error:
        raw, status = str(error).encode(), 0
    filename = label + ".txt"
    (out / filename).write_bytes(raw)
    records.append(
        dict(url=url, file=filename, status=status, sha256=hashlib.sha256(raw).hexdigest())
    )
    try:
        return json.loads(raw) if status == 200 else None
    except json.JSONDecodeError:
        return None


repos = []
for page in range(1, 4):
    result = fetch(
        f"https://api.github.com/users/{cfg['author']}/repos?per_page=100&page={page}",
        f"repositories-{page}",
    )
    if not isinstance(result, list):
        break
    repos.extend(result)
    if len(result) < 100:
        break
selected = [
    r
    for r in repos
    if re.search(cfg["pattern"], r["name"] + " " + str(r.get("description")), re.I)
    and not r["fork"]
][:12]
for name in cfg["probe_names"]:
    fetch(f"https://api.github.com/repos/{cfg['author']}/{name}", "probe-" + name)
for repo in selected:
    name, branch = repo["name"], repo["default_branch"]
    fetch(
        f"https://api.github.com/repos/{cfg['author']}/{name}/git/trees/{branch}?recursive=1",
        name + "-tree",
    )
    fetch(
        f"https://raw.githubusercontent.com/{cfg['author']}/{name}/{branch}/README.md",
        name + "-readme",
    )
fetch(
    f"https://api.github.com/repos/{cfg['author']}/openadmet_scripts/commits?"
    "path=CYP_Challenge&per_page=100",
    "cyp-commits",
)
fetch(
    f"https://api.github.com/repos/{cfg['author']}/openadmet_scripts/branches?per_page=100",
    "openadmet-branches",
)
summary = dict(
    config=cfg,
    repositories=[
        {k: r.get(k) for k in ["name", "description", "html_url", "fork", "default_branch"]}
        for r in repos
    ],
    selected=[r["name"] for r in selected],
    snapshots=records,
)
(out / "inventory.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
