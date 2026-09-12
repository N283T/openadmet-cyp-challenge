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


inp = Path(os.environ["CYP_INPUT_DIR"])


def source(name):
    paths = [p for p in inp.iterdir() if p.name.partition("-")[2] == name]
    assert len(paths) == 1
    return json.loads(paths[0].read_text())


commits = source("cyp-commits.txt")
summary = dict(
    cyp_commits=[
        dict(sha=c["sha"], message=c["commit"]["message"], date=c["commit"]["author"]["date"])
        for c in commits
    ]
)
for repo, branch in [
    ("jeremycheminf.github.io", "main"),
    ("colab_notebooks", "main"),
    ("streamlit_apps", "main"),
    ("chemprop", "master"),
]:
    tree = fetch(
        f"https://api.github.com/repos/jeremycheminf/{repo}/git/trees/{branch}?recursive=1",
        repo + "-tree",
    )
    if tree:
        assert not tree.get("truncated")
        paths = [t["path"] for t in tree["tree"]]
        summary[repo] = paths
        candidates = [
            t
            for t in tree["tree"]
            if t["type"] == "blob"
            and t.get("size", 0) < 200000
            and re.search(r"adme|pretrain|cyp|chemprop|README", t["path"], re.I)
            and t["path"].endswith((".md", ".py", ".json", ".yaml", ".toml"))
        ][:12]
        if repo == "chemprop":
            candidates = [t for t in candidates if t["path"] == "README.md"]
        for i, item in enumerate(candidates):
            fetch(
                f"https://raw.githubusercontent.com/jeremycheminf/{repo}/{branch}/{item['path']}",
                repo + f"-source-{i}",
            )
for i, commit in enumerate(commits[:10]):
    sha = commit["sha"]
    fetch(
        f"https://raw.githubusercontent.com/jeremycheminf/openadmet_scripts/{sha}/CYP_Challenge/README.md",
        f"history-readme-{i}",
    )
    tree = fetch(
        f"https://api.github.com/repos/jeremycheminf/openadmet_scripts/git/trees/{sha}?recursive=1",
        f"history-tree-{i}",
    )
    if tree:
        summary[f"history-paths-{i}"] = [t["path"] for t in tree["tree"]]
summary["snapshots"] = records
(out / "detail-summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
