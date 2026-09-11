"""Archive revision-pinned CYP reproduction source without executing it."""

import hashlib
import json
import os
from pathlib import Path

import httpx

out = Path(os.environ["CYP_OUTPUT_DIR"])
repo = "jeremycheminf/openadmet_scripts"
with httpx.Client(follow_redirects=True, timeout=90) as client:
    r = client.get(f"https://api.github.com/repos/{repo}/commits/main")
    r.raise_for_status()
    commit = r.json()
    revision = commit["sha"]
    r = client.get(f"https://api.github.com/repos/{repo}/git/trees/{revision}?recursive=1")
    r.raise_for_status()
    tree = r.json()
    (out / "tree.json").write_text(json.dumps(tree, indent=2))
    files = []
    for item in tree["tree"]:
        name = item["path"]
        if item["type"] != "blob" or not name.startswith("CYP_Challenge/"):
            continue
        if (
            not name.endswith((".py", ".md", ".toml", ".json", ".txt"))
            or item.get("size", 0) > 2000000
        ):
            continue
        url = f"https://raw.githubusercontent.com/{repo}/{revision}/{name}"
        r = client.get(url)
        r.raise_for_status()
        target = out / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(r.content)
        files.append(dict(path=name, url=url, sha256=hashlib.sha256(r.content).hexdigest()))
    (out / "provenance.json").write_text(
        json.dumps(
            dict(
                repository=repo,
                revision=revision,
                commit_date=commit["commit"]["committer"]["date"],
                files=files,
            ),
            indent=2,
        )
    )
    print(revision, len(files))
