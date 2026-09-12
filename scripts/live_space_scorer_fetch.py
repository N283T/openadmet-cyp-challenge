"""Pin public live Space source and tutorial history without requesting labels."""

import json
import os
from pathlib import Path

import httpx

out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
with httpx.Client(timeout=60, follow_redirects=True) as client:

    def get(url):
        r = client.get(url)
        r.raise_for_status()
        return r

    info = get(f"https://huggingface.co/api/spaces/{cfg['space']}").json()
    sha = info["sha"]
    (out / "space-info.json").write_text(json.dumps(info, indent=2))
    tree = get(
        f"https://huggingface.co/api/spaces/{cfg['space']}/tree/{sha}?recursive=true&expand=false"
    ).json()
    (out / "space-tree.json").write_text(json.dumps(tree, indent=2))
    for e in tree:
        p = e.get("path", "")
        if (
            e.get("type") == "file"
            and p.endswith((".py", ".txt", ".toml", ".yaml", ".yml", ".md"))
            and not p.startswith(("data/", "uploads/"))
        ):
            r = get(f"https://huggingface.co/spaces/{cfg['space']}/resolve/{sha}/{p}")
            dest = out / "space" / p
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(r.content)
    history = get(f"https://huggingface.co/api/spaces/{cfg['space']}/commits/main").json()
    (out / "space-commits.json").write_text(json.dumps(history, indent=2))
    gh = f"https://api.github.com/repos/{cfg['tutorial']}"
    ti = get(gh + "/commits/main").json()
    (out / "tutorial-commit.json").write_text(json.dumps(ti, indent=2))
    tt = get(gh + f"/git/trees/{ti['sha']}?recursive=1").json()
    (out / "tutorial-tree.json").write_text(json.dumps(tt, indent=2))
    for e in tt["tree"]:
        p = e["path"]
        if (
            e["type"] == "blob"
            and ("scor" in p.lower() or p.endswith("requirements.txt"))
            and p.endswith((".py", ".txt"))
        ):
            dest = out / "tutorial" / p
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(
                get(f"https://raw.githubusercontent.com/{cfg['tutorial']}/{ti['sha']}/{p}").content
            )
            hh = get(gh + f"/commits?path={p}&per_page=30").json()
            (out / f"history-{Path(p).name}.json").write_text(json.dumps(hh, indent=2))
    print(json.dumps(dict(space_commit=sha, tutorial_commit=ti["sha"])), flush=True)
