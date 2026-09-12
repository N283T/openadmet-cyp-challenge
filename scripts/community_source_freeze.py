"""Freeze upstream CYP recipe history and public-data loaders at GitHub commit."""

import json
import os
from pathlib import Path

import httpx

out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
repo = cfg["repo"]
api = f"https://api.github.com/repos/{repo}"
with httpx.Client(timeout=120, follow_redirects=True) as client:

    def get(url):
        r = client.get(url)
        r.raise_for_status()
        return r

    head = get(api + "/commits/main").json()
    sha = head["sha"]
    (out / "head.json").write_text(json.dumps(head, indent=2))
    print("commit", sha, flush=True)
    tree = get(api + f"/git/trees/{sha}?recursive=1").json()
    (out / "tree.json").write_text(json.dumps(tree))
    selected = []
    for t in tree["tree"]:
        p = t["path"]
        low = p.lower()
        if t["type"] == "blob" and (
            p.startswith("ml_pipelines/OpenADMET/cyp/")
            or ("cyp" in low and ("chembl" in low or "pubchem" in low))
            or ("public" in low and p.endswith(".py"))
            or p.endswith("multi_task.py")
        ):
            selected.append(p)
    for p in selected:
        f = out / "source" / p
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_bytes(get(f"https://raw.githubusercontent.com/{repo}/{sha}/{p}").content)
    history = get(
        api + "/commits?path=ml_pipelines/OpenADMET/cyp/cyp_chemprop_union.py&per_page=100"
    ).json()
    (out / "union-history.json").write_text(json.dumps(history, indent=2))
    for c in history:
        cs = c["sha"]
        r = get(
            f"https://raw.githubusercontent.com/{repo}/{cs}/ml_pipelines/OpenADMET/cyp/cyp_chemprop_union.py"
        )
        (out / f"union-{cs}.py").write_bytes(r.content)
        print(cs, c["commit"]["message"].splitlines()[0], flush=True)
    (out / "manifest.json").write_text(
        json.dumps(
            dict(
                repo=repo, commit=sha, files=selected, history_commits=[c["sha"] for c in history]
            ),
            indent=2,
        )
    )
