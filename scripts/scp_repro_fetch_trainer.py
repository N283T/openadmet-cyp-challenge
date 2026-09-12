"""Snapshot published Workbench Chemprop implementation and training defaults."""

import json
import os
from pathlib import Path

import httpx


def main():
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    with httpx.Client(timeout=60, follow_redirects=True) as client:
        api = "https://api.github.com/repos/SuperCowPowers/workbench"
        revision = "4cf07a8cb4366eb29111fb42396336b20b69ada0"
        response = client.get(api + f"/git/trees/{revision}?recursive=1")
        response.raise_for_status()
        tree = response.json()
        (out / "tree.json").write_text(json.dumps(tree, indent=2))
        paths = [
            r["path"]
            for r in tree["tree"]
            if r["type"] == "blob"
            and r["path"].endswith(".py")
            and ("chemprop" in r["path"].lower() or r["path"].endswith("multi_task.py"))
            and r["path"].startswith("src/")
        ]
        for path in paths:
            r = client.get(
                f"https://raw.githubusercontent.com/SuperCowPowers/workbench/{tree['sha']}/{path}"
            )
            r.raise_for_status()
            dest = out / path
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(r.content)
        template = "src/workbench/model_scripts/chemprop/chemprop.template"
        r = client.get(
            f"https://raw.githubusercontent.com/SuperCowPowers/workbench/{tree['sha']}/{template}"
        )
        r.raise_for_status()
        (out / "chemprop.template").write_bytes(r.content)
        paths.append(template)
        (out / "manifest.json").write_text(
            json.dumps(dict(commit=tree["sha"], files=paths), indent=2)
        )
        print(json.dumps(dict(commit=tree["sha"], files=len(paths))), flush=True)


if __name__ == "__main__":
    main()
