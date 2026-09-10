"""Persist the requested community article with the shared reference archive machinery."""

from datetime import UTC, datetime

import httpx
from archive_official_references import archive

from cyp_challenge.storage import ROOT, allocate, now, write_json


def main():
    identifier, log = allocate("reference-fetch")
    manifest = {"id": identifier, "created_at": now(), "entries": []}
    write_json(log / "intent.json", manifest)
    url = "https://supercowpowers.github.io/workbench/blogs/cyp_challenge/"
    try:
        with httpx.Client(timeout=60, follow_redirects=True) as client:
            path = archive(
                client,
                "cyp",
                "supercowpowers-working-notes",
                url,
                datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ"),
                collection="community",
            )
        manifest["entries"].append(path)
        index = ROOT / "references/community/README.md"
        if not index.exists():
            index.write_text(
                "# Community source archive\n\n"
                "Participant and third-party articles, separate from official sources.\n"
                "Snapshots preserve the author's claims; they are not independently verified.\n"
                "Each timestamped version includes article text, figures and hashed provenance.\n"
                "Raw source responses remain local and Git-ignored. Never replace old versions.\n"
                "Original copyright applies. Interpretations belong in project docs.\n\n"
                "Refresh: `pixi run python scripts/archive_community_references.py`.\n\n"
            )
        with index.open("a") as handle:
            relative = path.removeprefix("references/community/")
            handle.write(
                f"- [OpenADMET CYP Challenge — Working Notes]({relative}/snapshot.md) "
                f"— retrieved {manifest['created_at']}; [original]({url}).\n"
            )
        manifest["status"] = "completed"
    except Exception as exc:
        manifest.update(status="failed", error=repr(exc))
        raise
    finally:
        write_json(log / "result.json", manifest)


if __name__ == "__main__":
    main()
