"""Archive the public regression reports linked from the top of the 2026-10-29 Live board."""

import hashlib
from datetime import UTC, datetime

import httpx

from cyp_challenge.storage import ROOT, allocate, digest, inventory, now, write_json

STIR_BAR_REPO = "lachrymator/openadmet-cyp-challenge-public"
STIR_BAR_URL = f"https://github.com/{STIR_BAR_REPO}"
STIR_BAR_FIGURES = (
    "leaderboard_history.png",
    "cv_vs_blind.png",
    "standing_regression.png",
    "standing_tdi.png",
    "submission_delta.png",
    "chemical_space_tsne.png",
)
# Linked reports that held no method content at retrieval; recorded, not archived.
PLACEHOLDERS = (
    (
        "Rubyy",
        "https://huggingface.co/basedonAI/openadmet-cyp-inhibition-rubyy/raw/main/README.md",
        "Model-card template; says it will be updated when the competition closes.",
    ),
    (
        "hamza322",
        "https://docs.google.com/document/d/141YL5Gc6BfROdl5wubPtAyJQAywspCCqzcvYyvS0Yn8"
        "/export?format=txt",
        "Team names and 'preprint upcoming'; model and data sections say 'coming soon'.",
    ),
    (
        "rasayan-labs",
        "https://api.github.com/repos/whis9/rasayan-cyp",
        "Repository exists but is empty (size 0).",
    ),
)


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def archive_stir_bar(client: httpx.Client, stamp: str) -> tuple[str, dict]:
    path = ROOT / "references/community/cyp/stir-bar-writeup" / stamp
    commit = client.get(f"https://api.github.com/repos/{STIR_BAR_REPO}/commits/main")
    commit.raise_for_status()
    sha = commit.json()["sha"]
    raw_base = f"https://raw.githubusercontent.com/{STIR_BAR_REPO}/{sha}"
    readme = client.get(f"{raw_base}/README.md")
    readme.raise_for_status()
    path.mkdir(parents=True, exist_ok=False)
    (path / "source.md").write_bytes(readme.content)
    figures = []
    for name in STIR_BAR_FIGURES:
        response = client.get(f"{raw_base}/{name}")
        response.raise_for_status()
        (path / "assets").mkdir(exist_ok=True)
        (path / "assets" / name).write_bytes(response.content)
        figures.append({"path": name, "sha256": sha256_bytes(response.content)})
    metadata = {
        "url": STIR_BAR_URL,
        "title": "OpenADMET CYP Inhibition Challenge write-up",
        "participant": "stir_bar",
        "fetched_at": now(),
        "commit": sha,
        "commit_date": commit.json()["commit"]["committer"]["date"],
        "source_sha256": sha256_bytes(readme.content),
        "figures": figures,
        "extraction": "README at the pinned commit preserved verbatim; figures kept as assets.",
    }
    write_json(path / "intent.json", metadata)
    header = f"Source: {STIR_BAR_URL} | Commit: {sha} | Retrieved: {metadata['fetched_at']}\n\n"
    (path / "snapshot.md").write_text(header + readme.text.rstrip() + "\n")
    metadata.update(
        snapshot_sha256=digest(path / "snapshot.md"), files=inventory(path), status="completed"
    )
    write_json(path / "meta.json", metadata)
    return str(path.relative_to(ROOT)), metadata


def check_placeholders(client: httpx.Client) -> list[dict]:
    records = []
    for participant, url, observed in PLACEHOLDERS:
        response = client.get(url)
        response.raise_for_status()
        records.append(
            {
                "participant": participant,
                "url": url,
                "sha256": sha256_bytes(response.content),
                "bytes": len(response.content),
                "observed": observed,
            }
        )
    return records


def main() -> None:
    identifier, log = allocate("reference-fetch")
    manifest = {"id": identifier, "created_at": now(), "entries": []}
    write_json(log / "intent.json", manifest)
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    try:
        with httpx.Client(timeout=90, follow_redirects=True) as client:
            path, metadata = archive_stir_bar(client, stamp)
            placeholders = check_placeholders(client)
        manifest["entries"].append(path)
        manifest["placeholders"] = placeholders
        index = ROOT / "references/community/README.md"
        relative = path.removeprefix("references/community/")
        with index.open("a") as handle:
            handle.write(
                f"- [OpenADMET CYP Inhibition Challenge write-up — stir_bar]"
                f"({relative}/snapshot.md) — retrieved {metadata['fetched_at']}; "
                f"[original]({STIR_BAR_URL}).\n"
            )
            for record in placeholders:
                handle.write(
                    f"- {record['participant']} report link checked {metadata['fetched_at']}: "
                    f"no method content ({record['observed']}) [{record['url']}].\n"
                )
        manifest["status"] = "completed"
    except Exception as exc:
        manifest.update(status="failed", error=repr(exc))
        raise
    finally:
        write_json(log / "result.json", manifest)


if __name__ == "__main__":
    main()
