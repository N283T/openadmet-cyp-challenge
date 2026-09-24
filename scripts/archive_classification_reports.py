"""Archive the three public classification reports selected for CYP2D6 review."""

import hashlib
import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

import httpx
from bs4 import BeautifulSoup

from cyp_challenge.storage import ROOT, digest, inventory, now, write_json

DASHBOARD_URL = "https://openadmet-cyp-dashboard.vercel.app/operating_dashboard"
DOC_ID = "1ohRu6s3tod3sjRF2gV3mYdoe9cHzgn-dYE5ruSyBrWg"
DOC_URL = f"https://docs.google.com/document/d/{DOC_ID}/edit?usp=sharing"
DOC_EXPORT_URL = f"https://docs.google.com/document/d/{DOC_ID}/export?format=txt"
JEREMY_URL = "https://github.com/jeremycheminf/openadmet_scripts/tree/main/CYP_Challenge"
JEREMY_REPO = ROOT.parents[1] / "jeremycheminf/openadmet_scripts"
JEREMY_FILES = (
    "CYP_Challenge/README.md",
    "CYP_Challenge/scripts/02_baseline_lgb.py",
    "CYP_Challenge/scripts/05_run_tabicl_tdi.py",
    "CYP_Challenge/scripts/07_build_submission.py",
)


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def write_snapshot(path: Path, content: str, metadata: dict) -> None:
    path.mkdir(parents=True, exist_ok=False)
    write_json(path / "intent.json", metadata)
    (path / "snapshot.md").write_text(content)
    metadata.update(
        snapshot_sha256=digest(path / "snapshot.md"),
        files=inventory(path),
        status="completed",
    )
    write_json(path / "meta.json", metadata)


def archive_vaishnavi(client: httpx.Client, stamp: str) -> Path:
    path = ROOT / "references/community/cyp/vaishnavi53-methodology" / stamp
    response = client.get(DOC_EXPORT_URL)
    response.raise_for_status()
    raw = response.content
    (path.parent).mkdir(parents=True, exist_ok=True)
    path.mkdir(parents=True, exist_ok=False)
    (path / "source.txt").write_bytes(raw)
    text = raw.decode("utf-8-sig").strip()
    metadata = {
        "url": DOC_URL,
        "export_url": DOC_EXPORT_URL,
        "title": "Methodology Report_CYP",
        "participant": "vaishnavi53",
        "fetched_at": now(),
        "source_sha256": sha256_bytes(raw),
        "extraction": "Public Google Doc plain-text export; original order preserved.",
    }
    write_json(path / "intent.json", metadata)
    header = f"Source: {DOC_URL} | Retrieved: {metadata['fetched_at']}\n\n"
    (path / "snapshot.md").write_text(header + text + "\n")
    metadata.update(
        snapshot_sha256=digest(path / "snapshot.md"),
        files=inventory(path),
        status="completed",
    )
    write_json(path / "meta.json", metadata)
    return path


def archive_dashboard(client: httpx.Client, stamp: str) -> Path:
    path = ROOT / "references/community/cyp/450nm-campaign-console" / stamp
    response = client.get(DASHBOARD_URL)
    response.raise_for_status()
    raw = response.content
    soup = BeautifulSoup(response.text, "html.parser")
    payload_node = soup.select_one("#payload")
    if payload_node is None or not payload_node.string:
        raise ValueError("Dashboard JSON payload was not found")
    payload = json.loads(payload_node.string)
    boards = payload["subs"]["board"]["boards"]
    selected = []
    for experiment in payload["experiments"]:
        searchable = json.dumps(experiment, ensure_ascii=False).lower()
        if experiment.get("track") == "tdi" and "2d6" in searchable:
            selected.append(experiment)

    metadata = {
        "url": DASHBOARD_URL,
        "title": "CYP Campaign Console",
        "participant": "450nm",
        "fetched_at": now(),
        "dashboard_generated": payload.get("generated"),
        "source_sha256": sha256_bytes(raw),
        "extraction": (
            "Classification leaderboard rows plus every source experiment record "
            "whose track is tdi and whose fields mention 2D6; full HTML retained locally."
        ),
        "selected_experiment_ids": [item["id"] for item in selected],
    }
    path.mkdir(parents=True, exist_ok=False)
    (path / "source.html").write_bytes(raw)
    write_json(path / "intent.json", metadata)
    lines = [
        f"Source: {DASHBOARD_URL} | Retrieved: {metadata['fetched_at']}",
        "",
        "# CYP Campaign Console — classification and CYP2D6 TDI source extract",
        "",
        f"Dashboard generated: {payload.get('generated')}",
        "",
        f"Campaign status: {payload['campaign'].get('status', '')}",
        "",
        "## Classification leaderboard snapshot",
        "",
    ]
    for key in ("classification/overall", "classification/cyp2d6", "classification/cyp3a4"):
        board = boards[key]
        lines.extend(
            [
                f"### {key}",
                "",
                "```json",
                json.dumps({"metric": board.get("metric"), "top10": board["top10"]}, indent=2),
                "```",
                "",
            ]
        )
    lines.extend(["## CYP2D6 TDI experiment records", ""])
    for experiment in selected:
        lines.extend([f"### {experiment['id']} — {experiment['title']}", ""])
        for key in (
            "date",
            "track",
            "isoform",
            "verdict",
            "hypothesis",
            "experiment",
            "results",
            "interpretation",
            "artifact",
        ):
            value = experiment.get(key)
            if value not in (None, ""):
                lines.extend([f"**{key}:** {value}", ""])
    (path / "snapshot.md").write_text("\n".join(lines).rstrip() + "\n")
    metadata.update(
        snapshot_sha256=digest(path / "snapshot.md"),
        files=inventory(path),
        status="completed",
    )
    write_json(path / "meta.json", metadata)
    return path


def git_output(*args: str) -> bytes:
    return subprocess.run(
        ["git", "-C", str(JEREMY_REPO), *args], check=True, capture_output=True
    ).stdout


def archive_jeremy(stamp: str) -> Path:
    path = ROOT / "references/community/cyp/jeremy-reproduction-kit" / stamp
    commit = git_output("rev-parse", "HEAD").decode().strip()
    commit_date = git_output("show", "-s", "--format=%cI", commit).decode().strip()
    files = []
    for name in JEREMY_FILES:
        content = git_output("show", f"{commit}:{name}")
        files.append({"path": name, "sha256": sha256_bytes(content), "bytes": len(content)})
    readme = git_output("show", f"{commit}:CYP_Challenge/README.md").decode()
    metadata = {
        "url": JEREMY_URL,
        "title": "OpenADMET CYP Inhibition Challenge — submission reproduction kit",
        "participant": "jeremy",
        "fetched_at": now(),
        "local_repository": str(JEREMY_REPO),
        "commit": commit,
        "commit_date": commit_date,
        "tracked_source_files": files,
        "extraction": (
            "README preserved as the report; relevant executable files are pinned by "
            "commit and SHA-256 but not duplicated because the repository has no license file."
        ),
    }
    header = f"Source: {JEREMY_URL} | Commit: {commit} | Retrieved: {metadata['fetched_at']}\n\n"
    write_snapshot(path, header + readme.rstrip() + "\n", metadata)
    return path


def update_index(entries: list[tuple[str, Path, str]], retrieved: str) -> None:
    index = ROOT / "references/community/README.md"
    with index.open("a") as handle:
        for title, path, url in entries:
            relative = path.relative_to(ROOT / "references/community")
            handle.write(
                f"- [{title}]({relative}/snapshot.md) — retrieved {retrieved}; [original]({url}).\n"
            )


def main() -> None:
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    with httpx.Client(timeout=90, follow_redirects=True) as client:
        vaishnavi = archive_vaishnavi(client, stamp)
        dashboard = archive_dashboard(client, stamp)
    jeremy = archive_jeremy(stamp)
    retrieved = now()
    update_index(
        [
            ("Methodology Report_CYP — vaishnavi53", vaishnavi, DOC_URL),
            ("CYP Campaign Console — 450nm", dashboard, DASHBOARD_URL),
            ("CYP submission reproduction kit — jeremy", jeremy, JEREMY_URL),
        ],
        retrieved,
    )
    print(f"Archived classification reports at {stamp}")


if __name__ == "__main__":
    main()
