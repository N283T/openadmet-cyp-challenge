"""Archive the selected official pages as versioned Markdown, figures, and provenance."""

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urljoin, urlparse

import httpx
from bs4 import BeautifulSoup
from markdownify import markdownify
from psycopg.types.json import Jsonb

from cyp_challenge.database import connect
from cyp_challenge.storage import ROOT, allocate, digest, inventory, now, write_json

GHOST = "https://openadmet.ghost.io/"
SOURCES = [
    ("cyp", "challenge", "https://huggingface.co/spaces/openadmet/cyp-challenge"),
    ("cyp", "announcement", GHOST + "announcing-openadmets-cyp-inhibition-blind-challenge/"),
    ("cyp", "structural-alerts", GHOST + "do-our-cyp-structural-alerts-actually-work/"),
    ("cyp", "underway", GHOST + "openadmets-cyp-challenge-is-underway/"),
    ("cyp", "assay-dataset", "https://openadmet.github.io/octant-cyp-inhib-blog-post/"),
    (
        "pxr",
        "retrospective-part-1",
        GHOST + "dont-look-back-in-error-what-we-learned-predicting-pxr-induction-part-i/",
    ),
]


def fetch(client, url, target, requests):
    response = client.get(url)
    target.write_bytes(response.content)
    requests.append(
        {
            "url": url,
            "final_url": str(response.url),
            "status": response.status_code,
            "fetched_at": now(),
            "file": target.name,
            "sha256": digest(target),
        }
    )
    response.raise_for_status()
    return response


def extract_article(html, url):
    soup = BeautifulSoup(html, "html.parser")
    title = soup.find("h1").get_text(" ", strip=True)
    body = (
        soup.select_one("section.gh-content")
        or soup.select_one("main#quarto-document-content")
        or soup.select_one("article.md-content__inner")
    )
    if body is None:
        raise ValueError(f"Unknown article layout: {url}")
    metadata = {"title": title, "authors": [], "published": None, "modified": None}
    for node in soup.select('script[type="application/ld+json"]'):
        data = json.loads(node.string or node.get_text())
        if isinstance(data, dict) and data.get("@type") in ("Article", "BlogPosting"):
            metadata.update(published=data.get("datePublished"), modified=data.get("dateModified"))
            authors = data.get("author", [])
            authors = authors if isinstance(authors, list) else [authors]
            metadata["authors"] = [a.get("name") for a in authors if isinstance(a, dict)]
    explicit_authors = soup.select_one(".doc-authors")
    if explicit_authors:
        metadata["authors_display"] = explicit_authors.get_text(" ", strip=True)
    date = soup.select_one(".quarto-title-meta .date")
    if not metadata["published"] and date:
        metadata["published_display"] = date.get_text(" ", strip=True)
    # Keep the article header image/caption and subtitle, but omit navigation and subscription UI.
    prefix = "".join(str(x) for x in soup.select(".gh-article-excerpt, .gh-article-image"))
    body = BeautifulSoup(prefix + str(body), "html.parser")
    for node in body.select("script, style, .anchorjs-link, .code-copy-button, .headerlink"):
        node.decompose()
    for node in body.select("a[href]"):
        node["href"] = urljoin(url, node["href"])
    return body, metadata


def archive(client, group, slug, url, stamp, collection="official"):
    path = ROOT / "references" / collection / group / slug / stamp
    path.mkdir(parents=True, exist_ok=False)
    metadata = {"url": url, "group": group, "slug": slug, "fetched_at": now(), "requests": []}
    write_json(path / "intent.json", metadata)
    try:
        if slug == "challenge":
            response = fetch(
                client,
                "https://openadmet-cyp-challenge.hf.space/config",
                path / "source.json",
                metadata["requests"],
            )
            config = response.json()
            space = fetch(
                client,
                "https://huggingface.co/api/spaces/openadmet/cyp-challenge",
                path / "source-space.json",
                metadata["requests"],
            ).json()
            metadata.update(
                title="OpenADMET CYP Blind Challenge",
                authors=["OpenADMET"],
                space_revision=space["sha"],
                extraction="Live Gradio Markdown components; excludes leaderboard widgets",
            )
            chunks = [
                c["props"]["value"]
                for c in config["components"]
                if c["type"] == "markdown" and c.get("props", {}).get("value")
            ]
            content = "\n\n---\n\n".join(chunks)
        else:
            response = fetch(client, url, path / "source.html", metadata["requests"])
            body, fields = extract_article(response.text, str(response.url))
            metadata.update(
                fields, extraction="Article body converted to Markdown; figures saved locally"
            )
            assets = path / "assets"
            images = []
            for img in body.select("img"):
                src = img.get("src") or img.get("data-src")
                if not src:
                    continue
                image_url = urljoin(str(response.url), src)
                if urlparse(image_url).scheme not in ("http", "https"):
                    continue
                assets.mkdir(exist_ok=True)
                suffix = Path(urlparse(image_url).path).suffix.lower()
                if suffix not in (".png", ".jpg", ".jpeg", ".webp", ".gif", ".svg"):
                    suffix = ".img"
                filename = hashlib.sha256(image_url.encode()).hexdigest()[:16] + suffix
                destination = assets / filename
                if not destination.exists():
                    fetch(client, image_url, destination, metadata["requests"])
                img["src"] = f"assets/{filename}"
                images.append(
                    {
                        "url": image_url,
                        "file": f"assets/{filename}",
                        "sha256": digest(destination),
                        "alt": img.get("alt", ""),
                    }
                )
            metadata["images"] = images
            content = markdownify(str(body), heading_style="ATX", bullets="-", wrap=False)
        if len(content.strip()) < 500:
            raise ValueError("Unexpectedly short article extraction")
        header = f"Source: {url} | Retrieved: {metadata['fetched_at']}\n\n"
        (path / "snapshot.md").write_text(header + f"# {metadata['title']}\n\n" + content + "\n")
        metadata["snapshot_sha256"] = digest(path / "snapshot.md")
        metadata["status"] = "completed"
    except Exception as exc:
        metadata.update(status="failed", error=repr(exc))
        raise
    finally:
        metadata["files"] = inventory(path)
        write_json(path / "meta.json", metadata)
        with connect() as conn:
            conn.execute(
                "INSERT INTO reference_snapshots(id,url,path,metadata) VALUES (%s,%s,%s,%s)",
                (
                    f"{group}/{slug}/{stamp}"
                    if collection == "official"
                    else f"{collection}/{group}/{slug}/{stamp}",
                    url,
                    str(path.relative_to(ROOT)),
                    Jsonb(metadata),
                ),
            )
    print(
        f"{group}/{slug}: {len(content)} characters, {len(metadata.get('images', []))} images",
        flush=True,
    )
    return str(path.relative_to(ROOT))


def update_index():
    root = ROOT / "references/official"
    lines = [
        "# Official source archive",
        "",
        "Unabridged article-body snapshots for this private research repository.",
        "Original language is preserved. Interpretations belong in Issues, not snapshots.",
        "Each version has `snapshot.md`, `meta.json`, and local figures under `assets/`.",
        "Raw `source*` responses remain local and are excluded from Git.",
        "Updates create a new timestamped version; previous versions are never replaced.",
        "",
        "The Space snapshot includes its live Markdown instructions/FAQ, not dynamic LB tables.",
        "Use the API archive for leaderboard data. Snapshot dates are not publication dates.",
        "Original copyright applies; review redistribution before public release.",
        "",
        "Refresh: `pixi run python scripts/archive_official_references.py`.",
        "",
    ]
    for group in ("cyp", "pxr"):
        lines += [
            f"## {group.upper()}",
            "",
            "| Source | Snapshot (UTC) | Original |",
            "|---|---|---|",
        ]
        for file in sorted((root / group).glob("*/*/meta.json")):
            meta = json.loads(file.read_text())
            if meta["status"] != "completed":
                continue
            rel = file.parent.relative_to(root)
            title = meta["title"].replace("|", "\\|")
            lines.append(
                f"| [{title}]({rel}/snapshot.md) | {meta['fetched_at']} | [Web]({meta['url']}) |"
            )
        lines.append("")
    (root / "README.md").write_text("\n".join(lines))


def main():
    with connect():
        pass
    identifier, log = allocate("reference-fetch")
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    manifest = {"id": identifier, "created_at": now(), "entries": []}
    write_json(log / "intent.json", manifest)
    try:
        with httpx.Client(timeout=60, follow_redirects=True) as client:
            for group, slug, url in SOURCES:
                manifest["entries"].append(archive(client, group, slug, url, stamp))
        update_index()
        manifest["status"] = "completed"
    except Exception as exc:
        manifest.update(status="failed", error=repr(exc))
        raise
    finally:
        write_json(log / "result.json", manifest)


if __name__ == "__main__":
    main()
