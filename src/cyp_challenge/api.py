"""Persist public API contracts and HTTP downloads, including failures."""

import httpx

from .database import connect, register_asset
from .storage import ROOT, allocate, inventory, now, write_json

SPACE = "https://openadmet-cyp-challenge.hf.space"
HF = "https://huggingface.co/spaces/openadmet/cyp-challenge"


def inspect_api():
    with connect():
        pass
    identifier, path = allocate("api")
    manifest = {"id": identifier, "created_at": now(), "requests": []}
    write_json(path / "intent.json", manifest)
    urls = {
        "info.json": f"{SPACE}/gradio_api/info",
        "space.json": "https://huggingface.co/api/spaces/openadmet/cyp-challenge",
        "config.json": f"{SPACE}/config",
    }
    try:
        with httpx.Client(timeout=60, follow_redirects=True) as client:
            for filename, url in urls.items():
                response = client.get(url)
                (path / filename).write_bytes(response.content)
                manifest["requests"].append(
                    {
                        "url": url,
                        "status": response.status_code,
                        "file": filename,
                        "fetched_at": now(),
                    }
                )
                response.raise_for_status()
            revision = __import__("json").loads((path / "space.json").read_text())["sha"]
            manifest["space_revision"] = revision
            for filename in ("config.py", "submission.py", "leaderboards.py"):
                url = f"{HF}/raw/{revision}/{filename}"
                response = client.get(url)
                (path / filename).write_bytes(response.content)
                manifest["requests"].append(
                    {
                        "url": url,
                        "status": response.status_code,
                        "file": filename,
                        "fetched_at": now(),
                    }
                )
                response.raise_for_status()
        manifest["status"] = "completed"
    except Exception as exc:
        manifest.update(status="failed", error=repr(exc))
        raise
    finally:
        manifest["files"] = inventory(path)
        write_json(path / "manifest.json", manifest)
        register_asset(identifier, "api", path, manifest)
        print(path)


def register_raw():
    with connect():
        pass
    identifier, path = allocate("data-imports")
    manifest = {
        "id": identifier,
        "created_at": now(),
        "source": "Existing local CSVs; original download revision is unknown",
        "root": str(ROOT / "data/raw"),
        "files": inventory(ROOT / "data/raw"),
    }
    write_json(path / "manifest.json", manifest)
    register_asset(identifier, "raw", path, manifest)
    print(path)


def leaderboard(track, endpoint, board="live"):
    """Resolve a download endpoint from the selected board and retain the actual CSV."""

    from gradio_client import Client

    from .tracks import ENDPOINTS

    if endpoint != "MA" and endpoint not in ENDPOINTS[track]:
        raise ValueError("Endpoint does not belong to this track")
    if board not in ("live", "interim", "final"):
        raise ValueError("Board must be live, interim or final")
    with connect():
        pass
    identifier, path = allocate("api", track)
    manifest = {
        "id": identifier,
        "created_at": now(),
        "track": track,
        "endpoint": endpoint,
        "board": board,
        "operation": "leaderboard-download",
        "space": SPACE,
    }
    write_json(path / "intent.json", manifest)
    try:
        client = Client(
            SPACE, download_files=path / "downloads", analytics_enabled=False, verbose=False
        )
        config = client.config
        write_json(path / "config.json", config)
        write_json(path / "info.json", client.view_api(return_format="dict", print_info=False))
        components = {item["id"]: item for item in config["components"]}
        ancestors = {}

        def walk(node, labels):
            component = components.get(node["id"], {})
            if component.get("type") in ("tab", "tabitem"):
                labels = [*labels, component.get("props", {}).get("label", "")]
            ancestors[node["id"]] = labels
            for child in node.get("children", []):
                walk(child, labels)

        walk(config["layout"], [])
        track_label = {
            "regression": "Direct Inhibition Prediction (Regression)",
            "classification": "TDI Prediction (Classification)",
        }[track]
        board_label = {
            "live": "🏎️ Live Leaderboard",
            "interim": "⏱️ Interim Leaderboard",
            "final": "🏁 Final Leaderboard",
        }[board]
        endpoint_label = "🏆 Overall" if endpoint == "MA" else endpoint.replace("_", " ")
        matches = []
        for dep in config["dependencies"]:
            for output in dep["outputs"]:
                labels = ancestors.get(output, [])
                if (
                    components[output]["type"] == "downloadbutton"
                    and board_label in labels
                    and track_label in labels
                    and endpoint_label in labels
                    and not dep["inputs"]
                    and dep.get("api_name")
                ):
                    matches.append(dep["api_name"])
        if len(matches) != 1:
            raise ValueError(f"Cannot unambiguously resolve download endpoint: {matches}")
        manifest["api_name"] = "/" + matches[0]
        response = client.predict(api_name=manifest["api_name"])
        write_json(path / "response.json", {"result": response})
        from pathlib import Path

        if not isinstance(response, str) or not Path(response).is_file():
            raise ValueError("No downloadable leaderboard file returned")
        if not Path(response).resolve().is_relative_to(path.resolve()):
            raise ValueError("Client returned a file outside the persistent download directory")
        manifest["status"] = "completed"
    except Exception as exc:
        manifest.update(status="failed", error=repr(exc))
        raise
    finally:
        manifest["files"] = inventory(path)
        write_json(path / "manifest.json", manifest)
        register_asset(identifier, "api", path, manifest)
        print(path)
    return path
