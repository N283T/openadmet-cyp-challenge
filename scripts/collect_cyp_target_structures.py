"""Discover one CYP target with TogoMCP and snapshot its experimental mmCIFs."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path

import httpx
import polars as pl


def mcp_request(process: subprocess.Popen[str], request_id: int, method: str, params: dict) -> dict:
    assert process.stdin is not None and process.stdout is not None
    process.stdin.write(
        json.dumps({"jsonrpc": "2.0", "id": request_id, "method": method, "params": params})
        + "\n"
    )
    process.stdin.flush()
    for line in process.stdout:
        message = json.loads(line)
        if message.get("id") == request_id:
            if "error" in message:
                raise RuntimeError(f"MCP {method} failed: {message['error']}")
            return message["result"]
    raise RuntimeError(f"MCP server closed during {method}")


def discover(config: dict) -> tuple[dict, dict]:
    executable = shutil.which("togo-mcp-local")
    if executable is None:
        raise RuntimeError("togo-mcp-local is unavailable on PATH")
    process = subprocess.Popen(
        [executable],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        bufsize=1,
    )
    try:
        initialized = mcp_request(
            process,
            1,
            "initialize",
            {
                "protocolVersion": "2025-03-26",
                "capabilities": {},
                "clientInfo": {"name": "cyp-target-structure-collection", "version": "0.1"},
            },
        )
        assert process.stdin is not None
        process.stdin.write('{"jsonrpc":"2.0","method":"notifications/initialized"}\n')
        process.stdin.flush()
        result = mcp_request(
            process,
            2,
            "tools/call",
            {
                "name": "search_pdb_entity",
                "arguments": {
                    "db": "pdb",
                    "query": config["togo_query"],
                    "limit": config["search_limit"],
                },
            },
        )
        if result.get("isError"):
            raise RuntimeError(f"TogoMCP search failed: {result}")
        payload = json.loads(result["content"][0]["text"])
        if "error" in payload:
            raise RuntimeError(f"PDBj search failed: {payload['error']}")
        return payload, initialized["serverInfo"]
    finally:
        process.terminate()
        process.wait(timeout=10)


def download(row: dict, destination: Path, template: str) -> dict:
    pdb_id = row["pdb_id"]
    response = httpx.get(
        template.format(pdb_id=pdb_id.upper()), timeout=60, follow_redirects=True
    )
    response.raise_for_status()
    data = response.content
    if not data.startswith(f"data_{pdb_id.upper()}".encode()):
        raise ValueError(f"Unexpected mmCIF header for {pdb_id}")
    path = destination / f"{pdb_id}.cif"
    path.write_bytes(data)
    return {
        **row,
        "source_url": str(response.url),
        "file": f"pdb/{pdb_id}.cif",
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
    }


def main() -> None:
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    payload, server = discover(config)
    (output_dir / "togo_pdbj_search.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    )
    hits = payload["results"]
    if payload.get("total") is not None and payload["total"] > len(hits):
        raise ValueError("TogoMCP search limit truncated the result set")
    protein_name = config["pdb_protein_name"].casefold()
    heme_names = [name.casefold() for name in config["pdb_heme_names"]]
    selected = []
    for hit in hits:
        contents = (hit.get("ligands") or "").casefold()
        if protein_name in contents and any(name in contents for name in heme_names):
            selected.append(
                {
                    "pdb_id": hit["id"].lower(),
                    "title": hit.get("title"),
                    "method": hit.get("method"),
                    "resolution_angstrom": hit.get("resolution"),
                    "ligands_summary": hit.get("ligands"),
                    "year": hit.get("year"),
                    "doi": hit.get("doi"),
                }
            )
    selected.sort(key=lambda row: row["pdb_id"])
    if not selected or len({row["pdb_id"] for row in selected}) != len(selected):
        raise ValueError("TogoMCP selection is empty or contains duplicate PDB IDs")

    destination = output_dir / "pdb"
    destination.mkdir()
    completed = []
    with ThreadPoolExecutor(max_workers=config["max_workers"]) as pool:
        futures = [
            pool.submit(download, row, destination, config["download_url_template"])
            for row in selected
        ]
        for future in as_completed(futures):
            completed.append(future.result())
            print(f"Downloaded {len(completed)}/{len(selected)}", flush=True)
    manifest = pl.DataFrame(completed).sort("pdb_id")
    manifest.write_csv(output_dir / "source_manifest.csv")
    manifest.write_parquet(output_dir / "source_manifest.parquet")
    summary = {
        "collected_at_utc": datetime.now(UTC).isoformat(),
        "target": config["target"],
        "accession": config["accession"],
        "togo_server": server,
        "search_query": config["togo_query"],
        "search_total": payload.get("total"),
        "search_returned": len(hits),
        "selection_rule": "PDBj summary contains the exact target protein and HEM names",
        "selected": manifest.height,
        "downloaded": manifest.height,
        "total_bytes": int(manifest["bytes"].sum()),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
