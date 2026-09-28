"""Discover CYP3A4 structures with TogoMCP and archive their source mmCIFs."""

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
        init = mcp_request(
            process,
            1,
            "initialize",
            {
                "protocolVersion": "2025-03-26",
                "capabilities": {},
                "clientInfo": {"name": "cyp-structure-source-collection", "version": "0.1"},
            },
        )
        assert process.stdin is not None
        process.stdin.write('{"jsonrpc":"2.0","method":"notifications/initialized"}\n')
        process.stdin.flush()
        guide = mcp_request(
            process, 2, "tools/call", {"name": "TogoMCP_Usage_Guide", "arguments": {}}
        )
        if guide.get("isError"):
            raise RuntimeError("TogoMCP usage guide failed")
        result = mcp_request(
            process,
            3,
            "tools/call",
            {
                "name": "search_pdb_entity",
                "arguments": {
                    "db": "pdb",
                    "query": config["query"],
                    "limit": config["search_limit"],
                },
            },
        )
        if result.get("isError"):
            raise RuntimeError(f"TogoMCP search failed: {result}")
        payload = json.loads(result["content"][0]["text"])
        if "error" in payload:
            raise RuntimeError(f"PDBj search failed: {payload['error']}")
        return payload, init["serverInfo"]
    finally:
        process.terminate()
        process.wait(timeout=10)


def download(row: dict, destination: Path, url_template: str) -> dict:
    pdb_id = row["pdb_id"]
    url = url_template.format(pdb_id=pdb_id.upper())
    path = destination / f"{pdb_id}.cif"
    if path.exists():
        raise FileExistsError(path)
    for attempt in range(3):
        try:
            response = httpx.get(url, timeout=45, follow_redirects=True)
            response.raise_for_status()
            data = response.content
            if not data.startswith(f"data_{pdb_id.upper()}".encode()):
                raise ValueError(f"Unexpected mmCIF header for {pdb_id}")
            path.write_bytes(data)
            return {
                **row,
                "source_url": str(response.url),
                "file": f"pdb/{pdb_id}.cif",
                "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
            }
        except (httpx.HTTPError, ValueError):
            if attempt == 2:
                raise
    raise AssertionError("unreachable")


def main() -> None:
    output_dir = Path(os.environ["CYP_OUTPUT_DIR"])
    input_dir = Path(os.environ["CYP_INPUT_DIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    targets = pl.read_csv(input_dir / "structure_targets.csv")
    if targets.height != 20 or targets["Molecule_Name"].n_unique() != 20:
        raise ValueError("Expected 20 unique blinded structure targets")
    sequence = (input_dir / "cyp3a4_sequence.fasta").read_text()
    if "P08684" not in sequence.splitlines()[0]:
        raise ValueError("Expected UniProt P08684 CYP3A4 sequence")

    payload, server_info = discover(config)
    (output_dir / "togo_pdbj_search.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    )
    hits = payload["results"]
    if payload.get("total") is not None and payload["total"] > len(hits):
        raise ValueError("Search limit truncated the PDBj results")
    protein_name = config["protein_name"].casefold()
    heme_name = config["heme_name"].casefold()
    selected = []
    for hit in hits:
        ligands = (hit.get("ligands") or "").casefold()
        if protein_name in ligands and heme_name in ligands:
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
    selected = sorted(selected, key=lambda row: row["pdb_id"])
    if len({row["pdb_id"] for row in selected}) != len(selected):
        raise ValueError("Duplicate PDB IDs in discovery results")
    if not selected:
        raise ValueError("No CYP3A4/HEM structures found")

    destination = output_dir / "pdb"
    destination.mkdir()
    completed = []
    with ThreadPoolExecutor(max_workers=config["max_workers"]) as pool:
        futures = {
            pool.submit(download, row, destination, config["download_url_template"]): row["pdb_id"]
            for row in selected
        }
        for future in as_completed(futures):
            completed.append(future.result())
            if len(completed) % 20 == 0 or len(completed) == len(selected):
                print(f"Downloaded {len(completed)}/{len(selected)} mmCIF files", flush=True)

    manifest = pl.DataFrame(sorted(completed, key=lambda row: row["pdb_id"]))
    manifest.write_parquet(output_dir / "source_manifest.parquet")
    manifest.write_csv(output_dir / "source_manifest.csv")
    summary = {
        "collected_at_utc": datetime.now(UTC).isoformat(),
        "togo_server": server_info,
        "search_query": config["query"],
        "search_total": payload.get("total"),
        "search_returned": len(hits),
        "selected_cyp3a4_with_heme": len(selected),
        "downloaded": len(completed),
        "total_bytes": sum(row["bytes"] for row in completed),
        "target_count": targets.height,
        "official_dataset_revision": config["official_dataset_revision"],
        "official_tutorial_revision": config["official_tutorial_revision"],
        "selection_rule": "PDBj ligands summary contains exact CYP3A4 protein and HEM names",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
