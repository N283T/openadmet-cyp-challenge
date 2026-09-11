"""Retrieve the full panel by explicit SID batches, then keyed compound structures."""

import io
import json
import os
import shutil
import time
from pathlib import Path

import httpx
import polars as pl

from cyp_challenge.storage import digest, now, write_json


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    shutil.copyfile(next(inp.glob("*-aid1851-description.json")), out / "aid1851-description.json")
    base = cfg["base_url"]
    records = []
    with httpx.Client(timeout=180, follow_redirects=True) as client:

        def get(name, url):
            for attempt in range(4):
                started = now()
                try:
                    response = client.get(url)
                except httpx.TransportError as exc:
                    write_json(
                        out / f"{name}.attempt{attempt}.error.json",
                        dict(url=url, retrieved_at=started, error=repr(exc)),
                    )
                    if attempt == 3:
                        raise
                    time.sleep(3 * (attempt + 1))
                    continue
                path = out / f"{name}.attempt{attempt}"
                path.write_bytes(response.content)
                record = dict(
                    name=name,
                    url=url,
                    retrieved_at=started,
                    status=response.status_code,
                    sha256=digest(path),
                    bytes=path.stat().st_size,
                    attempt=attempt,
                )
                records.append(record)
                write_json(out / f"{name}.attempt{attempt}.metadata.json", record)
                time.sleep(cfg["request_interval_seconds"])
                if response.is_success:
                    return response.content
                if attempt == 3 or response.status_code not in [429, 500, 502, 503, 504]:
                    response.raise_for_status()
                time.sleep(3 * (attempt + 1))
            raise RuntimeError("No response")

        sids = sorted(set(map(int, get("sids.txt", f"{base}/assay/aid/1851/sids/TXT").split())))
        tables = []
        for start in range(0, len(sids), cfg["sid_batch_size"]):
            batch = sids[start : start + cfg["sid_batch_size"]]
            content = get(
                f"assay-{start:05d}.csv",
                f"{base}/assay/aid/1851/CSV?sid=" + ",".join(map(str, batch)),
            )
            table = pl.read_csv(io.BytesIO(content), infer_schema=False).filter(
                pl.col("PUBCHEM_SID").cast(pl.Int64, strict=False).is_not_null()
            )
            assert set(table["PUBCHEM_SID"].cast(pl.Int64)) == set(batch)
            tables.append(table)
            print(
                f"Assay SIDs {min(start + len(batch), len(sids))}/{len(sids)}; {table.height} rows",
                flush=True,
            )
        raw = pl.concat(tables)
        raw.write_parquet(out / "pubchem-raw.parquet")
        print("COLUMNS", raw.columns, flush=True)
        cids = sorted(set(raw["PUBCHEM_CID"].cast(pl.Int64, strict=False).drop_nulls()))
        structures = []
        for start in range(0, len(cids), cfg["cid_batch_size"]):
            batch = cids[start : start + cfg["cid_batch_size"]]
            url = (
                f"{base}/compound/cid/"
                + ",".join(map(str, batch))
                + "/property/IsomericSMILES,InChIKey/CSV"
            )
            content = get(f"structures-{start:05d}.csv", url)
            table = pl.read_csv(io.BytesIO(content))
            assert set(table["CID"]) == set(batch)
            structures.append(table)
            if start % (cfg["cid_batch_size"] * 10) == 0:
                print(f"Compound structures {start + len(batch)}/{len(cids)}", flush=True)
        combined = pl.concat(structures)
        assert combined["CID"].n_unique() == combined.height == len(cids)
        combined.write_parquet(out / "pubchem-structures.parquet")
    write_json(out / "requests.json", records)
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            sids=len(sids),
            cids=len(cids),
            rows=raw.height,
            columns=raw.columns,
            structure_columns=combined.columns,
        ),
    )


if __name__ == "__main__":
    main()
