"""Acquire five allowlisted local AZ ADME assays under a read-only snapshot."""

import json
import os
from pathlib import Path

import polars as pl
import psycopg
from psycopg.rows import dict_row

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
query = next(inp.glob("*-adme_az_extract.sql")).read_text()
with psycopg.connect(
    dbname=cfg["dbname"],
    host=cfg["host"],
    port=cfg["port"],
    row_factory=dict_row,
    options=(
        "-c default_transaction_read_only=on -c statement_timeout=90000 "
        "-c max_parallel_workers_per_gather=0"
    ),
) as conn:
    conn.execute("SET TRANSACTION ISOLATION LEVEL REPEATABLE READ READ ONLY")
    version = conn.execute("SELECT name,creation_date::text FROM version").fetchall()
    source = conn.execute("SELECT * FROM source WHERE src_id=27").fetchall()
    rows = conn.execute(query).fetchall()
d = pl.DataFrame(rows, infer_schema_length=None)
assert d["activity_id"].n_unique() == d.height and set(d["source_assay_id"]) == set(
    cfg["assay_ids"]
)
d.write_parquet(out / "az-raw.parquet")
(out / "source.json").write_text(
    json.dumps(
        dict(
            version=version,
            source=source,
            configuration=cfg,
            rows=d.height,
            molecules=d["molregno"].n_unique(),
            read_only=True,
        ),
        indent=2,
    )
)
print(d.select(pl.len(), pl.col("molregno").n_unique().alias("molecules")))
