"""Export uncurated CYP measurements with source identity from the existing local DB."""

import json
import os
from pathlib import Path

import polars as pl
import psycopg
from psycopg.rows import dict_row

from cyp_challenge.storage import now, write_json


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    query = next(inp.glob("*-extract_chembl.sql")).read_text()
    with psycopg.connect(
        dbname=cfg["dbname"],
        host=cfg["host"],
        port=cfg["port"],
        row_factory=dict_row,
        options="-c default_transaction_read_only=on -c statement_timeout=300000",
    ) as conn:
        conn.execute("SET TRANSACTION ISOLATION LEVEL REPEATABLE READ READ ONLY")
        version = conn.execute(
            "SELECT name, creation_date::text AS creation_date, comments FROM version"
        ).fetchall()
        sources = conn.execute("SELECT * FROM source ORDER BY src_id").fetchall()
        rows = conn.execute(query).fetchall()
    df = pl.DataFrame(rows, infer_schema_length=None)
    assert df.height == df["activity_id"].n_unique()
    df.write_parquet(out / "chembl-raw.parquet")
    pl.DataFrame(sources).write_csv(out / "chembl-sources.csv")
    df.select(
        "assay_chembl_id",
        "target_name",
        "assay_description",
        "assay_type",
        "confidence_score",
        "assay_src_id",
        "src_assay_id",
        "document_chembl_id",
        "pubmed_id",
        "doi",
        "title",
    ).unique().write_parquet(out / "assay-documents.parquet")
    df.group_by("target_name", "standard_type", "assay_src_id").agg(
        pl.len().alias("rows"),
        pl.col("molregno").n_unique().alias("molecules"),
    ).sort("target_name", "standard_type", "assay_src_id").write_csv(out / "counts.csv")
    write_json(
        out / "source.json",
        dict(
            config=cfg,
            version=version,
            retrieved_at=now(),
            rows=df.height,
            molecules=df["molregno"].n_unique(),
            read_only=True,
        ),
    )
    print(f"Exported {df.height} activity rows, {df['molregno'].n_unique()} molecules", flush=True)


if __name__ == "__main__":
    main()
