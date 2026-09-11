"""Record live local ChEMBL metadata and grouped evidence as immutable run outputs."""

import datetime as dt
import json
import os
from pathlib import Path

import polars as pl
import psycopg
from psycopg.rows import dict_row

from cyp_challenge.storage import write_json


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    query = next(inp.glob("*-" + cfg["sql_input"])).read_text()
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
        targets = conn.execute(
            "SELECT tid, chembl_id, pref_name, organism, target_type FROM target_dictionary "
            "WHERE organism = 'Homo sapiens' AND pref_name ILIKE 'Cytochrome P450%' "
            "ORDER BY pref_name"
        ).fetchall()
        rows = conn.execute(query).fetchall()
    table = pl.DataFrame(rows, infer_schema_length=None)
    assert table.height > 0 and table["target_chembl_id"].n_unique() == 4
    table.write_parquet(out / "activity-inventory.parquet")
    table.write_csv(out / "activity-inventory.csv")
    pl.DataFrame(targets).write_csv(out / "human-cyp-targets.csv")
    # These are strata sums: molecule/assay counts must not be summed across strata.
    summary = (
        table.group_by("pref_name", "standard_type")
        .agg(
            pl.col("rows").sum(),
            pl.col("with_pchembl").sum(),
        )
        .sort("pref_name", "rows", descending=[False, True])
    )
    summary.write_csv(out / "measurement-types.csv")
    (out / "executed-query.sql").write_text(query)
    write_json(
        out / "source.json",
        dict(
            config=cfg,
            database_version=version,
            retrieved_at=dt.datetime.now(dt.UTC).isoformat(),
            transaction="repeatable read, read only",
            grouped_rows=table.height,
            caveat="Live DB query/version and result snapshot retained; not curated training data.",
        ),
    )
    print(summary.filter(pl.col("standard_type").is_in(["IC50", "Ki", "Inhibition"])))


if __name__ == "__main__":
    main()
