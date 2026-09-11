"""Verify curation input hashes and append external target lineage atomically to the DB."""

import json
import os
from pathlib import Path

import polars as pl
import psycopg
from load_molecular_data import copy_rows
from psycopg.types.json import Jsonb

from cyp_challenge.external_training import snapshot_input
from cyp_challenge.storage import digest, write_json


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    dataset = cfg["curation_run"]

    def read(name):
        return pl.read_parquet(snapshot_input(inp, name + ".parquet"))

    compounds, targets = read("external-compounds"), read("external-targets")
    summary = json.loads(snapshot_input(inp, "summary.json").read_text())
    assert compounds.height == summary["compounds"] and targets.height == summary["labels"]
    sources = []
    for kind in ["chembl", "pubchem"]:
        f = read(kind + "-curated-records").filter(
            pl.col("connectivity").is_not_null() & ~pl.col("challenge_overlap")
        )
        if kind == "chembl":
            f = f.with_columns((pl.lit("chembl_ic50_") + pl.col("isoform")).alias("endpoint"))
            idcol, valcol = "activity_id", "pchembl_value"
            meta = [
                "assay_chembl_id",
                "document_chembl_id",
                "doi",
                "pubmed_id",
                "standard_type",
                "standard_relation",
                "standard_value",
                "standard_units",
            ]
        else:
            f = f.filter(pl.col("value").is_finite()).with_columns(
                (pl.lit("pubchem_max_response_") + pl.col("isoform")).alias("endpoint")
            )
            idcol, valcol = "PUBCHEM_SID", "raw_value"
            meta = ["PUBCHEM_CID", "Panel Name", "PUBCHEM_ACTIVITY_OUTCOME", "value"]
        for r in f.select("compound_id", "endpoint", idcol, valcol, *meta).iter_rows(named=True):
            sources.append(
                (
                    dataset,
                    r["compound_id"],
                    r["endpoint"],
                    str(r[idcol]),
                    r[valcol],
                    Jsonb({k: r[k] for k in meta}),
                )
            )
    with psycopg.connect(**cfg["database"]) as conn:
        conn.execute("SELECT pg_advisory_xact_lock(847293)")
        parent = conn.execute(
            "SELECT status, manifest FROM runs WHERE id=%s", (dataset,)
        ).fetchone()
        assert parent and parent[0] == "completed"
        hashes = {item["sha256"] for item in parent[1]["files"]}
        assert all(digest(p) in hashes for p in inp.iterdir())
        assert not conn.execute(
            "SELECT 1 FROM external_datasets WHERE run_id=%s", (dataset,)
        ).fetchone()
        conn.execute(
            "INSERT INTO external_datasets VALUES (%s,%s,%s,%s,now())",
            (
                dataset,
                summary["config"]["chembl_run"],
                summary["config"]["pubchem_run"],
                Jsonb(summary),
            ),
        )
        with conn.cursor() as cursor:
            cursor.executemany(
                "INSERT INTO compounds(id,standardized_smiles) VALUES (%s,%s) "
                "ON CONFLICT (id) DO NOTHING",
                compounds.select("compound_id", "standardized_smiles").iter_rows(),
            )
        copy_rows(
            conn,
            "external_compounds",
            ["dataset_id", "compound_id", "connectivity"],
            ((dataset, *r) for r in compounds.select("compound_id", "connectivity").iter_rows()),
        )
        columns = ["compound_id", "endpoint", "source", "value", "source_rows", "assays"]
        copy_rows(
            conn,
            "external_targets",
            ["dataset_id", *columns],
            ((dataset, *r) for r in targets.select(columns).iter_rows()),
        )
        copy_rows(
            conn,
            "external_target_sources",
            ["dataset_id", "compound_id", "endpoint", "source_record_id", "raw_value", "metadata"],
            sources,
        )
        count = conn.execute(
            "SELECT count(*) FROM external_targets WHERE dataset_id=%s", (dataset,)
        ).fetchone()[0]
        assert count == targets.height
    write_json(
        out / "registered.json",
        dict(
            dataset=dataset,
            compounds=compounds.height,
            targets=targets.height,
            source_records=len(sources),
        ),
    )


if __name__ == "__main__":
    main()
