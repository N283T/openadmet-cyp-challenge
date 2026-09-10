"""Verify a completed frozen-feature run and append its vectors to PostgreSQL."""

import argparse
import json

import polars as pl
from load_molecular_data import copy_rows, verified_run
from psycopg.types.json import Jsonb

from cyp_challenge.database import connect
from cyp_challenge.embedding_validation import validate_embeddings
from cyp_challenge.storage import allocate, digest, now, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", required=True)
    args = parser.parse_args()
    identifier, path = allocate("embedding-loads")
    result = {"id": identifier, "created_at": now(), "run": args.run}
    write_json(path / "intent.json", result)
    try:
        with connect() as conn:
            conn.execute("SELECT pg_advisory_xact_lock(847292)")
            output, manifest = verified_run(conn, args.run)
            summary = json.loads((output / "summary.json").read_text())
            provenance_file = next((output.parent / "inputs").glob("*-provenance.json"))
            provenance = json.loads(provenance_file.read_text())
            assert provenance["checkpoint_sha256"] == summary["checkpoint_sha256"]
            assert digest(provenance_file) in {item["sha256"] for item in manifest["inputs"]}
            summary["checkpoint_provenance"] = provenance
            dataset_id = summary["config"]["dataset_run"]
            dataset, _ = verified_run(conn, dataset_id)
            hashes = {item["sha256"] for item in manifest["inputs"]}
            assert digest(dataset / "compounds.parquet") in hashes
            assert summary["checkpoint_sha256"] in hashes
            assert digest(output / "embeddings.parquet") == summary["output_sha256"]
            frame = pl.read_parquet(output / "embeddings.parquet")
            ids = pl.read_parquet(dataset / "compounds.parquet")["compound_id"].to_list()
            validate_embeddings(frame, ids, summary["columns"])
            if not conn.execute(
                "SELECT 1 FROM embedding_sets WHERE run_id=%s", (args.run,)
            ).fetchone():
                conn.execute(
                    "INSERT INTO embedding_sets VALUES (%s,%s,%s,%s,%s)",
                    (args.run, dataset_id, "chemeleon_raw", summary["dimension"], Jsonb(summary)),
                )
                copy_rows(
                    conn,
                    "compound_embeddings",
                    ["feature_set_id", "compound_id", "dimension", "embedding"],
                    (
                        (args.run, row[0], summary["dimension"], list(row[1:]))
                        for row in frame.select("compound_id", *summary["columns"]).iter_rows()
                    ),
                )
            rows = conn.execute(
                "SELECT count(*) FROM compound_embeddings WHERE feature_set_id=%s", (args.run,)
            ).fetchone()[0]
            assert rows == len(ids)
            result.update(rows=rows, dimension=summary["dimension"])
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error=repr(exc))
        raise
    finally:
        write_json(path / "result.json", result)
        print(path)


if __name__ == "__main__":
    main()
