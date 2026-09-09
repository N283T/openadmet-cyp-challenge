"""Verify and atomically append completed preprocessing/descriptor runs to PostgreSQL."""

import argparse
import json
from pathlib import Path

import polars as pl
from psycopg import sql
from psycopg.types.json import Jsonb

from cyp_challenge.database import connect
from cyp_challenge.storage import allocate, digest, now, write_json


def verified_run(conn, identifier):
    row = conn.execute(
        "SELECT status,manifest_path,manifest FROM runs WHERE id=%s", (identifier,)
    ).fetchone()
    if not row or row[0] != "completed":
        raise ValueError(f"Expected a completed run: {identifier}")
    path = Path(row[1]).parent
    for file in row[2]["files"]:
        if digest(path / file["path"]) != file["sha256"]:
            raise ValueError(f"Run artifact changed: {file['path']}")
    return path / "outputs", row[2]


def copy_rows(conn, table, columns, rows):
    query = sql.SQL("COPY {} ({}) FROM STDIN").format(
        sql.Identifier(table), sql.SQL(",").join(map(sql.Identifier, columns))
    )
    with conn.cursor().copy(query) as copy:
        for row in rows:
            copy.write_row(row)


def load_dataset(conn, identifier):
    output, _ = verified_run(conn, identifier)
    if conn.execute("SELECT 1 FROM molecular_datasets WHERE run_id=%s", (identifier,)).fetchone():
        print("Dataset already loaded; no changes")
        return
    summary = json.loads((output / "summary.json").read_text())
    compounds = pl.read_parquet(output / "compounds.parquet")
    conn.execute(
        "INSERT INTO molecular_datasets(run_id,summary) VALUES (%s,%s)",
        (identifier, Jsonb(summary)),
    )
    # Existing identities are shared across dataset versions, never overwritten.
    with conn.cursor() as cursor:
        cursor.executemany(
            "INSERT INTO compounds(id,standardized_smiles) VALUES (%s,%s) "
            "ON CONFLICT (id) DO NOTHING",
            compounds.iter_rows(),
        )
    copy_rows(
        conn,
        "dataset_compounds",
        ["dataset_id", "compound_id"],
        ((identifier, cid) for cid in compounds["compound_id"]),
    )
    columns = [
        "record_id",
        "source_file",
        "role",
        "molecule_name",
        "raw_smiles",
        "compound_id",
        "payload",
    ]
    with (output / "source_records.jsonl").open() as stream:

        def records():
            for line in stream:
                record = json.loads(line)
                yield (
                    identifier,
                    *(Jsonb(record[k]) if k == "payload" else record[k] for k in columns),
                )

        copy_rows(conn, "source_records", ["dataset_id", *columns], records())
    for table, filename, columns in [
        ("record_origins", "record_origins", ["source_file", "source_row", "record_id"]),
        (
            "measurements",
            "measurements",
            ["measurement_id", "compound_id", "endpoint", "value", "conf_low", "conf_high", "std"],
        ),
        ("measurement_sources", "measurement_sources", ["measurement_id", "record_id"]),
    ]:
        frame = pl.read_parquet(output / (filename + ".parquet")).select(columns)
        copy_rows(
            conn, table, ["dataset_id", *columns], ((identifier, *r) for r in frame.iter_rows())
        )
    assert (
        conn.execute(
            "SELECT count(*) FROM dataset_compounds WHERE dataset_id=%s", (identifier,)
        ).fetchone()[0]
        == compounds.height
    )
    print(f"Loaded dataset {identifier}: {compounds.height} unique compounds", flush=True)


def load_features(conn, identifier, dataset_id):
    output, manifest = verified_run(conn, identifier)
    dataset_output, _ = verified_run(conn, dataset_id)
    expected = digest(dataset_output / "compounds.parquet")
    if not any(item["sha256"] == expected for item in manifest["inputs"]):
        raise ValueError("Descriptor inputs do not match the requested dataset")
    if conn.execute("SELECT 1 FROM feature_sets WHERE run_id=%s", (identifier,)).fetchone():
        print("Feature set already loaded; no changes")
        return
    catalog = json.loads((output / "descriptor_catalog.json").read_text())
    summary = json.loads((output / "summary.json").read_text())
    conn.execute(
        "INSERT INTO feature_sets(run_id,dataset_id,catalog,summary) VALUES (%s,%s,%s,%s)",
        (identifier, dataset_id, Jsonb(catalog), Jsonb(summary)),
    )
    expected_ids = set(pl.read_parquet(dataset_output / "compounds.parquet")["compound_id"])
    for family, table in [("rdkit", "compound_descriptors_full"), ("mordred", "compound_mordred")]:
        seen = set()
        for file in sorted(output.glob(f"{family}-*.parquet")):
            frame = pl.read_parquet(file)
            if set(frame.columns) != {"compound_id", *catalog["names"][family]}:
                raise ValueError("Descriptor columns differ from catalog")
            ids = frame["compound_id"].to_list()
            if seen.intersection(ids) or len(ids) != len(set(ids)):
                raise ValueError("Duplicate descriptor compound")
            seen.update(ids)

            def values(frame=frame):
                for row in frame.iter_rows(named=True):
                    cid = row.pop("compound_id")
                    yield identifier, cid, Jsonb(row)

            copy_rows(conn, table, ["feature_set_id", "compound_id", "descriptors"], values())
        if seen != expected_ids:
            raise ValueError(f"Incomplete/extra descriptor compounds: {family}")
        print(f"Loaded {family}: {len(seen)} compounds", flush=True)
    for file in sorted(output.glob("errors-*.parquet")):
        rows = pl.read_parquet(file).iter_rows()
        copy_rows(
            conn,
            "descriptor_errors",
            ["feature_set_id", "compound_id", "errors"],
            ((identifier, cid, Jsonb(json.loads(errors))) for cid, errors in rows),
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-run", required=True)
    parser.add_argument("--feature-run")
    args = parser.parse_args()
    identifier, path = allocate("data-loads")
    manifest = {"id": identifier, "created_at": now(), **vars(args)}
    write_json(path / "intent.json", manifest)
    try:
        with connect() as conn:
            conn.execute("SELECT pg_advisory_xact_lock(847292)")
            load_dataset(conn, args.dataset_run)
            if args.feature_run:
                load_features(conn, args.feature_run, args.dataset_run)
        manifest["status"] = "completed"
    except Exception as exc:
        manifest.update(status="failed", error=repr(exc))
        raise
    finally:
        write_json(path / "result.json", manifest)
        print(path)


if __name__ == "__main__":
    main()
