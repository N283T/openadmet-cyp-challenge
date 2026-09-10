"""Verify completed stacked readouts and register their predictions and parents."""

import argparse
import json

import polars as pl
from load_molecular_data import copy_rows, verified_run
from psycopg.types.json import Jsonb

from cyp_challenge.database import connect
from cyp_challenge.storage import allocate, digest, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", required=True)
    args = parser.parse_args()
    _, path = allocate("stacking-loads")
    result = {"runs": args.run, "status": "started"}
    write_json(path / "intent.json", result)
    try:
        with connect() as conn:
            conn.execute("SELECT pg_advisory_xact_lock(847292)")
            for run in args.run:
                out, manifest = verified_run(conn, run)
                summary = json.loads((out / "summary.json").read_text())
                cfg = summary["config"]
                hashes = {i["sha256"] for i in manifest["inputs"]}
                for parent in cfg["auxiliary"]:
                    source, _ = verified_run(conn, parent["run"])
                    for p in [
                        source / "summary.json",
                        source / "membership.parquet",
                        *source.glob("predictions-*.parquet"),
                    ]:
                        assert digest(p) in hashes
                foundation, _ = verified_run(conn, cfg["foundation_run"])
                assert all(digest(p) in hashes for p in foundation.glob("chemeleon-*.parquet"))
                split, _ = verified_run(conn, cfg["split_run"])
                assert digest(split / "assignments.parquet") in hashes
                dataset, _ = verified_run(conn, cfg["dataset_run"])
                assert all(
                    digest(dataset / p) in hashes
                    for p in ["measurements.parquet", "aliases.parquet"]
                )
                assert all(h in hashes for h in summary["checkpoint_sha256"].values())
                predictions = pl.read_parquet(out / "oof.parquet")
                folds = pl.read_parquet(split / "assignments.parquet").filter(
                    pl.col("strategy") == "umap"
                )
                expected = predictions.join(
                    folds.select("compound_id", pl.col("fold").alias("expected")),
                    on="compound_id",
                    how="left",
                    validate="m:1",
                )
                assert expected["expected"].null_count() == 0
                assert (expected["fold"] == expected["expected"]).all()
                assert set(predictions["fold"]) == {cfg["fold"]}
                if conn.execute(
                    "SELECT 1 FROM stacked_readouts WHERE run_id=%s", (run,)
                ).fetchone():
                    continue
                conn.execute(
                    "INSERT INTO stacked_readouts VALUES (%s,%s,%s,%s)",
                    (run, cfg["split_run"], cfg["fold"], Jsonb(summary)),
                )
                columns = ["compound_id", "endpoint", "fold", "y_true", "prediction"]
                copy_rows(
                    conn,
                    "stacked_predictions",
                    ["run_id", *columns],
                    ((run, *r) for r in predictions.select(columns).iter_rows()),
                )
                for parent in cfg["auxiliary"]:
                    conn.execute(
                        "INSERT INTO stacked_auxiliary_parents VALUES (%s,%s)", (run, parent["run"])
                    )
                copy_rows(
                    conn,
                    "metrics",
                    ["run_id", "endpoint", "split_id", "fold", "metric", "value"],
                    (
                        (
                            run,
                            r["endpoint"],
                            cfg["split_run"] + ":umap",
                            str(cfg["fold"]),
                            "pca_log2fc_tabpfn:" + m,
                            r[m],
                        )
                        for r in summary["metrics"]
                        for m in ["ST-RAE", "MAE", "Spearman", "MCC", "AUROC", "AP"]
                        if r.get(m) is not None
                    ),
                )
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error=repr(exc))
        raise
    finally:
        write_json(path / "result.json", result)
        print(path)


if __name__ == "__main__":
    main()
