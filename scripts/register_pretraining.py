"""Validate completed fold runs and append auxiliary membership, metrics and OOF predictions."""

import argparse
import json

import numpy as np
import polars as pl
from load_molecular_data import copy_rows, verified_run
from psycopg.types.json import Jsonb

from cyp_challenge.database import connect
from cyp_challenge.storage import allocate, digest, now, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", required=True)
    args = parser.parse_args()
    identifier, path = allocate("pretraining-loads")
    result = {"id": identifier, "created_at": now(), "runs": args.run}
    write_json(path / "intent.json", result)
    try:
        with connect() as conn:
            conn.execute("SELECT pg_advisory_xact_lock(847292)")
            counts = []
            for run in args.run:
                output, manifest = verified_run(conn, run)
                summary = json.loads((output / "summary.json").read_text())
                config = summary["config"]
                split_output, _ = verified_run(conn, config["split_run"])
                dataset_output, _ = verified_run(conn, config["dataset_run"])
                input_hashes = {i["sha256"] for i in manifest["inputs"]}
                assert digest(split_output / "assignments.parquet") in input_hashes
                assert digest(dataset_output / "compounds.parquet") in input_hashes
                assert config["checkpoint_sha256"] in input_hashes
                assert digest(output / "best.pt") == summary["best_checkpoint_sha256"]
                history = pl.read_ndjson(output / "history.jsonl")
                assert history["val_mse_z"].min() == summary["best_validation_mse_z"]
                membership = pl.read_parquet(output / "membership.parquet")
                assignments = (
                    pl.read_parquet(split_output / "assignments.parquet")
                    .filter(pl.col("strategy") == "umap")
                    .select("compound_id", "fold")
                )
                joined = membership.join(assignments, on="compound_id", how="left", validate="1:1")
                assert joined["fold"].null_count() == 0
                assert (
                    (joined["fold"] == config["outer_fold"]) == (joined["role"] == "outer_holdout")
                ).all()
                predictions = pl.read_parquet(output / "outer_predictions.parquet")
                assert set(predictions["compound_id"]) == set(
                    membership.filter(pl.col("role") == "outer_holdout")["compound_id"]
                )
                assert set(predictions["fold"]) == {config["outer_fold"]}
                assert np.isfinite(predictions.select("y_true", "prediction").to_numpy()).all()
                if not conn.execute(
                    "SELECT 1 FROM pretraining_sets WHERE run_id=%s", (run,)
                ).fetchone():
                    conn.execute(
                        "INSERT INTO pretraining_sets VALUES (%s,%s,%s,%s,%s,%s)",
                        (
                            run,
                            config["dataset_run"],
                            config["split_run"],
                            config["outer_fold"],
                            config["initialization"],
                            Jsonb(summary),
                        ),
                    )
                    copy_rows(
                        conn,
                        "pretraining_membership",
                        ["run_id", "compound_id", "role"],
                        ((run, *r) for r in membership.select("compound_id", "role").iter_rows()),
                    )
                    columns = ["compound_id", "endpoint", "fold", "y_true", "prediction"]
                    copy_rows(
                        conn,
                        "auxiliary_predictions",
                        ["run_id", *columns],
                        ((run, *r) for r in predictions.select(columns).iter_rows()),
                    )
                    copy_rows(
                        conn,
                        "metrics",
                        ["run_id", "endpoint", "split_id", "fold", "metric", "value"],
                        (
                            (
                                run,
                                row["head"],
                                config["split_run"] + ":umap",
                                str(config["outer_fold"]),
                                "log2fc:" + metric,
                                row[metric],
                            )
                            for row in summary["outer_metrics"]
                            for metric in ["rmse", "mae", "spearman"]
                        ),
                    )
                count = conn.execute(
                    "SELECT count(*) FROM auxiliary_predictions WHERE run_id=%s", (run,)
                ).fetchone()[0]
                assert count == predictions.height
                counts.append({"run": run, "prediction_rows": count})
            result["counts"] = counts
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error=repr(exc))
        raise
    finally:
        write_json(path / "result.json", result)
        print(path)


if __name__ == "__main__":
    main()
