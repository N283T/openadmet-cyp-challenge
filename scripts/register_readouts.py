"""Register fold-matched readouts after checking encoder and split lineage."""

import argparse
import json

import polars as pl
from load_molecular_data import copy_rows, verified_run
from psycopg.types.json import Jsonb

from cyp_challenge.database import connect
from cyp_challenge.storage import allocate, digest, now, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", required=True)
    args = parser.parse_args()
    identifier, path = allocate("readout-loads")
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
                encoder, _ = verified_run(conn, config["encoder_run"])
                split, _ = verified_run(conn, config["split_run"])
                dataset, _ = verified_run(conn, config["dataset_run"])
                hashes = {i["sha256"] for i in manifest["inputs"]}
                for parent, files in [
                    (
                        encoder,
                        [
                            "summary.json",
                            "membership.parquet",
                            *[p.name for p in encoder.glob("embeddings-*.parquet")],
                        ],
                    ),
                    (split, ["assignments.parquet"]),
                    (dataset, ["measurements.parquet", "aliases.parquet"]),
                ]:
                    assert all(digest(parent / file) in hashes for file in files)
                assert all(h in hashes for h in summary["checkpoint_sha256"].values())
                encoder_config = json.loads((encoder / "summary.json").read_text())["config"]
                assert encoder_config["outer_fold"] == config["fold"]
                assert encoder_config["split_run"] == config["split_run"]
                predictions = pl.read_parquet(output / "oof.parquet")
                folds = pl.read_parquet(split / "assignments.parquet").filter(
                    pl.col("strategy") == "umap"
                )
                joined = predictions.join(
                    folds.select("compound_id", pl.col("fold").alias("expected_fold")),
                    on="compound_id",
                    how="left",
                    validate="m:1",
                )
                assert joined["expected_fold"].null_count() == 0
                assert (joined["fold"] == joined["expected_fold"]).all()
                assert set(predictions["encoder_run"]) == {config["encoder_run"]}
                if not conn.execute(
                    "SELECT 1 FROM readout_sets WHERE run_id=%s", (run,)
                ).fetchone():
                    conn.execute(
                        "INSERT INTO readout_sets VALUES (%s,%s,%s,%s,%s,%s)",
                        (
                            run,
                            config["dataset_run"],
                            config["encoder_run"],
                            config["split_run"],
                            config["fold"],
                            Jsonb(summary),
                        ),
                    )
                    columns = ["compound_id", "endpoint", "fold", "y_true", "prediction"]
                    copy_rows(
                        conn,
                        "readout_predictions",
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
                                row["endpoint"],
                                config["split_run"] + ":umap",
                                str(config["fold"]),
                                "tabpfn_embedding:" + metric,
                                row[metric],
                            )
                            for row in summary["metrics"]
                            for metric in ["ST-RAE", "MAE", "Spearman", "MCC", "AUROC", "AP"]
                            if row.get(metric) is not None
                        ),
                    )
                count = conn.execute(
                    "SELECT count(*) FROM readout_predictions WHERE run_id=%s", (run,)
                ).fetchone()[0]
                assert count == predictions.height
                counts.append({"run": run, "rows": count})
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
