"""Validate hashes/lineage and append split/OOF/metrics to the local registry."""

import argparse
import json

import polars as pl
from load_molecular_data import copy_rows, verified_run
from psycopg.types.json import Jsonb

from cyp_challenge.database import connect
from cyp_challenge.storage import allocate, digest, now, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-run", required=True)
    parser.add_argument("--audit-run", required=True)
    args = parser.parse_args()
    identifier, path = allocate("baseline-loads")
    result = {"id": identifier, "created_at": now(), **vars(args)}
    write_json(path / "intent.json", result)
    try:
        with connect() as conn:
            conn.execute("SELECT pg_advisory_xact_lock(847292)")
            output, manifest = verified_run(conn, args.baseline_run)
            summary = json.loads((output / "summary.json").read_text())
            audit_output, audit_manifest = verified_run(conn, args.audit_run)
            audit_hashes = {item["sha256"] for item in audit_manifest["inputs"]}
            assert all(
                digest(output / f) in audit_hashes
                for f in ["oof.parquet", "metrics.parquet", "summary.json"]
            )
            assert json.loads((audit_output / "audit.json").read_text())["passed"] is True
            summary["audit_run"] = args.audit_run
            config = summary["config"]
            split_id = config["split_run"]
            split_output, split_manifest = verified_run(conn, split_id)
            dataset_output, _ = verified_run(conn, config["dataset_run"])
            feature_output, _ = verified_run(conn, config["feature_run"])
            for parent, files, child in [
                (split_output, ["assignments.parquet", "fingerprints.npz"], manifest),
                (dataset_output, ["measurements.parquet", "aliases.parquet"], manifest),
                (dataset_output, ["compounds.parquet", "aliases.parquet"], split_manifest),
                (
                    feature_output,
                    [p.name for p in feature_output.glob("mordred-*.parquet")],
                    manifest,
                ),
            ]:
                hashes = {item["sha256"] for item in child["inputs"]}
                assert all(digest(parent / file) in hashes for file in files)
            assignments = pl.read_parquet(split_output / "assignments.parquet")
            if not conn.execute("SELECT 1 FROM split_sets WHERE run_id=%s", (split_id,)).fetchone():
                conn.execute(
                    "INSERT INTO split_sets VALUES (%s,%s,%s)",
                    (
                        split_id,
                        config["dataset_run"],
                        Jsonb(json.loads((split_output / "summary.json").read_text())),
                    ),
                )
                copy_rows(
                    conn,
                    "split_assignments",
                    ["split_run_id", "strategy", "compound_id", "fold", "group_label"],
                    (
                        (split_id, *r)
                        for r in assignments.select(
                            "strategy", "compound_id", "fold", "group"
                        ).iter_rows()
                    ),
                )
            existing = conn.execute(
                "SELECT summary FROM baseline_sets WHERE run_id=%s", (args.baseline_run,)
            ).fetchone()
            if existing and existing[0].get("audit_run") != args.audit_run:
                raise ValueError("Baseline already registered with a different immutable audit")
            if not existing:
                oof = pl.read_parquet(output / "oof.parquet")
                joined = oof.join(
                    assignments.select(
                        "strategy", "compound_id", pl.col("fold").alias("expected_fold")
                    ),
                    on=["strategy", "compound_id"],
                    how="left",
                    validate="m:1",
                )
                assert joined["expected_fold"].null_count() == 0
                assert (joined["fold"] == joined["expected_fold"]).all()
                assert oof["split_id"].unique().to_list() == [split_id]
                conn.execute(
                    "INSERT INTO baseline_sets VALUES (%s,%s,%s,%s)",
                    (args.baseline_run, split_id, config["feature_run"], Jsonb(summary)),
                )
                cols = [
                    "strategy",
                    "features",
                    "compound_id",
                    "endpoint",
                    "fold",
                    "y_true",
                    "prediction",
                ]
                copy_rows(
                    conn,
                    "baseline_predictions",
                    ["run_id", *cols],
                    ((args.baseline_run, *r) for r in oof.select(cols).iter_rows()),
                )
                rows = []
                for row in pl.read_parquet(audit_output / "metrics.parquet").iter_rows(named=True):
                    for metric in ["ST-RAE", "MAE", "Spearman", "MCC", "AUROC", "AP"]:
                        if row.get(metric) is not None:
                            rows.append(
                                (
                                    args.audit_run,
                                    row["endpoint"],
                                    split_id + ":" + row["strategy"],
                                    row["fold"],
                                    row["features"] + ":" + metric,
                                    row[metric],
                                )
                            )
                copy_rows(
                    conn,
                    "metrics",
                    ["run_id", "endpoint", "split_id", "fold", "metric", "value"],
                    rows,
                )
            result["split_rows"] = conn.execute(
                "SELECT count(*) FROM split_assignments WHERE split_run_id=%s", (split_id,)
            ).fetchone()[0]
            result["prediction_rows"] = conn.execute(
                "SELECT count(*) FROM baseline_predictions WHERE run_id=%s", (args.baseline_run,)
            ).fetchone()[0]
            result["metric_rows"] = conn.execute(
                "SELECT count(*) FROM metrics WHERE run_id=%s", (args.audit_run,)
            ).fetchone()[0]
            assert result["split_rows"] == assignments.height
            assert result["prediction_rows"] == pl.read_parquet(output / "oof.parquet").height
            assert (
                result["metric_rows"]
                == pl.read_parquet(audit_output / "metrics.parquet").height * 3
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
