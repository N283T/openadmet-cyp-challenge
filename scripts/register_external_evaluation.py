"""Append a verified external-data model comparison to existing evaluation tables."""

import json
import os
from pathlib import Path

import polars as pl
import psycopg
from load_molecular_data import copy_rows
from psycopg.types.json import Jsonb

from cyp_challenge.storage import digest, write_json


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    audit = json.loads(next(inp.glob("*-audit.json")).read_text())
    assert audit["passed"]
    audit_id = cfg["audit_run"]
    with psycopg.connect(**cfg["database"]) as conn:
        conn.execute("SELECT pg_advisory_xact_lock(847294)")
        row = conn.execute("SELECT status, manifest FROM runs WHERE id=%s", (audit_id,)).fetchone()
        assert row[0] == "completed"
        audit_hashes = {r["sha256"] for r in row[1]["files"]}
        assert all(digest(p) in audit_hashes for p in inp.iterdir())
        score = pl.read_csv(next(inp.glob("*-metrics.csv")))
        for variant in ["control", "external"]:
            predictions = pl.read_parquet(next(inp.glob(f"*-{variant}-oof.parquet")))
            for entry in audit["config"]["variants"][variant]:
                run_id, fold = entry["run_id"], entry["fold"]
                parent = conn.execute(
                    "SELECT status,manifest_path,manifest FROM runs WHERE id=%s", (run_id,)
                ).fetchone()
                assert parent[0] == "completed"
                summary_path = Path(parent[1]).parent / "outputs/summary.json"
                # Audit input hashes bind the original model summary and predictions.
                audit_input_hashes = {r["sha256"] for r in row[1]["inputs"]}
                assert digest(summary_path) in audit_input_hashes
                assert digest(summary_path.parent / "oof.parquet") in audit_input_hashes
                summary = json.loads(summary_path.read_text())
                model_cfg = summary["config"]
                conn.execute(
                    "INSERT INTO model_evaluations VALUES (%s,%s,%s,%s,%s,%s,%s)",
                    (
                        run_id,
                        audit_id,
                        model_cfg["dataset_run"],
                        model_cfg["split_run"],
                        "external_v1_" + variant,
                        fold,
                        Jsonb(summary),
                    ),
                )
                columns = ["compound_id", "endpoint", "fold", "y_true", "prediction"]
                f = predictions.filter(pl.col("fold") == fold)
                copy_rows(
                    conn,
                    "evaluation_predictions",
                    ["run_id", *columns],
                    ((run_id, *r) for r in f.select(columns).iter_rows()),
                )
                m = score.filter((pl.col("variant") == variant) & (pl.col("fold") == fold))
                copy_rows(
                    conn,
                    "metrics",
                    ["run_id", "endpoint", "split_id", "fold", "metric", "value"],
                    (
                        (
                            run_id,
                            r["endpoint"],
                            model_cfg["split_run"] + ":umap",
                            str(fold),
                            metric,
                            r[metric],
                        )
                        for r in m.iter_rows(named=True)
                        for metric in ["ST_RAE", "MAE", "Spearman"]
                    ),
                )
        for r in audit["pooled"]:
            copy_rows(
                conn,
                "metrics",
                ["run_id", "endpoint", "split_id", "fold", "metric", "value"],
                (
                    (
                        audit_id,
                        r["endpoint"],
                        "20260909T135130268493Z-edba6fc0:umap",
                        "-1",
                        r["variant"] + ":" + metric,
                        r[metric],
                    )
                    for metric in ["ST_RAE", "MAE", "Spearman"]
                ),
            )
    write_json(
        out / "registered.json", dict(audit_run=audit_id, model_evaluations=10, status="registered")
    )


if __name__ == "__main__":
    main()
