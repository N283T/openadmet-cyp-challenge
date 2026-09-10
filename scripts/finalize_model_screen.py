"""Audit each completed five-fold screen and register validated predictions."""

import argparse
import json
import re
import subprocess
from pathlib import Path

import polars as pl
from load_molecular_data import copy_rows, verified_run
from psycopg.types.json import Jsonb

from cyp_challenge.database import connect
from cyp_challenge.storage import allocate, digest, write_json

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("variant")
    args = parser.parse_args()
    runs = re.findall(
        r"run_id=(\S+)",
        (ROOT / f"artifacts/launches/screen-{args.variant}-20260910.log").read_text(),
    )
    assert len(runs) == len(set(runs)) == 5
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        f"Audit {args.variant} UMAP OOF",
        "--config",
        "configs/readout/audit.json",
        "--input",
        "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
        "--input",
        "artifacts/metric-reference/custom_scoring_functions.py",
    ]
    for run in runs:
        for name in ["oof.parquet", "summary.json"]:
            command += ["--input", str(ROOT / "artifacts/runs/shared" / run / "outputs" / name)]
    result = subprocess.run(
        command + ["--", "python", "scripts/audit_readouts.py"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    print(result.stdout, flush=True)
    audit = re.search(r"run_id=(\S+)", result.stdout)[1]
    _, path = allocate("evaluation-loads")
    record = dict(runs=runs, audit=audit, status="started")
    write_json(path / "intent.json", record)
    try:
        with connect() as conn:
            conn.execute("SELECT pg_advisory_xact_lock(847292)")
            audit_out, audit_manifest = verified_run(conn, audit)
            assert json.loads((audit_out / "audit.json").read_text())["passed"]
            audit_hashes = {i["sha256"] for i in audit_manifest["inputs"]}
            for run in runs:
                out, manifest = verified_run(conn, run)
                assert all(digest(out / f) in audit_hashes for f in ["oof.parquet", "summary.json"])
                summary = json.loads((out / "summary.json").read_text())
                cfg = summary["config"]
                hashes = {i["sha256"] for i in manifest["inputs"]}
                split, _ = verified_run(conn, cfg["split_run"])
                data, _ = verified_run(conn, cfg["dataset_run"])
                assert digest(split / "assignments.parquet") in hashes
                assert digest(data / "measurements.parquet") in hashes
                if args.variant != "direct_ft":
                    feat, _ = verified_run(conn, cfg["feature_run"])
                    assert all(digest(p) in hashes for p in feat.glob("mordred-*.parquet"))
                    if args.variant.startswith("stack_mordred"):
                        parent, _ = verified_run(conn, cfg["stack_run"])
                        assert all(digest(p) in hashes for p in parent.glob("*-features.parquet"))
                else:
                    assert summary["checkpoint_sha256"] in hashes
                    assert digest(data / "compounds.parquet") in hashes
                pred = pl.read_parquet(out / "oof.parquet")
                folds = pl.read_parquet(split / "assignments.parquet").filter(
                    pl.col("strategy") == "umap"
                )
                joined = pred.join(
                    folds.select("compound_id", pl.col("fold").alias("expected")),
                    on="compound_id",
                    how="left",
                    validate="m:1",
                )
                assert (
                    joined["expected"].null_count() == 0
                    and (joined["expected"] == joined["fold"]).all()
                )
                conn.execute(
                    "INSERT INTO model_evaluations VALUES (%s,%s,%s,%s,%s,%s,%s)",
                    (
                        run,
                        audit,
                        cfg["dataset_run"],
                        cfg["split_run"],
                        args.variant,
                        cfg["fold"],
                        Jsonb(summary),
                    ),
                )
                cols = ["compound_id", "endpoint", "fold", "y_true", "prediction"]
                copy_rows(
                    conn,
                    "evaluation_predictions",
                    ["run_id", *cols],
                    ((run, *r) for r in pred.select(cols).iter_rows()),
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
                            args.variant + ":" + m,
                            r[m],
                        )
                        for r in summary["metrics"]
                        for m in ["ST-RAE", "MAE", "Spearman", "MCC", "AUROC", "AP"]
                        if r.get(m) is not None
                    ),
                )
        record["status"] = "completed"
    except Exception as exc:
        record.update(status="failed", error=repr(exc))
        raise
    finally:
        write_json(path / "result.json", record)
        print(path)


if __name__ == "__main__":
    main()
