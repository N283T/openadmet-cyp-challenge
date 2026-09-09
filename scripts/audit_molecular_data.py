"""Persist DB/Parquet reconciliation for one immutable dataset and feature set."""

import argparse
import json
from collections import Counter
from pathlib import Path

import polars as pl
from psycopg import sql

from cyp_challenge.database import connect
from cyp_challenge.storage import allocate, now, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-run", required=True)
    parser.add_argument("--feature-run", required=True)
    args = parser.parse_args()
    identifier, path = allocate("data-audits")
    report = {"id": identifier, "created_at": now(), **vars(args), "checks": {}, "families": {}}
    write_json(path / "intent.json", report)
    try:
        with connect() as conn:
            dataset = conn.execute(
                "SELECT summary FROM molecular_datasets WHERE run_id=%s", (args.dataset_run,)
            ).fetchone()[0]
            feature = conn.execute(
                "SELECT catalog,summary,dataset_id FROM feature_sets WHERE run_id=%s",
                (args.feature_run,),
            ).fetchone()
            catalog, summary, dataset_id = feature
            report["checks"]["dataset_lineage"] = dataset_id == args.dataset_run
            report["checks"]["expected_compounds"] = (
                summary["compounds"] == dataset["unique_compounds"]
            )
            count = conn.execute(
                "SELECT count(*) FROM source_records WHERE dataset_id=%s AND role='test'",
                (args.dataset_run,),
            ).fetchone()[0]
            report["checks"]["all_750_test_identities_retained"] = count == 750
            source_count = conn.execute(
                "SELECT count(*) FROM record_origins WHERE dataset_id=%s", (args.dataset_run,)
            ).fetchone()[0]
            report["checks"]["all_source_rows_traceable"] = source_count == sum(
                file["original_rows"] for file in dataset["files"].values()
            )
            output = (
                Path(
                    conn.execute(
                        "SELECT manifest_path FROM runs WHERE id=%s", (args.feature_run,)
                    ).fetchone()[0]
                ).parent
                / "outputs"
            )
            for family, table in [
                ("rdkit", "compound_descriptors_full"),
                ("mordred", "compound_mordred"),
            ]:
                expected_names = set(catalog["names"][family])
                missing = Counter()
                rows = 0
                cursor_name = f"audit_{family}"
                with conn.cursor(name=cursor_name) as cursor:
                    cursor.execute(
                        sql.SQL(
                            "SELECT compound_id,descriptors FROM {} WHERE feature_set_id=%s"
                        ).format(sql.Identifier(table)),
                        (args.feature_run,),
                    )
                    for cid, values in cursor:
                        if set(values) != expected_names:
                            raise ValueError(f"Missing descriptor keys for {cid}")
                        rows += 1
                        missing.update(name for name, value in values.items() if value is None)
                report["checks"][f"{family}_row_count"] = rows == summary["compounds"]
                report["checks"][f"{family}_missing_counts"] = all(
                    missing[name] == summary["missing_counts"][family][name]
                    for name in expected_names
                )
                # Retain all constant/all-null columns; report them without selecting features.
                lazy = pl.scan_parquet(str(output / f"{family}-*.parquet"))
                unique = (
                    lazy.select(pl.exclude("compound_id").drop_nulls().n_unique())
                    .collect()
                    .row(0, named=True)
                )
                report["families"][family] = {
                    "rows": rows,
                    "columns": len(expected_names),
                    "null_values": sum(missing.values()),
                    "all_null_columns": sorted(
                        name for name in expected_names if missing[name] == rows
                    ),
                    "constant_nonnull_columns": sorted(
                        name for name, count in unique.items() if count == 1
                    ),
                }
            error_count = conn.execute(
                "SELECT count(*) FROM descriptor_errors WHERE feature_set_id=%s",
                (args.feature_run,),
            ).fetchone()[0]
            report["checks"]["error_records_cover_all_compounds"] = (
                error_count == summary["compounds"]
            )
        if not all(report["checks"].values()):
            raise ValueError("DB/Parquet reconciliation failed")
        report["status"] = "passed"
    except Exception as exc:
        report.update(status="failed", error=repr(exc))
        raise
    finally:
        write_json(path / "report.json", report)
        print(path)
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
