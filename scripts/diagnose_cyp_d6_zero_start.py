"""Inspect label uncertainty and regime-specific errors of independent D6 models."""

import json
import os
import shutil
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score


def main():
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    audit = pl.read_parquet(source / "source-audit.parquet")
    measurements = pl.read_parquet(source / "measurements.parquet")
    interval = None
    for condition, name in [("direct_inhibition", "direct"), ("TDI_condition", "tdi")]:
        rows = measurements.filter(pl.col("endpoint") == f"CYP2D6_pIC50_{condition}").select(
            "compound_id",
            pl.col("conf_low").alias(f"{name}_low"),
            pl.col("conf_high").alias(f"{name}_high"),
        )
        interval = (
            rows
            if interval is None
            else interval.join(rows, on="compound_id", how="full", coalesce=True)
        )
    audit = audit.join(interval, on="compound_id", how="left").with_columns(
        (pl.col("tdi_low") - pl.col("direct_high")).alias("shift_low"),
        (pl.col("tdi_high") - pl.col("direct_low")).alias("shift_high"),
    )
    # These are conservative combinations of marginal credible intervals, not
    # posterior probabilities of a paired shift (joint samples are unavailable).
    audit = audit.with_columns(
        pl.when(
            ~pl.col("shift_low").is_finite().fill_null(False)
            | ~pl.col("shift_high").is_finite().fill_null(False)
        )
        .then(pl.lit("missing_interval"))
        .when(pl.col("shift_low") > np.log10(2))
        .then(pl.lit("shift_positive_interval"))
        .when(pl.col("shift_high") <= np.log10(2))
        .then(pl.lit("shift_negative_interval"))
        .otherwise(pl.lit("shift_boundary_overlap"))
        .alias("interval_status")
    )
    audit.write_parquet(out / "interval-audit.parquet")
    labelled = audit.filter(pl.col("label").is_finite())
    counts = labelled.group_by("regime", "interval_status").agg(
        pl.len().alias("n"), pl.col("label").sum().alias("positive")
    )
    counts.write_csv(out / "regime-interval-counts.csv")
    predictions = pl.read_parquet(source / "oof.parquet").join(
        audit.select("compound_id", "interval_status"), on="compound_id", how="left"
    )
    records = []
    for split in ["random", "scaffold"]:
        for model in cfg.get(
            "models", ["morgan_logistic", "morgan_rdkit_classifier", "paired_measurement_margin"]
        ):
            base = predictions.filter(
                (pl.col("split_id") == split)
                & (pl.col("model") == model)
                & pl.col("label").is_finite()
            )
            for field in ["regime", "interval_status"]:
                for group in base[field].unique().sort():
                    frame = base.filter(pl.col(field) == group)
                    y, score, calls = (
                        frame["label"].to_numpy(),
                        frame["score"].to_numpy(),
                        frame["call"].to_numpy(),
                    )
                    two = len(np.unique(y)) == 2
                    records.append(
                        dict(
                            split_id=split,
                            model=model,
                            grouping=field,
                            group=group,
                            n=len(frame),
                            positives=int(y.sum()),
                            predicted_positive=int(calls.sum()),
                            recall=float(calls[y == 1].mean()) if (y == 1).any() else None,
                            false_positive_rate=float(calls[y == 0].mean())
                            if (y == 0).any()
                            else None,
                            ap=float(average_precision_score(y, score)) if two else None,
                            auc=float(roc_auc_score(y, score)) if two else None,
                            mcc=float(matthews_corrcoef(y, calls)) if two else None,
                        )
                    )
    pl.DataFrame(records, infer_schema_length=None).write_csv(out / "subgroup-metrics.csv")
    summary = json.loads((source / "summary.json").read_text())
    paired = labelled.filter(pl.col("regime") == "measured_direct_active")
    interval_count = paired.group_by("interval_status").len().sort("interval_status")
    text = [
        "# Independent D6 restart: pilot result and diagnosis",
        "",
        f"Source model run: `classification/{cfg['source_run']}`.",
        "",
        "Every fitted model starts from source molecular structure or paired assay measurements. "
        "No historical ensemble, scores, call count, or prior selected weights are inputs.",
        "",
        "| Family | Random AP | Random nested MCC | Scaffold AP | Scaffold nested MCC |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for model in ["morgan_logistic", "morgan_rdkit_classifier", "paired_measurement_margin"]:
        rows = {r["split_id"]: r for r in summary["aggregate"] if r["model"] == model}
        text.append(
            f"| {model} | {rows['random']['ap']:.4f} | "
            f"{rows['random']['bootstrap_mcc']:.4f} | {rows['scaffold']['ap']:.4f} | "
            f"{rows['scaffold']['bootstrap_mcc']:.4f} |"
        )
    text += [
        "",
        "MCC uses outer calls with inner-selected thresholds. It is not the "
        "optimistic maximum MCC after fitting a threshold to all pooled OOF labels.",
        "",
        "## Measurement uncertainty",
        "",
        "Among the 1,364 supplied labels with direct pIC50 above 4:",
        "",
    ]
    for row in interval_count.iter_rows(named=True):
        text.append(f"- {row['interval_status']}: {row['len']} molecules.")
    text += [
        "",
        "Shift intervals are formed as TDI lower minus direct upper, and TDI upper "
        "minus direct lower. They are conservative combinations of marginal intervals, "
        "not joint posterior shift probabilities. Do not discard or relabel official labels "
        "from this diagnostic. Inferred positives need separate interpretation because their "
        "direct-arm potency is below the assay's reliable range.",
        "",
        "## Operating-point transfer",
        "",
        "| Family | Validation split | Selected threshold | Blind positive calls |",
        "| --- | --- | ---: | ---: |",
    ]
    for row in summary["operating_points"]:
        text.append(
            f"| {row['model']} | {row['split_id']} | "
            f"{row['threshold']:.4f} | {row['blind_positive']} |"
        )
    text += [
        "",
        "These are unsubmitted diagnostic predictions. Counts emerge from the model and "
        "threshold; none was targeted. Their large split sensitivity is evidence against "
        "treating one OOF-optimal point threshold as established.",
        "",
        "## Decision",
        "",
        "Continue the independent restart. Logistic molecular features provide a weak "
        "ranking signal; the point-estimate paired-measurement recipe is a negative result. "
        "Neither justifies adopting a new submission. This does not authorize restoring "
        "the historical model or returning to a fixed call budget.",
        "",
        "Next hypotheses must address paired-shift uncertainty, separation of conventional "
        "and inferred positives, and acquisition-population transfer. More trees or another "
        "small mixture with historical scores is not the next experiment.",
        "",
    ]
    (out / "report.md").write_text("\n".join(text))
    for name in [
        "summary.json",
        "metrics.csv",
        "operating-points.csv",
        "oof.parquet",
        "assignments.parquet",
        "source-audit.parquet",
        "representation.json",
    ]:
        shutil.copy2(source / name, out / name)
    print(interval_count)
    print(pl.DataFrame(summary["aggregate"]))


if __name__ == "__main__":
    main()
