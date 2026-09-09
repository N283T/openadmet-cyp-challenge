"""Recompute OOF scores against pinned official ST-RAE and report macro comparisons."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

from cyp_challenge.storage import write_json


def main():
    inputs, output = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])

    def read(name):
        return pl.read_parquet(next(inputs.glob("*-" + name + ".parquet")))

    scorer_path = next(inputs.glob("*-custom_scoring_functions.py"))
    spec = importlib.util.spec_from_file_location("official_scoring", scorer_path)
    official = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(official)
    oof, scores, labels = read("oof"), read("metrics"), read("measurements")
    assert (
        oof.select("strategy", "features", "compound_id", "endpoint").unique().height == oof.height
    )
    joined = oof.join(
        labels.select("compound_id", "endpoint", "value", "conf_low", "conf_high"),
        on=["compound_id", "endpoint"],
        how="left",
        validate="m:1",
    )
    assert joined["value"].null_count() == 0
    assert (joined["y_true"] == joined["value"]).all()
    summary = json.loads(next(inputs.glob("*-summary.json")).read_text())
    rebuilt = []
    missing_metric_columns = sorted(
        {"ST-RAE", "MAE", "Spearman", "MCC", "AUROC", "AP"} - set(scores.columns)
    )
    for key, frame in joined.group_by("endpoint", "strategy", "features"):
        endpoint, strategy, features = key
        expected_ids = set(labels.filter(pl.col("endpoint") == endpoint)["compound_id"])
        assert set(frame["compound_id"]) == expected_ids
        for fold in [-1, 0, 1, 2, 3, 4]:
            f = frame if fold == -1 else frame.filter(pl.col("fold") == fold)
            y, p = f["y_true"].to_numpy(), f["prediction"].to_numpy()
            if endpoint.endswith("is_TDI"):
                metric = "MCC"
                result = matthews_corrcoef(y, p >= 0.5)
                extra = {
                    "AUROC": float(roc_auc_score(y, p)),
                    "AP": float(average_precision_score(y, p)),
                }
            else:
                metric = "ST-RAE"
                result = official.rae_soft_threshold_absolute_error(
                    y,
                    p,
                    y_true_upper=f["conf_high"].to_numpy(),
                    y_true_lower=f["conf_low"].to_numpy(),
                )
                extra = {
                    "MAE": float(np.abs(y - p).mean()),
                    "Spearman": float(spearmanr(y, p).statistic),
                }
            rebuilt.append(
                {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "features": features,
                    "fold": fold,
                    "n": len(y),
                    metric: float(result),
                    **extra,
                }
            )
            original = scores.filter(
                (pl.col("endpoint") == endpoint)
                & (pl.col("strategy") == strategy)
                & (pl.col("features") == features)
                & (pl.col("fold") == fold)
            )
            assert original.height == 1
            if metric in original.columns:
                assert np.isclose(original[metric].item(), result, rtol=1e-12, atol=1e-12)
            if fold == -1:
                saved = [
                    r
                    for r in summary["pooled_metrics"]
                    if r["endpoint"] == endpoint
                    and r["strategy"] == strategy
                    and r["features"] == features
                ]
                assert len(saved) == 1
                assert np.isclose(saved[0][metric], result, rtol=1e-12, atol=1e-12)
    scores = pl.DataFrame(rebuilt, infer_schema_length=None)
    scores.write_parquet(output / "metrics.parquet")
    macro = (
        scores.filter(pl.col("fold") == -1)
        .group_by("strategy", "features")
        .agg(
            pl.col("ST-RAE").mean(),
            pl.col("MAE").mean(),
            pl.col("MCC").mean(),
            pl.col("AUROC").mean(),
        )
        .sort("strategy", "features")
    )
    diagnostics = read("diagnostics")
    proximity = (
        diagnostics.group_by("endpoint", "strategy")
        .agg(
            ((pl.col("nn_mean") * pl.col("n")).sum() / pl.col("n").sum()).alias("nn_mean"),
            ((pl.col("nn_gt_05") * pl.col("n")).sum() / pl.col("n").sum()).alias("nn_gt_05"),
            pl.col("n").min().alias("min_fold_n"),
            pl.col("n").max().alias("max_fold_n"),
        )
        .sort("endpoint", "strategy")
    )
    macro.write_csv(output / "macro.csv")
    proximity.write_csv(output / "proximity.csv")
    scores.filter(pl.col("fold") == -1).write_csv(output / "endpoint_scores.csv")
    fold_ranges = (
        scores.filter(pl.col("fold") >= 0)
        .group_by("endpoint", "strategy", "features")
        .agg(
            pl.col("ST-RAE").min().alias("ST-RAE_min"),
            pl.col("ST-RAE").max().alias("ST-RAE_max"),
            pl.col("MCC").min().alias("MCC_min"),
            pl.col("MCC").max().alias("MCC_max"),
        )
    )
    fold_ranges.write_csv(output / "fold_ranges.csv")
    write_json(
        output / "audit.json",
        {
            "passed": True,
            "original_missing_metric_columns": missing_metric_columns,
            "oof_rows": oof.height,
            "scored_groups": scores.height,
            "official_revision": next(inputs.glob("*-revision.txt")).read_text().strip(),
            "macro": macro.to_dicts(),
        },
    )
    lines = [
        "# LightGBM split comparison",
        "",
        "Pooled OOF, endpoint macro averages; ST-RAE lower is better, MCC higher is better.",
        "",
        "| Split | Features | ST-RAE | MAE | MCC | AUROC |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in macro.iter_rows(named=True):
        lines.append(
            "| "
            + " | ".join(
                [row["strategy"], row["features"]]
                + [f"{row[m]:.4f}" for m in ["ST-RAE", "MAE", "MCC", "AUROC"]]
            )
            + " |"
        )
    lines += [
        "",
        "## Mean nearest-neighbor Tanimoto",
        "",
        "| Endpoint | Split | Mean max Tanimoto | Fraction >0.5 | Min fold n | Max fold n |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in proximity.iter_rows(named=True):
        lines.append(
            f"| {row['endpoint']} | {row['strategy']} | {row['nn_mean']:.4f} | "
            f"{row['nn_gt_05']:.4f} | {row['min_fold_n']} | {row['max_fold_n']} |"
        )
    lines += [
        "",
        "One seed, one learner, no tuning; no claim of uncertainty or leaderboard calibration.",
        "Validation uses 80% reference labels; the test comparison uses all endpoint labels.",
    ]
    (output / "report.md").write_text("\n".join(lines) + "\n")
    print(macro)
    print(proximity)


if __name__ == "__main__":
    main()
