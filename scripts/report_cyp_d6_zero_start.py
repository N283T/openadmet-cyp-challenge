"""Compare independent restart families using paired scaffold-group bootstrap."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score


def weighted_mcc(y, calls, weights):
    tp = weights @ ((y == 1) & calls)
    fp = weights @ ((y == 0) & calls)
    fn = weights @ ((y == 1) & ~calls)
    tn = weights @ ((y == 0) & ~calls)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return np.divide(
        tp * tn - fp * fn, denominator, out=np.zeros_like(denominator), where=denominator > 0
    )


def main():
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    frames, summaries = [], []
    for version, run in cfg["sources"].items():
        frame = pl.read_parquet(source / f"{version}-oof.parquet").with_columns(
            pl.lit(version).alias("version")
        )
        frames.append(frame)
        summary = json.loads((source / f"{version}-summary.json").read_text())
        summaries.append((version, run, summary))
    predictions = pl.concat(frames)
    predictions.write_parquet(out / "all-independent-oof.parquet")
    comparisons = []
    for split in cfg.get("split_strategies", ["random", "scaffold"]):
        assignments = pl.read_parquet(
            source / cfg.get("assignment_aliases", {}).get(split, "assignments.parquet")
        )
        network = predictions.filter(
            (pl.col("split_id") == split)
            & (pl.col("model") == "fresh_multiassay_network")
            & pl.col("label").is_finite()
        ).sort("compound_id")
        groups = network.select("compound_id").join(
            assignments.filter(pl.col("split_id") == split).select("compound_id", "group"),
            on="compound_id",
            how="left",
        )
        unique, codes = np.unique(groups["group"].to_numpy(), return_inverse=True)
        rng = np.random.default_rng(cfg["seed"])
        draws = rng.integers(0, len(unique), (cfg["bootstrap_samples"], len(unique)))
        group_weights = np.zeros((cfg["bootstrap_samples"], len(unique)))
        np.add.at(group_weights, (np.arange(cfg["bootstrap_samples"])[:, None], draws), 1)
        weights = group_weights[:, codes]
        labels = network["label"].to_numpy()
        network_ap = np.array(
            [
                average_precision_score(labels, network["score"].to_numpy(), sample_weight=w)
                for w in weights
            ]
        )
        network_mcc = weighted_mcc(labels, network["call"].to_numpy(), weights)
        baseline_results = {}
        for baseline in ["morgan_logistic", "fresh_d6_network"]:
            control = predictions.filter(
                (pl.col("split_id") == split)
                & (pl.col("model") == baseline)
                & pl.col("label").is_finite()
            ).sort("compound_id")
            assert control["compound_id"].to_list() == network["compound_id"].to_list()
            np.testing.assert_array_equal(labels, control["label"].to_numpy())
            baseline_ap = np.array(
                [
                    average_precision_score(labels, control["score"].to_numpy(), sample_weight=w)
                    for w in weights
                ]
            )
            baseline_mcc = weighted_mcc(labels, control["call"].to_numpy(), weights)
            baseline_results[baseline] = (baseline_ap, baseline_mcc)
            for metric, delta in [
                ("AP", network_ap - baseline_ap),
                ("nested_calls_MCC", network_mcc - baseline_mcc),
            ]:
                lo, hi = np.quantile(delta, [0.025, 0.975])
                comparisons.append(
                    dict(
                        split_id=split,
                        comparison=f"multiassay minus {baseline}",
                        metric=metric,
                        groups=len(unique),
                        bootstrap_mean=float(delta.mean()),
                        lower95=float(lo),
                        upper95=float(hi),
                    )
                )
        for index, metric in enumerate(["AP", "nested_calls_MCC"]):
            delta = (
                baseline_results["fresh_d6_network"][index]
                - baseline_results["morgan_logistic"][index]
            )
            lo, hi = np.quantile(delta, [0.025, 0.975])
            comparisons.append(
                dict(
                    split_id=split,
                    comparison="D6 network minus morgan_logistic",
                    metric=metric,
                    groups=len(unique),
                    bootstrap_mean=float(delta.mean()),
                    lower95=float(lo),
                    upper95=float(hi),
                )
            )
    pl.DataFrame(comparisons).write_csv(out / "paired-comparisons.csv")
    text = [
        "# Independent CYP2D6 restart: consolidated development evidence",
        "",
        "All models below were trained from source molecules/assays with no historical "
        "ensemble, fitted encoder, prediction or call budget. Four outer/three inner "
        "folds exclude each held-out molecule's complete assay record from fitting. "
        "Outer calls use inner-selected thresholds, evaluated against all supplied labels.",
        "",
        "| Family | Validation | AP | Nested MCC | AUROC |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for _version, _run, summary in summaries:
        for row in summary["aggregate"]:
            text.append(
                f"| {row['model']} | {row['split_id']} | {row['ap']:.4f} | "
                f"{row['bootstrap_mcc']:.4f} | {row['auc']:.4f} |"
            )
    text += [
        "",
        "## Paired development uncertainty",
        "",
        "Both sides use the same labelled molecules and the same scaffold-group resamples. "
        "Intervals are conditional on saved predictions; they do not include retraining "
        "uncertainty, model-selection effects, or transfer to the unknown Live scored population.",
        "",
        "| Split | Comparison | Metric | Mean delta | 95% interval |",
        "| --- | --- | --- | ---: | --- |",
    ]
    for row in comparisons:
        text.append(
            f"| {row['split_id']} | {row['comparison']} | {row['metric']} | "
            f"{row['bootstrap_mean']:+.4f} | [{row['lower95']:+.4f}, {row['upper95']:+.4f}] |"
        )
    text += [
        "",
        "## Unsubmitted diagnostic operating points",
        "",
        "| Family | Validation | Threshold | Blind positive calls |",
        "| --- | --- | ---: | ---: |",
    ]
    for _version, _run, summary in summaries:
        for row in summary["operating_points"]:
            text.append(
                f"| {row['model']} | {row['split_id']} | {row['threshold']:.4f} | "
                f"{row['blind_positive']} |"
            )
    text += [
        "",
        "No count was selected in advance. These point thresholds are fitted on "
        "pooled development OOF scores and are not separate unbiased tests. No submission "
        "payload was prepared and no upload was made.",
        "",
        "## Provenance",
        "",
    ]
    for version, run, _summary in summaries:
        text.append(f"- {version}: `classification/{run}`.")
    (out / "report.md").write_text("\n".join(text) + "\n")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                comparisons=comparisons,
                sources=cfg["sources"],
                results=[s for _, _, s in summaries],
            ),
            indent=2,
        )
        + "\n"
    )
    print(pl.DataFrame(comparisons))


if __name__ == "__main__":
    main()
