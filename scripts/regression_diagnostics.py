"""Diagnose immutable regression OOF and test predictions without fitting corrections."""

import importlib.util
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from scipy.stats import pearsonr, spearmanr

from cyp_challenge.storage import write_json
from cyp_challenge.tracks import ENDPOINTS


def distribution(values):
    x = np.asarray(values, dtype=float)
    assert len(x) and np.isfinite(x).all()
    return dict(
        n=len(x),
        mean=float(x.mean()),
        sd=float(x.std(ddof=0)),
        **{
            name: float(v)
            for name, v in zip(
                ["min", "q05", "q25", "median", "q75", "q95", "max"],
                np.quantile(x, [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1]),
                strict=True,
            )
        },
        fraction_below_4=float((x < 4).mean()),
        fraction_above_6=float((x >= 6).mean()),
    )


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    eps = ENDPOINTS["regression"]

    def read(name):
        return pl.read_parquet(next(inp.glob(f"*-{name}.parquet")))

    def regression(frame):
        return frame.filter(pl.col("endpoint").is_in(eps))

    spec = importlib.util.spec_from_file_location(
        "official", next(inp.glob("*-custom_scoring_functions.py"))
    )
    official = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(official)
    labels, oof = regression(read("measurements")), regression(read("oof"))
    test, members = regression(read("ensemble_probabilities")), regression(read("fold_predictions"))
    assignments = read("assignments").filter(pl.col("strategy") == cfg["strategy"])
    for frame in [labels, oof, test]:
        assert frame.select("compound_id", "endpoint").n_unique() == frame.height
    assert members.select("compound_id", "endpoint", "fold").n_unique() == members.height
    assert oof.height == labels.height
    assert set(oof["fold"]) == set(range(5))
    joined = labels.join(oof, on=["compound_id", "endpoint"], validate="1:1")
    assert joined.height == labels.height and (joined["value"] == joined["y_true"]).all()
    joined = joined.join(
        assignments.select("compound_id", pl.col("fold").alias("expected_fold")),
        on="compound_id",
        how="left",
        validate="m:1",
    )
    assert joined["expected_fold"].null_count() == 0
    assert (joined["fold"] == joined["expected_fold"]).all()
    assert set(joined["compound_id"]).isdisjoint(test["compound_id"])
    assert test.height == 750 * len(eps) and members.height == 5 * test.height
    averaged = (
        members.group_by("compound_id", "endpoint")
        .agg(pl.col("prediction").mean().alias("recomputed"), pl.len().alias("members"))
        .join(test, on=["compound_id", "endpoint"], validate="1:1")
    )
    assert (averaged["members"] == 5).all()
    assert np.allclose(averaged["recomputed"], averaged["prediction"], atol=1e-7, rtol=1e-7)
    for col in ["y_true", "prediction", "conf_low", "conf_high"]:
        assert joined[col].null_count() == 0 and np.isfinite(joined[col].to_numpy()).all()
    assert (joined["conf_low"] <= joined["conf_high"]).all()
    joined = joined.with_columns(
        (pl.col("prediction") - pl.col("y_true")).alias("residual"),
        (pl.col("conf_high") - pl.col("conf_low")).alias("ci_width"),
        (
            (pl.col("prediction") - pl.col("conf_high")).clip(lower_bound=0)
            + (pl.col("conf_low") - pl.col("prediction")).clip(lower_bound=0)
        ).alias("soft_error"),
        pl.lit(cfg["split_run"] + ":umap").alias("split_id"),
    )
    stats, scores, bins, disagreement, correlations = [], [], [], [], []
    fig, axes = plt.subplots(4, 3, figsize=(16, 17), layout="constrained")

    def score(frame):
        y, p = frame["y_true"].to_numpy(), frame["prediction"].to_numpy()
        lo, hi = frame["conf_low"].to_numpy(), frame["conf_high"].to_numpy()
        denominator = float((np.maximum(y.mean() - hi, 0) + np.maximum(lo - y.mean(), 0)).mean())
        return dict(
            n=len(y),
            ST_RAE=float(
                official.rae_soft_threshold_absolute_error(y, p, y_true_upper=hi, y_true_lower=lo)
            ),
            MAE=float(np.abs(p - y).mean()),
            bias=float((p - y).mean()),
            RMSE=float(np.sqrt(np.mean((p - y) ** 2))),
            Spearman=float(spearmanr(y, p).statistic),
            Pearson=float(pearsonr(y, p).statistic),
            R2=float(1 - np.sum((p - y) ** 2) / np.sum((y - y.mean()) ** 2)),
            pred_sd_over_label_sd=float(p.std() / y.std()),
            mean_soft_error=float(frame["soft_error"].mean()),
            st_denominator_per_row=denominator,
            ci_coverage=float(((p >= lo) & (p <= hi)).mean()),
            below_ci=float((p < lo).mean()),
            above_ci=float((p > hi).mean()),
            median_ci_width=float(np.median(hi - lo)),
        )

    for row, ep in enumerate(eps):
        f = joined.filter(pl.col("endpoint") == ep)
        t = test.filter(pl.col("endpoint") == ep).sort("compound_id")
        m = members.filter(pl.col("endpoint") == ep)
        for fold in [-1, 0, 1, 2, 3, 4]:
            held = f if fold == -1 else f.filter(pl.col("fold") == fold)
            scores.append(dict(endpoint=ep, fold=fold, **score(held)))
            for name, values in [("labels", held["y_true"]), ("oof", held["prediction"])]:
                stats.append(dict(endpoint=ep, fold=fold, population=name, **distribution(values)))
            if fold >= 0:
                train = f.filter(pl.col("fold") != fold)
                stats.append(
                    dict(
                        endpoint=ep,
                        fold=fold,
                        population="outer_train_labels",
                        **distribution(train["y_true"]),
                    )
                )
                stats.append(
                    dict(
                        endpoint=ep,
                        fold=fold,
                        population="test_member",
                        **distribution(m.filter(pl.col("fold") == fold)["prediction"]),
                    )
                )
        stats.append(
            dict(endpoint=ep, fold=-1, population="test_ensemble", **distribution(t["prediction"]))
        )
        for kind in ["prediction_quantile", "label_range"]:
            x = f["prediction" if kind == "prediction_quantile" else "y_true"].to_numpy()
            edges = (
                np.unique(np.quantile(x, np.linspace(0, 1, cfg["prediction_bins"] + 1)[1:-1]))
                if kind == "prediction_quantile"
                else np.array(cfg["label_bin_edges"])
            )
            grouped = f.with_columns(pl.Series("bin", np.searchsorted(edges, x, side="right")))
            for b in sorted(grouped["bin"].unique()):
                part = grouped.filter(pl.col("bin") == b)
                bins.append(
                    dict(
                        endpoint=ep,
                        kind=kind,
                        bin=int(b),
                        n=part.height,
                        lower=float(edges[b - 1]) if b else None,
                        upper=float(edges[b]) if b < len(edges) else None,
                        label_mean=float(part["y_true"].mean()),
                        prediction_mean=float(part["prediction"].mean()),
                        bias=float(part["residual"].mean()),
                        MAE=float(part["residual"].abs().mean()),
                        mean_soft_error=float(part["soft_error"].mean()),
                    )
                )
        wide = m.pivot(on="fold", index="compound_id", values="prediction").sort("compound_id")
        matrix = wide.select([str(i) for i in range(5)]).to_numpy()
        assert np.isfinite(matrix).all()
        sd = matrix.std(axis=1)
        disagreement.append(dict(endpoint=ep, **distribution(sd)))
        for i in range(5):
            for j in range(i + 1, 5):
                correlations.append(
                    dict(
                        endpoint=ep,
                        fold_i=i,
                        fold_j=j,
                        Pearson=float(pearsonr(matrix[:, i], matrix[:, j]).statistic),
                        Spearman=float(spearmanr(matrix[:, i], matrix[:, j]).statistic),
                    )
                )
        wide.with_columns(
            pl.Series("member_sd", sd),
            pl.lit(ep).alias("endpoint"),
            pl.lit(cfg["split_run"] + ":umap").alias("split_id"),
        ).write_parquet(out / f"{ep}-test-disagreement.parquet")
        y, p = f["y_true"].to_numpy(), f["prediction"].to_numpy()
        name = ep.split("_")[0]
        ax = axes[row, 0]
        for values, label, color in [
            (y, "Training labels", "#424242"),
            (p, "OOF: one held-out model", "#2878b5"),
            (t["prediction"].to_numpy(), "Test: five-model mean", "#d95f02"),
        ]:
            xx = np.sort(values)
            ax.step(xx, np.arange(1, len(xx) + 1) / len(xx), where="post", label=label, color=color)
        ax.set(title=f"{name}: distributions", xlabel="pIC50", ylabel="Cumulative fraction")
        ax.legend(fontsize=8)
        ax = axes[row, 1]
        ax.hexbin(y, p, gridsize=35, mincnt=1, cmap="Blues")
        limits = [min(y.min(), p.min()), max(y.max(), p.max())]
        ax.plot(limits, limits, "--", color="grey")
        ax.set(title=f"{name}: OOF predictions", xlabel="Observed pIC50", ylabel="Predicted pIC50")
        ax = axes[row, 2]
        b = [r for r in bins if r["endpoint"] == ep and r["kind"] == "prediction_quantile"]
        ax.plot([r["prediction_mean"] for r in b], [r["bias"] for r in b], "o-")
        ax.axhline(0, color="grey", linestyle="--")
        ax.set(
            title=f"{name}: OOF bias by predicted quintile",
            xlabel="Mean predicted pIC50",
            ylabel="Mean (prediction − label)",
        )
    fig.savefig(out / "regression-diagnostics.png", dpi=150)
    plt.close(fig)
    for name, rows in [
        ("distributions", stats),
        ("metrics", scores),
        ("residual_bins", bins),
        ("test_disagreement", disagreement),
        ("test_member_correlations", correlations),
    ]:
        pl.DataFrame(rows, infer_schema_length=None).write_csv(out / f"{name}.csv")
    joined.write_parquet(out / "regression-oof-diagnostics.parquet")
    saved = read("metrics").filter(pl.col("endpoint").is_in(eps))
    for r in scores:
        old = saved.filter(
            (pl.col("endpoint") == r["endpoint"]) & (pl.col("fold") == r["fold"])
        ).row(0, named=True)
        for oldkey, key in [("ST-RAE", "ST_RAE"), ("MAE", "MAE"), ("Spearman", "Spearman")]:
            assert np.isclose(old[oldkey], r[key], rtol=1e-10, atol=1e-10)
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            audit_passed=True,
            pooled=[r for r in scores if r["fold"] == -1],
            limitations=[
                "No blinded labels: test bias or correction cannot be identified.",
                "OOF uses one model per row; test uses five. Member distributions are saved.",
                "Fold variation includes different chemistry; it is not seed uncertainty.",
                "True-label residual groups show shrinkage, not an actionable calibration rule.",
                "No calibration, training, or submission performed.",
            ],
        ),
    )
    print(
        pl.DataFrame([r for r in scores if r["fold"] == -1]).select(
            "endpoint",
            "n",
            "ST_RAE",
            "MAE",
            "bias",
            "Spearman",
            "pred_sd_over_label_sd",
            "ci_coverage",
        )
    )


if __name__ == "__main__":
    main()
