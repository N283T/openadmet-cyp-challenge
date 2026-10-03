"""D6 recipe comparison under metrics that do not use pseudo-negative labels as the yardstick.

The test-like pool comparison scores pseudo-negative rows, so a model trained on them is
rewarded on its own hypothesis. This script rescores the same frozen held-fold predictions
of run 20261002T215345688167Z-7920b49a under two metrics that avoid that.

Metric O scores only officially labelled rows and reweights the two classes so that the
scored prevalence equals a target prior near the 7.06% implied by the public aggregates.
It never credits a model for a pseudo-negative row, and it is therefore a lower bound on
the pool model's value: it also never credits correct down-ranking of weak-primary rows.

Metric M treats the blind scored set as a mixture of a strong-primary-like stratum, sampled
by the officially labelled rows with their real labels, and a weak-primary-like stratum,
sampled by the pseudo-negative rows. The mixture fraction phi is identified rather than
assumed: the public prevalence identity gives 0.0706 = phi * p_strong + (1 - phi) * p_weak
with p_strong = 324/1497 observed and p_weak bounded by the three-enzyme measurement, which
pins phi near one third. Metric M at phi = 1 reduces to the unweighted official metric.

Calls come from every development molecule of the held fold, as the protocol requires. The
rate grid runs far past the 10% boundary that truncated every earlier selection. Nothing is
fitted, no candidate or prepared submission is created, and no API is contacted.
"""

import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl

ENDPOINT = "CYP2D6_is_TDI"
CONTROL = "official_lgbm"
L4 = "l4_rank_mean"
HISTORICAL_RATE = 49 / 750


def mcc_from_counts(tp, fp, fn, tn):
    tp, fp, fn, tn = (np.asarray(v, dtype=float) for v in (tp, fp, fn, tn))
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return np.divide(tp * tn - fp * fn, denom, out=np.zeros_like(denom), where=denom > 0)


def weighted_counts(y, calls, weights):
    labels = y.astype(bool)
    return (
        float((weights * (labels & calls)).sum()),
        float((weights * (~labels & calls)).sum()),
        float((weights * (labels & ~calls)).sum()),
        float((weights * (~labels & ~calls)).sum()),
    )


def weighted_mcc(y, calls, weights):
    return float(mcc_from_counts(*weighted_counts(y, calls, weights)))


def effective_n(weights):
    positive = weights[weights > 0]
    return float(positive.sum() ** 2 / (positive**2).sum()) if positive.size else 0.0


def official_weights(y, mask, prior):
    """Metric O: class weights on officially labelled rows giving scored prevalence `prior`."""
    weights = np.zeros(len(y))
    labelled = mask & np.isfinite(y)
    positive = labelled & (y == 1)
    negative = labelled & (y == 0)
    n = int(labelled.sum())
    weights[positive] = prior * n / max(int(positive.sum()), 1)
    weights[negative] = (1 - prior) * n / max(int(negative.sum()), 1)
    return weights


def mixture_weights(official, pseudo, phi):
    """Metric M: stratum weights giving the strong-primary-like stratum mass `phi`."""
    weights = np.zeros(len(official))
    n = int(official.sum()) + int(pseudo.sum())
    weights[official] = phi * n / max(int(official.sum()), 1)
    weights[pseudo] = (1 - phi) * n / max(int(pseudo.sum()), 1)
    return weights


def top_calls(score, fold, rate):
    calls = np.zeros(len(score), dtype=bool)
    for f in np.unique(fold):
        rows = np.flatnonzero(fold == f)
        k = int(round(rate * len(rows)))
        calls[rows[np.argsort(-score[rows], kind="stable")[:k]]] = True
    return calls


def group_bootstrap(y, groups, weights, calls_a, calls_b, samples, seed):
    """Paired group bootstrap of pooled weighted MCC(a) - MCC(b)."""
    keep = np.isfinite(y) & (weights > 0)
    labels = y[keep].astype(bool)
    weight = weights[keep]
    _, index = np.unique(groups[keep], return_inverse=True)
    count = index.max() + 1
    counts = []
    for calls in (calls_a[keep], calls_b[keep]):
        counts.append(
            [
                np.bincount(index, weights=weight * (labels & calls), minlength=count),
                np.bincount(index, weights=weight * (~labels & calls), minlength=count),
                np.bincount(index, weights=weight * (labels & ~calls), minlength=count),
                np.bincount(index, weights=weight * (~labels & ~calls), minlength=count),
            ]
        )
    draws = np.random.default_rng(seed).integers(0, count, (samples, count))
    resample = np.zeros((samples, count))
    for b in range(samples):
        resample[b] = np.bincount(draws[b], minlength=count)
    delta = mcc_from_counts(*(resample @ c for c in counts[0])) - mcc_from_counts(
        *(resample @ c for c in counts[1])
    )
    point = mcc_from_counts(*(c.sum() for c in counts[0])) - mcc_from_counts(
        *(c.sum() for c in counts[1])
    )
    return dict(
        delta=float(point),
        boot_mean=float(delta.mean()),
        ci_low=float(np.quantile(delta, 0.025)),
        ci_high=float(np.quantile(delta, 0.975)),
        p_delta_gt_0=float((delta > 0).mean()),
        groups=int(count),
        molecules=int(keep.sum()),
    )


def fmt(value, digits=3):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "–"
    return f"{value:.{digits}f}" if isinstance(value, float) else str(value)


def markdown_table(rows, columns):
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(c)) for c in columns) + " |")
    return "\n".join(lines)


def main():
    started = time.time()
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    oof = pl.read_parquet(source / cfg["inputs"]["testlike_oof"])
    assignments = pl.read_parquet(source / cfg["inputs"]["assignments"])
    assert oof["endpoint"].unique().to_list() == [ENDPOINT]
    rates = sorted({*cfg["rates"], HISTORICAL_RATE})
    # Snap configured comparison rates onto the grid so float spellings cannot diverge.
    bootstrap_rates = [min(rates, key=lambda r: abs(r - value)) for value in cfg["bootstrap_rates"]]
    assert all(
        abs(a - b) < 1e-9 for a, b in zip(bootstrap_rates, cfg["bootstrap_rates"], strict=True)
    )
    # OOF analogue of the deployed L4 payload: equal percentile-rank mean of the legacy
    # September 15 refit and pool_lgbm, formed inside each held fold.
    components = cfg["l4_components"]
    ranked = (
        oof.filter(pl.col("model").is_in(components))
        .with_columns(
            (
                (pl.col("score").rank("ordinal").over(["split_id", "fold", "model"]) - 0.5)
                / pl.len().over(["split_id", "fold", "model"])
            ).alias("pct")
        )
        .group_by(["split_id", "fold", "compound_id"])
        .agg(pl.col("pct").mean().alias("score"), pl.len().alias("members"))
    )
    assert ranked["members"].min() == len(components)
    template = oof.filter(pl.col("model") == CONTROL).drop("score")
    l4 = template.join(
        ranked.drop("members"), on=["split_id", "fold", "compound_id"], how="left"
    ).with_columns(pl.lit(L4).alias("model"))
    assert l4["score"].null_count() == 0
    oof = pl.concat([oof, l4.select(oof.columns)])
    models = oof["model"].unique(maintain_order=True).to_list()
    assert CONTROL in models and L4 in models

    # Identify the mixture fraction from the published prevalence instead of assuming it.
    official_rows = oof.filter(
        (pl.col("split_id") == oof["split_id"][0]) & (pl.col("model") == CONTROL)
    )
    p_strong = float(official_rows.filter(pl.col("pool_role") == "official")["y_official"].mean())
    identification = []
    for p_weak in cfg["weak_prevalence_scenarios"]:
        for prevalence in cfg["scored_prevalence_scenarios"]:
            phi = (prevalence - p_weak) / (p_strong - p_weak)
            identification.append(
                dict(
                    scored_prevalence=prevalence,
                    p_strong=p_strong,
                    p_weak=p_weak,
                    phi=phi,
                    interpretable=0.0 <= phi <= 1.0,
                )
            )
    print(json.dumps(identification), flush=True)
    pl.DataFrame(identification).write_csv(out / "phi_identification.csv")
    phi_grid = sorted(
        {
            *cfg["phi_grid"],
            *[round(r["phi"], 4) for r in identification if r["interpretable"]],
        }
    )

    curve_rows, boot_rows, optimum_rows = [], [], []
    for split in oof["split_id"].unique().to_list():
        groups = (
            oof.filter((pl.col("split_id") == split) & (pl.col("model") == CONTROL))
            .sort("compound_id")
            .join(
                assignments.filter(pl.col("strategy") == split).select("compound_id", "group"),
                on="compound_id",
                how="left",
            )
        )
        assert groups["group"].null_count() == 0
        frames, calls_grid = {}, {}
        for model in models:
            frame = (
                oof.filter((pl.col("split_id") == split) & (pl.col("model") == model))
                .sort("compound_id")
                .with_columns(pl.col(["y_official", "y_testlike"]).fill_null(np.nan))
            )
            assert frame["compound_id"].equals(groups["compound_id"])
            frames[model] = frame
            score = frame["score"].to_numpy()
            fold = frame["fold"].to_numpy()
            calls_grid[model] = {r: top_calls(score, fold, r) for r in rates}
        role = frames[CONTROL]["pool_role"].to_numpy()
        y_official = frames[CONTROL]["y_official"].to_numpy()
        y_testlike = frames[CONTROL]["y_testlike"].to_numpy()
        group_ids = groups["group"].to_numpy().astype(str)
        official_mask = role == "official"
        pseudo_mask = role == "pseudo_negative"

        # Metric O: official rows only, class-reweighted to each target prior.
        metric_weights = {}
        for prior in cfg["scored_prevalence_scenarios"]:
            metric_weights[("O", prior)] = (
                official_weights(y_official, official_mask, prior),
                y_official,
            )
        # Metric M: stratified mixture of real official labels and pseudo-negative rows.
        for phi in phi_grid:
            metric_weights[("M", phi)] = (
                mixture_weights(official_mask, pseudo_mask, phi),
                y_testlike,
            )
        for (metric, parameter), (weights, labels) in metric_weights.items():
            for model in models:
                for rate in rates:
                    calls = calls_grid[model][rate]
                    tp, fp, fn, tn = weighted_counts(labels, calls, weights)
                    curve_rows.append(
                        dict(
                            split_id=split,
                            metric=metric,
                            parameter=parameter,
                            model=model,
                            rate=rate,
                            n_calls=int(calls.sum()),
                            scored_calls=int((calls & (weights > 0)).sum()),
                            mcc=float(mcc_from_counts(tp, fp, fn, tn)),
                            weighted_prevalence=(tp + fn) / (tp + fp + fn + tn),
                            effective_n=effective_n(weights),
                        )
                    )
                subset = [
                    r
                    for r in curve_rows
                    if r["split_id"] == split
                    and r["metric"] == metric
                    and r["parameter"] == parameter
                    and r["model"] == model
                ]
                best = max(subset, key=lambda r: (r["mcc"], -r["rate"]))
                optimum_rows.append(
                    dict(
                        split_id=split,
                        metric=metric,
                        parameter=parameter,
                        model=model,
                        best_rate=best["rate"],
                        best_mcc=best["mcc"],
                        best_calls_of_750=int(round(best["rate"] * 750)),
                        at_grid_edge=best["rate"] in (min(rates), max(rates)),
                    )
                )
            for model in models:
                if model == CONTROL:
                    continue
                for rate in bootstrap_rates:
                    boot_rows.append(
                        dict(
                            split_id=split,
                            metric=metric,
                            parameter=parameter,
                            model=model,
                            control=CONTROL,
                            rate=rate,
                            **group_bootstrap(
                                labels,
                                group_ids,
                                weights,
                                calls_grid[model][rate],
                                calls_grid[CONTROL][rate],
                                cfg["bootstrap_samples"],
                                cfg["seed"],
                            ),
                        )
                    )
    curves = pl.DataFrame(curve_rows)
    curves.write_csv(out / "metric_curves.csv")
    pl.DataFrame(optimum_rows).write_csv(out / "rate_optima.csv")
    pl.DataFrame(boot_rows).write_csv(out / "bootstrap_deltas.csv")

    # Delta versus the mixture fraction, and every sign change over the non-degenerate grid.
    positive_phi = [phi for phi in phi_grid if phi > 0]
    delta_rows, crossover = [], []
    for split in curves["split_id"].unique().to_list():
        for model in models:
            if model == CONTROL:
                continue
            for rate in bootstrap_rates:
                series = []
                for phi in positive_phi:
                    rows = curves.filter(
                        (pl.col("split_id") == split)
                        & (pl.col("metric") == "M")
                        & (pl.col("parameter") == phi)
                        & (pl.col("rate") == rate)
                    )
                    delta = (
                        rows.filter(pl.col("model") == model)["mcc"].item()
                        - rows.filter(pl.col("model") == CONTROL)["mcc"].item()
                    )
                    series.append((phi, delta))
                    delta_rows.append(
                        dict(split_id=split, model=model, rate=rate, phi=phi, delta=delta)
                    )
                crossings = [
                    (a[0] + b[0]) / 2
                    for a, b in zip(series, series[1:], strict=False)
                    if (a[1] > 0) != (b[1] > 0)
                ]
                identified = min(series, key=lambda item: abs(item[0] - cfg["identified_phi"]))
                crossover.append(
                    dict(
                        split_id=split,
                        model=model,
                        rate=rate,
                        delta_at_identified_phi=identified[1],
                        phi_used=identified[0],
                        delta_at_phi_1=series[-1][1],
                        crossings=",".join(f"{c:.3f}" for c in crossings) or "none",
                        n_crossings=len(crossings),
                    )
                )
    pl.DataFrame(delta_rows).write_csv(out / "delta_vs_phi.csv")
    pl.DataFrame(crossover).write_csv(out / "crossover.csv")

    summary = dict(
        config=cfg,
        endpoint=ENDPOINT,
        p_strong=p_strong,
        phi_identification=identification,
        phi_grid=phi_grid,
        rate_optima=optimum_rows,
        crossover=crossover,
        bootstrap=boot_rows,
        runtime_seconds=round(time.time() - started, 1),
        caveats=[
            "Metric O conditions on the strong-primary region, so it cannot reward correct "
            "down-ranking of weak-primary rows and understates the pool model",
            "Metric M assumes the blind scored set is a two-stratum mixture and that each "
            "stratum is sampled by its development counterpart",
            "phi is identified from the published prevalence, which is itself a reconstruction",
            "Scores are frozen from one run; nothing is refitted here",
            "The L4 analogue is rebuilt from two frozen OOF columns, not refitted",
        ],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    report = [
        "# D6 recipes under metrics that do not score pseudo-negative labels",
        "",
        "## Mixture fraction identified from the published prevalence",
        "",
        markdown_table(
            identification, ["scored_prevalence", "p_strong", "p_weak", "phi", "interpretable"]
        ),
        "",
        "## Rate optima",
        "",
        markdown_table(
            optimum_rows,
            [
                "split_id",
                "metric",
                "parameter",
                "model",
                "best_rate",
                "best_mcc",
                "best_calls_of_750",
                "at_grid_edge",
            ],
        ),
        "",
        "## Paired group bootstrap versus the official-only control",
        "",
        markdown_table(
            boot_rows,
            [
                "split_id",
                "metric",
                "parameter",
                "model",
                "rate",
                "delta",
                "ci_low",
                "ci_high",
                "p_delta_gt_0",
            ],
        ),
        "",
        "## Crossover mixture fraction",
        "",
        markdown_table(
            crossover,
            [
                "split_id",
                "model",
                "rate",
                "delta_at_identified_phi",
                "delta_at_phi_1",
                "crossings",
            ],
        ),
        "",
    ]
    (out / "report.md").write_text("\n".join(report))
    print(json.dumps(dict(runtime_seconds=summary["runtime_seconds"])), flush=True)


if __name__ == "__main__":
    main()
