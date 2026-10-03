"""Diagnose the D6 activity gate before trusting or discarding it as a call filter.

The deployed D6 vector calls 49 of 750 blind rows positive and 40 of those sit below
gate probability 0.5. The cross-enzyme measurement of run 20261003T023544961394Z-95af3eea
found 0 positives among 638 labelled rows at primary log2fc >= -0.25. Those two facts are
only in tension if a low gate probability placed a row in that relabelled region, so this
script measures the three links in that chain separately.

Part A compares the blind p_active distribution with the development distribution, overall
and by whether the development row is a screened primary hit, a screened non-hit or
unscreened. Part B measures the gate's reliability directly on development held folds:
binned predicted p_active against the observed rate of log2fc < the gate threshold.
Part C isolates the unscreened development rows, which the blind rows resemble in carrying
no primary measurement. Part D measures what the gate says about the endpoint rather than
about the primary screen: the officially labelled positive rate by p_active bin, which is
the quantity the filter proposal depends on. Part E decomposes p_active into primary-screen
regions and composes the region-wise positive rates into an implied blind positive rate.
Part F quantifies how strongly p_active separates the official stratum from the
pseudo-negative stratum, which bounds how circular a pseudo-negative-labelled metric is
when it is used to judge a gate filter.

Nothing is fitted; every score is frozen. No candidate, prepared submission or upload.
"""

import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

HELD_FOLD_NOTE = "p_active values are held-fold predictions; no row scored its own fit"


def wilson(successes, total, z=1.959963984540054):
    if total == 0:
        return (float("nan"), float("nan"))
    p = successes / total
    denom = 1 + z**2 / total
    centre = (p + z**2 / (2 * total)) / denom
    half = z * np.sqrt(p * (1 - p) / total + z**2 / (4 * total**2)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def describe(name, values, thresholds):
    values = np.asarray(values, dtype=float)
    row = dict(
        population=name,
        n=int(values.size),
        mean=float(values.mean()),
        sd=float(values.std(ddof=1)) if values.size > 1 else float("nan"),
    )
    for q in (0.1, 0.25, 0.5, 0.75, 0.9):
        row[f"q{int(q * 100)}"] = float(np.quantile(values, q))
    for tau in thresholds:
        row[f"frac_ge_{tau:g}"] = float((values >= tau).mean())
    return row


def auroc(score, label):
    """Mann-Whitney AUROC of `score` against the boolean `label`."""
    score, label = np.asarray(score, dtype=float), np.asarray(label, dtype=bool)
    if label.all() or not label.any():
        return float("nan")
    ranks = stats.rankdata(score)
    positives = int(label.sum())
    negatives = int((~label).sum())
    return float((ranks[label].sum() - positives * (positives + 1) / 2) / (positives * negatives))


def bin_index(values, edges):
    return np.clip(np.searchsorted(edges, values, side="right") - 1, 0, len(edges) - 2)


def load_primary(path, aliases, enzyme):
    mapping = aliases.select("molecule_name", "compound_id").unique()
    assert not mapping["molecule_name"].is_duplicated().any()
    raw = pl.read_csv(path).rename({"Molecule_Name": "molecule_name"})
    linked = raw.join(mapping, on="molecule_name", how="left", validate="m:1")
    assert linked["compound_id"].null_count() == 0
    concentrations = linked["concentration_M"].unique().to_list()
    assert len(concentrations) == 1 and 0.000049 <= concentrations[0] <= 0.000051
    return (
        linked.filter(pl.col("enzyme") == enzyme)
        .group_by("compound_id")
        .agg(pl.col("log2fc_estimate").median().alias("log2fc"))
    )


def fmt(value, digits=4):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "–"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def markdown_table(rows, columns, digits=4):
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(c), digits) for c in columns) + " |")
    return "\n".join(lines)


def main():
    started = time.time()
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    inputs = cfg["inputs"]
    gate_model = cfg["gate_model"]
    thresholds = cfg["thresholds"]
    edges = np.asarray(cfg["reliability_bins"], dtype=float)
    coarse = np.asarray(cfg["coarse_bins"], dtype=float)

    oof = pl.read_parquet(source / inputs["testlike_oof"]).filter(pl.col("model") == gate_model)
    assert oof["endpoint"].unique().to_list() == [cfg["endpoint"]]
    blind = pl.read_parquet(source / inputs["blind_scores"]).filter(pl.col("model") == gate_model)
    assert blind.height == 750 and blind["p_active"].is_finite().all()
    aliases = pl.read_parquet(source / inputs["aliases"])
    primary = load_primary(source / inputs["primary_screen"], aliases, cfg["enzyme"])
    oof = oof.join(primary, on="compound_id", how="left")

    splits = sorted(oof["split_id"].unique().to_list())
    blind_pa = blind["p_active"].to_numpy()

    # ---- Part A: distribution of p_active, blind against development strata ----
    distribution_rows = [describe("blind_750", blind_pa, thresholds)]
    shift_rows = []
    for split in splits:
        frame = oof.filter(pl.col("split_id") == split)
        pa = frame["p_active"].to_numpy()
        log2fc = frame["log2fc"].to_numpy()
        role = frame["pool_role"].to_numpy()
        screened = np.isfinite(log2fc)
        groups = {
            f"{split}: development_all": np.ones(len(pa), dtype=bool),
            f"{split}: screened_hit (log2fc < -0.5)": screened & (log2fc < -0.5),
            f"{split}: screened_non_hit (log2fc >= -0.5)": screened & (log2fc >= -0.5),
            f"{split}: unscreened (no log2fc)": ~screened,
            f"{split}: role official": role == "official",
            f"{split}: role pseudo_negative": role == "pseudo_negative",
            f"{split}: role borderline": role == "borderline",
            f"{split}: role unlabelled_hit": role == "unlabelled_hit",
            f"{split}: role unscreened": role == "unscreened",
        }
        for name, mask in groups.items():
            distribution_rows.append(describe(name, pa[mask], thresholds))
            ks = stats.ks_2samp(blind_pa, pa[mask])
            shift_rows.append(
                dict(
                    reference=name,
                    n_reference=int(mask.sum()),
                    mean_blind=float(blind_pa.mean()),
                    mean_reference=float(pa[mask].mean()),
                    mean_difference=float(blind_pa.mean() - pa[mask].mean()),
                    ks_statistic=float(ks.statistic),
                    ks_p=float(ks.pvalue),
                    auroc_blind_above_reference=auroc(
                        np.concatenate([blind_pa, pa[mask]]),
                        np.concatenate(
                            [np.ones(len(blind_pa), bool), np.zeros(int(mask.sum()), bool)]
                        ),
                    ),
                )
            )
    pl.DataFrame(distribution_rows).write_csv(out / "p_active_distributions.csv")
    pl.DataFrame(shift_rows).write_csv(out / "p_active_shift.csv")

    # ---- Part B: gate reliability on development held folds ----
    reliability_rows, reliability_summary = [], []
    for split in splits:
        frame = oof.filter((pl.col("split_id") == split) & pl.col("log2fc").is_not_null())
        pa = frame["p_active"].to_numpy()
        active = (frame["log2fc"].to_numpy() < cfg["gate_threshold"]).astype(float)
        index = bin_index(pa, edges)
        gaps, weights = [], []
        for b in range(len(edges) - 1):
            mask = index == b
            if not mask.any():
                continue
            predicted = float(pa[mask].mean())
            observed = float(active[mask].mean())
            low, high = wilson(float(active[mask].sum()), int(mask.sum()))
            reliability_rows.append(
                dict(
                    split_id=split,
                    bin_low=float(edges[b]),
                    bin_high=float(edges[b + 1]),
                    n=int(mask.sum()),
                    predicted=predicted,
                    observed=observed,
                    reliability=observed - predicted,
                    observed_ci_low=low,
                    observed_ci_high=high,
                    predicted_inside_ci=bool(low <= predicted <= high),
                )
            )
            gaps.append(abs(observed - predicted))
            weights.append(int(mask.sum()))
        gaps, weights = np.asarray(gaps), np.asarray(weights, dtype=float)
        reliability_summary.append(
            dict(
                split_id=split,
                n_screened=int(len(pa)),
                base_rate=float(active.mean()),
                mean_p_active=float(pa.mean()),
                expected_calibration_error=float((gaps * weights).sum() / weights.sum()),
                max_bin_deviation=float(gaps.max()),
                brier=float(np.mean((pa - active) ** 2)),
                brier_of_base_rate=float(np.mean((active.mean() - active) ** 2)),
                auroc=auroc(pa, active.astype(bool)),
                bins_with_predicted_inside_ci=int(
                    sum(
                        1
                        for r in reliability_rows
                        if r["split_id"] == split and r["predicted_inside_ci"]
                    )
                ),
                bins=int(len(gaps)),
            )
        )
    pl.DataFrame(reliability_rows).write_csv(out / "gate_reliability.csv")
    pl.DataFrame(reliability_summary).write_csv(out / "gate_reliability_summary.csv")

    # ---- Part C: unscreened development rows against screened ones ----
    unscreened_rows = []
    for split in splits:
        frame = oof.filter(pl.col("split_id") == split)
        pa = frame["p_active"].to_numpy()
        screened = np.isfinite(frame["log2fc"].to_numpy())
        test = stats.mannwhitneyu(pa[~screened], pa[screened], alternative="two-sided")
        unscreened_rows.append(
            dict(
                split_id=split,
                n_unscreened=int((~screened).sum()),
                n_screened=int(screened.sum()),
                mean_unscreened=float(pa[~screened].mean()),
                mean_screened=float(pa[screened].mean()),
                difference=float(pa[~screened].mean() - pa[screened].mean()),
                auroc_unscreened_above_screened=auroc(pa, ~screened),
                mannwhitney_p=float(test.pvalue),
                ks_blind_vs_unscreened=float(stats.ks_2samp(blind_pa, pa[~screened]).statistic),
                ks_blind_vs_screened=float(stats.ks_2samp(blind_pa, pa[screened]).statistic),
            )
        )
    pl.DataFrame(unscreened_rows).write_csv(out / "unscreened_comparison.csv")

    # ---- Part D: what the gate says about the endpoint, not about the screen ----
    endpoint_rows, endpoint_summary = [], []
    for split in splits:
        frame = oof.filter((pl.col("split_id") == split) & (pl.col("pool_role") == "official"))
        pa = frame["p_active"].to_numpy()
        y = frame["y_official"].to_numpy()
        assert np.isfinite(y).all()
        for name, binning in (("fine", edges), ("coarse", coarse)):
            index = bin_index(pa, binning)
            for b in range(len(binning) - 1):
                mask = index == b
                if not mask.any():
                    continue
                low, high = wilson(float(y[mask].sum()), int(mask.sum()))
                endpoint_rows.append(
                    dict(
                        split_id=split,
                        binning=name,
                        bin_low=float(binning[b]),
                        bin_high=float(binning[b + 1]),
                        n=int(mask.sum()),
                        mean_p_active=float(pa[mask].mean()),
                        positives=int(y[mask].sum()),
                        positive_rate=float(y[mask].mean()),
                        ci_low=low,
                        ci_high=high,
                    )
                )
        spearman = stats.spearmanr(pa, y)
        endpoint_summary.append(
            dict(
                split_id=split,
                n_official=int(len(y)),
                positives=int(y.sum()),
                base_rate=float(y.mean()),
                auroc_p_active_for_tdi=auroc(pa, y.astype(bool)),
                spearman=float(spearman.statistic),
                spearman_p=float(spearman.pvalue),
                rate_below_0_2=float(y[pa < 0.2].mean()),
                n_below_0_2=int((pa < 0.2).sum()),
                rate_above_0_5=float(y[pa >= 0.5].mean()),
                n_above_0_5=int((pa >= 0.5).sum()),
            )
        )
    pl.DataFrame(endpoint_rows).write_csv(out / "tdi_rate_by_p_active.csv")
    pl.DataFrame(endpoint_summary).write_csv(out / "tdi_rate_summary.csv")

    # ---- Part E: region decomposition and the implied blind positive rate ----
    regions = cfg["log2fc_regions"]
    decomposition_rows, implied_rows = [], []
    for split in splits:
        frame = oof.filter((pl.col("split_id") == split) & pl.col("log2fc").is_not_null())
        pa = frame["p_active"].to_numpy()
        log2fc = frame["log2fc"].to_numpy()
        index = bin_index(pa, coarse)
        share = {}
        for b in range(len(coarse) - 1):
            mask = index == b
            row = dict(
                split_id=split,
                bin_low=float(coarse[b]),
                bin_high=float(coarse[b + 1]),
                n=int(mask.sum()),
                mean_p_active=float(pa[mask].mean()),
            )
            fractions = []
            for region in regions:
                inside = mask & (log2fc >= region["low"]) & (log2fc < region["high"])
                fraction = float(inside.sum() / max(int(mask.sum()), 1))
                row[f"share_{region['name']}"] = fraction
                fractions.append(fraction)
            share[b] = fractions
            for label, relabelled in cfg["relabelled_rate_scenarios"].items():
                rates = [
                    relabelled if region["name"] == "relabelled" else region["rate"]
                    for region in regions
                ]
                row[f"implied_positive_rate_{label}"] = float(np.dot(fractions, rates))
            decomposition_rows.append(row)
        # Apply the development region shares to the blind p_active distribution.
        blind_index = bin_index(blind_pa, coarse)
        for label, relabelled in cfg["relabelled_rate_scenarios"].items():
            rates = [
                relabelled if region["name"] == "relabelled" else region["rate"]
                for region in regions
            ]
            per_row = np.array([np.dot(share[b], rates) for b in blind_index])
            implied_rows.append(
                dict(
                    split_id=split,
                    scenario=label,
                    relabelled_rate=relabelled,
                    implied_blind_prevalence=float(per_row.mean()),
                    published_scored_prevalence=cfg["published_scored_prevalence"],
                    ratio_to_published=float(per_row.mean() / cfg["published_scored_prevalence"]),
                    implied_rate_blind_below_0_2=float(per_row[blind_pa < 0.2].mean()),
                    implied_rate_blind_0_2_to_0_5=float(
                        per_row[(blind_pa >= 0.2) & (blind_pa < 0.5)].mean()
                    ),
                    implied_rate_blind_above_0_5=float(per_row[blind_pa >= 0.5].mean()),
                )
            )
    pl.DataFrame(decomposition_rows).write_csv(out / "region_decomposition.csv")
    pl.DataFrame(implied_rows).write_csv(out / "implied_blind_rate.csv")

    # ---- Part F: how circular a pseudo-negative metric is for judging a gate filter ----
    circularity_rows = []
    for split in splits:
        frame = oof.filter(pl.col("split_id") == split)
        pa = frame["p_active"].to_numpy()
        role = frame["pool_role"].to_numpy()
        keep = np.isin(role, ["official", "pseudo_negative"])
        circularity_rows.append(
            dict(
                split_id=split,
                n_official=int((role == "official").sum()),
                n_pseudo_negative=int((role == "pseudo_negative").sum()),
                auroc_p_active_official_vs_pseudo=auroc(pa[keep], (role[keep] == "official")),
                mean_official=float(pa[role == "official"].mean()),
                mean_pseudo_negative=float(pa[role == "pseudo_negative"].mean()),
                **{
                    f"official_share_above_{tau:g}": float(
                        ((role == "official") & (pa >= tau)).sum()
                        / max(int(((pa >= tau) & keep).sum()), 1)
                    )
                    for tau in thresholds
                },
            )
        )
    pl.DataFrame(circularity_rows).write_csv(out / "metric_circularity.csv")

    summary = dict(
        config=cfg,
        blind=dict(
            n=int(blind.height),
            mean_p_active=float(blind_pa.mean()),
            **{f"n_ge_{tau:g}": int((blind_pa >= tau).sum()) for tau in thresholds},
        ),
        distributions=distribution_rows,
        shift=shift_rows,
        reliability_summary=reliability_summary,
        unscreened=unscreened_rows,
        endpoint_summary=endpoint_summary,
        implied=implied_rows,
        circularity=circularity_rows,
        runtime_seconds=round(time.time() - started, 1),
        notes=[
            HELD_FOLD_NOTE,
            "p_active targets P(primary log2fc < -1); the cross-enzyme pseudo-negative rule "
            "concerns log2fc >= -0.25, so a low p_active does not place a row in that region",
            "Region-wise positive rates are taken from the frozen measurements of run "
            "20261003T023544961394Z-95af3eea and from D6 labelled rows; the relabelled region "
            "has no D6 labels at all and is varied as a scenario",
            "Blind rows were never primary screened, so every blind statement is extrapolation",
        ],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")

    report = [
        "# Diagnosing the D6 activity gate on blind chemistry",
        "",
        "## p_active distribution, blind against development strata",
        "",
        markdown_table(
            distribution_rows,
            ["population", "n", "mean", "q10", "q50", "q90"]
            + [f"frac_ge_{t:g}" for t in thresholds],
        ),
        "",
        "## Distribution shift of the blind set against each reference",
        "",
        markdown_table(
            shift_rows,
            [
                "reference",
                "n_reference",
                "mean_reference",
                "mean_difference",
                "ks_statistic",
                "ks_p",
                "auroc_blind_above_reference",
            ],
        ),
        "",
        "## Gate reliability on development held folds",
        "",
        markdown_table(
            reliability_rows,
            [
                "split_id",
                "bin_low",
                "bin_high",
                "n",
                "predicted",
                "observed",
                "reliability",
                "observed_ci_low",
                "observed_ci_high",
                "predicted_inside_ci",
            ],
        ),
        "",
        markdown_table(
            reliability_summary,
            [
                "split_id",
                "n_screened",
                "base_rate",
                "mean_p_active",
                "expected_calibration_error",
                "max_bin_deviation",
                "brier",
                "brier_of_base_rate",
                "auroc",
                "bins_with_predicted_inside_ci",
                "bins",
            ],
        ),
        "",
        "## Unscreened development rows",
        "",
        markdown_table(
            unscreened_rows,
            [
                "split_id",
                "n_unscreened",
                "mean_unscreened",
                "mean_screened",
                "difference",
                "auroc_unscreened_above_screened",
                "mannwhitney_p",
                "ks_blind_vs_unscreened",
                "ks_blind_vs_screened",
            ],
        ),
        "",
        "## Officially labelled positive rate by p_active bin",
        "",
        markdown_table(
            [r for r in endpoint_rows if r["binning"] == "coarse"],
            [
                "split_id",
                "bin_low",
                "bin_high",
                "n",
                "mean_p_active",
                "positives",
                "positive_rate",
                "ci_low",
                "ci_high",
            ],
        ),
        "",
        markdown_table(
            endpoint_summary,
            [
                "split_id",
                "n_official",
                "positives",
                "base_rate",
                "auroc_p_active_for_tdi",
                "spearman",
                "spearman_p",
                "n_below_0_2",
                "rate_below_0_2",
                "n_above_0_5",
                "rate_above_0_5",
            ],
        ),
        "",
        "## Primary-screen region decomposition of p_active",
        "",
        markdown_table(
            decomposition_rows,
            ["split_id", "bin_low", "bin_high", "n"]
            + [f"share_{r['name']}" for r in regions]
            + [f"implied_positive_rate_{k}" for k in cfg["relabelled_rate_scenarios"]],
        ),
        "",
        "## Implied blind positive rate from the decomposition",
        "",
        markdown_table(
            implied_rows,
            [
                "split_id",
                "scenario",
                "relabelled_rate",
                "implied_blind_prevalence",
                "ratio_to_published",
                "implied_rate_blind_below_0_2",
                "implied_rate_blind_0_2_to_0_5",
                "implied_rate_blind_above_0_5",
            ],
        ),
        "",
        "## Circularity of a pseudo-negative-labelled metric for a gate filter",
        "",
        markdown_table(
            circularity_rows,
            [
                "split_id",
                "n_official",
                "n_pseudo_negative",
                "auroc_p_active_official_vs_pseudo",
                "mean_official",
                "mean_pseudo_negative",
            ]
            + [f"official_share_above_{t:g}" for t in thresholds],
        ),
        "",
    ]
    (out / "report.md").write_text("\n".join(report) + "\n")
    print(json.dumps(dict(runtime_seconds=summary["runtime_seconds"]), indent=2), flush=True)


if __name__ == "__main__":
    main()
