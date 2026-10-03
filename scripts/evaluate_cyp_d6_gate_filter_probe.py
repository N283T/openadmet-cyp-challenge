"""Offline evaluation of a gate-filtered D6 call vector against the deployed recipe.

The proposal keeps the deployed score, the equal percentile-rank mean of the September 15
blind score and `pool_lgbm`, and removes from eligibility any row whose Morgan activity-gate
p_active falls below a threshold, then takes the same number of calls. The call count is
held fixed so the comparison is a pure identity swap.

Part 1 rebuilds the deployed recipe's out-of-fold analogue `l4_rank_mean` from the frozen
held-fold scores of run 20261002T215345688167Z-7920b49a, forms the filtered variants at four
thresholds plus the two multiplier forms already tested against the legacy control, and scores
them on an extended rate grid under both metrics of run 20261003T024837735720Z-008f8427:
Metric O scores officially labelled rows only with the classes reweighted to a target scored
prevalence, and Metric M is the composition-matched two-stratum mixture at the identified
mixture fraction. Uncertainty is a paired group bootstrap against the unfiltered control at the
historical 49/750 operating point and at each variant's nested-selected rate.

Part 2 applies the filtered ordering to the eight scored D6 submission states and asks whether
it separates the removed from the added calls in the direction the public metric moved, under
the frozen concordance design of run 20261002T221918697661Z-2d8d20f5.

Part 3 is a blind swap inventory at 49 calls: overlap, additions and removals against the
deployed and September 15 vectors, the p_active distribution of each group, and the
similarity-weighted top-10 labelled positive rate reused from run
20261002T223515148996Z-ccb05cf3.

Nothing is fitted and no candidate, prepared submission or upload is created here.
"""

import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl

HISTORICAL_RATE = 49 / 750
FILTER_PRIORITY = 10.0


def percentile_rank(values):
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values))
    ranks[order] = np.arange(len(values))
    return (ranks + 0.5) / len(values)


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


def permutation_difference(removed, added, permutations, seed):
    """Two-sided permutation test of mean(added) - mean(removed)."""
    values = np.concatenate([removed, added])
    n_removed = len(removed)
    observed = float(added.mean() - removed.mean())
    rng = np.random.default_rng(seed)
    draws = values[np.argsort(rng.random((permutations, len(values))), axis=1)]
    diff = draws[:, n_removed:].mean(axis=1) - draws[:, :n_removed].mean(axis=1)
    extreme = int((np.abs(diff) >= abs(observed) - 1e-12).sum())
    return observed, (1 + extreme) / (1 + permutations)


def variant_scores(cfg, base_rank, p_active):
    """Ordering values for the control, the gate filters and the two multiplier forms."""
    scores = {cfg["control"]: base_rank.copy()}
    for tau in cfg["thresholds"]:
        scores[f"l4_filter_tau{tau:g}"] = base_rank + FILTER_PRIORITY * (p_active >= tau)
    scores["l4_x_p_active"] = base_rank * p_active
    scores["l4_x_sqrt_p_active"] = base_rank * np.sqrt(p_active)
    assert set(scores) == set(cfg["variants"]), (sorted(scores), sorted(cfg["variants"]))
    return scores


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


def build_l4_oof(oof, components, control_model):
    """OOF analogue of the deployed payload: equal per-fold percentile-rank mean."""
    ranked = (
        oof.filter(pl.col("model").is_in(components))
        .with_columns(
            (
                (pl.col("score").rank("ordinal").over(["split_id", "fold", "model"]) - 0.5)
                / pl.len().over(["split_id", "fold", "model"])
            ).alias("pct")
        )
        .group_by(["split_id", "fold", "compound_id"])
        .agg(pl.col("pct").mean().alias("l4"), pl.len().alias("members"))
    )
    assert ranked["members"].min() == len(components)
    return ranked.drop("members")


def part1(source, cfg, out):
    """Rate curves, nested selection and paired bootstrap under Metric O and Metric M."""
    oof = pl.read_parquet(source / cfg["inputs"]["testlike_oof"])
    assignments = pl.read_parquet(source / cfg["inputs"]["assignments"])
    assert oof["endpoint"].unique().to_list() == [cfg["endpoint"]]
    rates = sorted({*cfg["rates"], HISTORICAL_RATE})
    bootstrap_rates = [min(rates, key=lambda r: abs(r - v)) for v in cfg["bootstrap_rates"]]
    assert all(
        abs(a - b) < 1e-9 for a, b in zip(bootstrap_rates, cfg["bootstrap_rates"], strict=True)
    )
    l4 = build_l4_oof(oof, cfg["l4_components"], cfg["gate_model"])
    gate = oof.filter(pl.col("model") == cfg["gate_model"]).join(
        l4, on=["split_id", "fold", "compound_id"], how="left"
    )
    assert gate["l4"].null_count() == 0

    curve_rows, boot_rows, optimum_rows, nested_rows = [], [], [], []
    for split in sorted(gate["split_id"].unique().to_list()):
        frame = (
            gate.filter(pl.col("split_id") == split)
            .sort("compound_id")
            .join(
                assignments.filter(pl.col("strategy") == split).select("compound_id", "group"),
                on="compound_id",
                how="left",
                maintain_order="left",
            )
            .with_columns(pl.col(["y_official", "y_testlike"]).fill_null(np.nan))
        )
        assert frame["group"].null_count() == 0 and frame.height == 6146
        fold = frame["fold"].to_numpy()
        base_rank = np.empty(frame.height)
        for f in np.unique(fold):
            rows = np.flatnonzero(fold == f)
            base_rank[rows] = percentile_rank(frame["l4"].to_numpy()[rows])
        p_active = frame["p_active"].to_numpy()
        scores = variant_scores(cfg, base_rank, p_active)
        role = frame["pool_role"].to_numpy()
        y_official = frame["y_official"].to_numpy()
        y_testlike = frame["y_testlike"].to_numpy()
        group_ids = frame["group"].to_numpy().astype(str)
        official_mask = role == "official"
        pseudo_mask = role == "pseudo_negative"

        metrics = {}
        for prior in cfg["metric_o_priors"]:
            metrics[("O", prior)] = (official_weights(y_official, official_mask, prior), y_official)
        for phi in cfg["metric_m_phi"]:
            metrics[("M", phi)] = (mixture_weights(official_mask, pseudo_mask, phi), y_testlike)

        calls_grid = {
            model: {r: top_calls(score, fold, r) for r in rates} for model, score in scores.items()
        }
        for (metric, parameter), (weights, labels) in metrics.items():
            # Nested rate selection: the rate is chosen on the other four folds, same metric.
            nested_calls = {}
            for model in scores:
                selected = np.zeros(frame.height, dtype=bool)
                chosen = []
                for f in np.unique(fold):
                    other = fold != f
                    inner = {
                        r: weighted_mcc(labels[other], calls_grid[model][r][other], weights[other])
                        for r in rates
                    }
                    rate = max(rates, key=lambda r: (inner[r], -r))
                    chosen.append(rate)
                    held = fold == f
                    selected[held] = calls_grid[model][rate][held]
                nested_calls[model] = selected
                nested_rows.append(
                    dict(
                        split_id=split,
                        metric=metric,
                        parameter=parameter,
                        model=model,
                        n_calls=int(selected.sum()),
                        mcc=weighted_mcc(labels, selected, weights),
                        selected_rates=",".join(f"{r:g}" for r in chosen),
                    )
                )
            for model in scores:
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
                            blind_equivalent_calls=int(round(rate * 750)),
                            mcc=float(mcc_from_counts(tp, fp, fn, tn)),
                            weighted_precision=tp / (tp + fp) if tp + fp > 0 else None,
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
                        best_calls_of_750=best["blind_equivalent_calls"],
                        mcc_at_49=next(
                            r["mcc"] for r in subset if abs(r["rate"] - HISTORICAL_RATE) < 1e-9
                        ),
                    )
                )
                if model == cfg["control"]:
                    continue
                for rate in bootstrap_rates:
                    boot_rows.append(
                        dict(
                            split_id=split,
                            metric=metric,
                            parameter=parameter,
                            model=model,
                            control=cfg["control"],
                            operating_point=f"{rate:g}",
                            **group_bootstrap(
                                labels,
                                group_ids,
                                weights,
                                calls_grid[model][rate],
                                calls_grid[cfg["control"]][rate],
                                cfg["bootstrap_samples"],
                                cfg["seed"],
                            ),
                        )
                    )
                boot_rows.append(
                    dict(
                        split_id=split,
                        metric=metric,
                        parameter=parameter,
                        model=model,
                        control=cfg["control"],
                        operating_point="nested",
                        **group_bootstrap(
                            labels,
                            group_ids,
                            weights,
                            nested_calls[model],
                            nested_calls[cfg["control"]],
                            cfg["bootstrap_samples"],
                            cfg["seed"],
                        ),
                    )
                )
    pl.DataFrame(curve_rows).write_csv(out / "metric_curves.csv")
    pl.DataFrame(optimum_rows).write_csv(out / "rate_optima.csv")
    pl.DataFrame(nested_rows).write_csv(out / "nested_selection.csv")
    pl.DataFrame(boot_rows).write_csv(out / "bootstrap_deltas.csv")
    return curve_rows, optimum_rows, nested_rows, boot_rows


def blind_table(source, cfg):
    """Blind ordering values of every variant over the 750 blind compounds."""
    test = pl.read_parquet(source / cfg["inputs"]["test_scores"])
    legacy = (
        pl.read_parquet(source / cfg["inputs"]["sep15_blind_scores"])
        .filter(pl.col("model") == cfg["legacy_blind_model"])
        .select("compound_id", pl.col("score").alias("legacy"))
    )
    pool = test.filter(pl.col("model") == cfg["pool_blind_model"]).select(
        "compound_id", pl.col("score").alias("pool")
    )
    gate = test.filter(pl.col("model") == cfg["gate_model"]).select("compound_id", "p_active")
    names = pl.read_parquet(source / cfg["inputs"]["test_mapping"]).select(
        "compound_id", pl.col("molecule_name").alias("Molecule_Name")
    )
    table = (
        legacy.join(pool, on="compound_id", validate="1:1")
        .join(gate, on="compound_id", validate="1:1")
        .join(names, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    assert table.height == 750
    legacy_rank = percentile_rank(table["legacy"].to_numpy())
    pool_rank = percentile_rank(table["pool"].to_numpy())
    base_rank = 0.5 * (legacy_rank + pool_rank)
    p_active = table["p_active"].to_numpy()
    scores = variant_scores(cfg, base_rank, p_active)
    table = table.with_columns(
        pl.Series("legacy_rank", legacy_rank),
        pl.Series("pool_rank", pool_rank),
        pl.Series("l4", base_rank),
        *[pl.Series(f"score_{m}", percentile_rank(v)) for m, v in scores.items()],
    )
    return table, scores


def load_payloads(source, cfg, names):
    sets, rows = {}, []
    for date, alias in cfg["payload_inputs"].items():
        payload = pl.read_parquet(source / alias).join(
            names, on="Molecule_Name", how="left", validate="1:1"
        )
        assert payload.height == 750 and payload["compound_id"].null_count() == 0, alias
        sets[date] = set(payload.filter(pl.col(cfg["endpoint"]))["compound_id"].to_list())
        rows.append(dict(payload=date, alias=alias, positives=len(sets[date])))
    return sets, rows


def part2(source, cfg, table, out):
    """Historical transition concordance of the filtered orderings."""
    names = table.select("compound_id", "Molecule_Name")
    sets, payload_rows = load_payloads(source, cfg, names)
    columns = [f"score_{m}" for m in cfg["variants"]]
    values = {
        c: dict(zip(table["compound_id"].to_list(), table[c].to_list(), strict=True))
        for c in columns
    }
    transition_rows = []
    for t_index, (a, b) in enumerate(cfg["transitions"]):
        removed = sorted(sets[a] - sets[b])
        added = sorted(sets[b] - sets[a])
        lb_delta = cfg["live_mcc"][b] - cfg["live_mcc"][a]
        for c_index, column in enumerate(columns):
            removed_values = np.array([values[column][c] for c in removed])
            added_values = np.array([values[column][c] for c in added])
            observed, p_value = permutation_difference(
                removed_values, added_values, cfg["permutations"], [cfg["seed"], t_index, c_index]
            )
            transition_rows.append(
                dict(
                    transition=f"{a}->{b}",
                    lb_delta=lb_delta,
                    column=column,
                    n_removed=len(removed),
                    n_added=len(added),
                    mean_removed=float(removed_values.mean()),
                    mean_added=float(added_values.mean()),
                    difference=observed,
                    permutation_p=p_value,
                    sign_agreement=bool(np.sign(observed) == np.sign(lb_delta)),
                )
            )
    transitions = pl.DataFrame(transition_rows)
    transitions.write_csv(out / "transitions.csv")
    pl.DataFrame(payload_rows).write_csv(out / "payload_summary.csv")
    agreement = (
        transitions.group_by("column")
        .agg(
            pl.col("sign_agreement").sum().alias("agreeing_transitions"),
            pl.len().alias("transitions"),
            (pl.col("permutation_p") < 0.05).sum().alias("significant_transitions"),
            pl.col("permutation_p").min().alias("min_permutation_p"),
        )
        .sort("agreeing_transitions", descending=True)
    )
    agreement = agreement.with_columns(
        pl.lit(True).alias("contains_deployed_components"),
        pl.struct("column")
        .map_elements(
            lambda s: bool(
                transitions.filter(
                    (pl.col("column") == s["column"]) & (pl.col("transition") == "sep15->oct03")
                )["sign_agreement"].item()
            ),
            return_dtype=pl.Boolean,
        )
        .alias("agrees_on_sep15_to_oct03"),
    )
    agreement.write_csv(out / "concordance_agreement.csv")
    return transition_rows, agreement.to_dicts(), sets


def part3(source, cfg, table, sets, out):
    """Blind swap inventory of each variant at the deployed call count."""
    neighbours = pl.read_parquet(source / cfg["inputs"]["neighbour_rates"])
    table = table.join(neighbours, on="compound_id", how="left", validate="1:1")
    assert table["top10_weighted_positive_rate"].null_count() == 0
    deployed = np.array([c in sets[cfg["deployed_payload"]] for c in table["compound_id"]])
    sep15 = np.array([c in sets[cfg["sep15_payload"]] for c in table["compound_id"]])
    p_active = table["p_active"].to_numpy()
    neighbour = table["top10_weighted_positive_rate"].to_numpy()
    rows, call_rows = [], []
    for count in cfg["blind_call_counts"]:
        for model in cfg["variants"]:
            score = table[f"score_{model}"].to_numpy()
            order = np.argsort(-score, kind="stable")
            calls = np.zeros(len(score), dtype=bool)
            calls[order[:count]] = True
            for reference, name in ((deployed, "deployed_oct03"), (sep15, "sep15")):
                groups = dict(
                    retained=calls & reference,
                    added=calls & ~reference,
                    removed=~calls & reference,
                )
                row = dict(
                    model=model,
                    n_calls=count,
                    reference=name,
                    overlap=int(groups["retained"].sum()),
                    added=int(groups["added"].sum()),
                    removed=int(groups["removed"].sum()),
                    changed=int((calls != reference).sum()),
                )
                for group, mask in groups.items():
                    row[f"mean_p_active_{group}"] = (
                        float(p_active[mask].mean()) if mask.any() else None
                    )
                    row[f"min_p_active_{group}"] = (
                        float(p_active[mask].min()) if mask.any() else None
                    )
                    row[f"mean_neighbour_rate_{group}"] = (
                        float(neighbour[mask].mean()) if mask.any() else None
                    )
                rows.append(row)
            call_rows.append(
                table.select("compound_id", "Molecule_Name", "p_active", "l4")
                .with_columns(
                    pl.lit(model).alias("model"),
                    pl.lit(count).alias("n_calls"),
                    pl.Series("call", calls),
                    pl.Series("deployed_call", deployed),
                    pl.Series("sep15_call", sep15),
                    table["top10_weighted_positive_rate"],
                )
                .filter(pl.col("call") | pl.col("deployed_call"))
            )
    pl.DataFrame(rows).write_csv(out / "blind_swaps.csv")
    pl.concat(call_rows).write_parquet(out / "blind_calls.parquet")
    return rows


def main():
    started = time.time()
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    curve_rows, optimum_rows, nested_rows, boot_rows = part1(source, cfg, out)
    table, _ = blind_table(source, cfg)
    transition_rows, agreement, sets = part2(source, cfg, table, out)
    swap_rows = part3(source, cfg, table, sets, out)

    # Verdict against the rule declared in the config before any number was read.
    phi = cfg["identified_phi"]
    prior = cfg["primary_metric_o_prior"]
    verdict_rows = []
    for model in cfg["variants"]:
        if model == cfg["control"]:
            continue
        m_rows = [
            r
            for r in boot_rows
            if r["model"] == model
            and r["metric"] == "M"
            and abs(r["parameter"] - phi) < 1e-9
            and r["operating_point"] == f"{HISTORICAL_RATE:g}"
        ]
        o_rows = [
            r
            for r in boot_rows
            if r["model"] == model
            and r["metric"] == "O"
            and abs(r["parameter"] - prior) < 1e-9
            and r["operating_point"] == f"{HISTORICAL_RATE:g}"
        ]
        assert len(m_rows) == 2 and len(o_rows) == 2
        verdict_rows.append(
            dict(
                model=model,
                metric_m_positive_both_splits=all(r["ci_low"] > 0 for r in m_rows),
                metric_m_min_ci_low=min(r["ci_low"] for r in m_rows),
                metric_m_min_delta=min(r["delta"] for r in m_rows),
                metric_o_not_negative_either_split=all(r["ci_high"] >= 0 for r in o_rows),
                metric_o_min_delta=min(r["delta"] for r in o_rows),
                metric_o_max_ci_high=max(r["ci_high"] for r in o_rows),
                metric_o_min_ci_high=min(r["ci_high"] for r in o_rows),
                candidate_worthy=bool(
                    all(r["ci_low"] > 0 for r in m_rows) and all(r["ci_high"] >= 0 for r in o_rows)
                ),
            )
        )
    pl.DataFrame(verdict_rows).write_csv(out / "verdict.csv")

    summary = dict(
        config=cfg,
        rate_optima=optimum_rows,
        nested_selection=nested_rows,
        bootstrap=boot_rows,
        concordance_agreement=agreement,
        blind_swaps=swap_rows,
        verdict=verdict_rows,
        runtime_seconds=round(time.time() - started, 1),
        caveats=[
            "Metric M scores pseudo-negative rows whose labels are a function of primary log2fc, "
            "and p_active is a calibrated predictor of primary log2fc, so a gate filter raises "
            "Metric M precision partly by construction; see the gate diagnosis run",
            "Metric O conditions on the strong-primary region and cannot reward correct "
            "down-ranking of weak-primary rows",
            "The filtered blind ordering contains both components of the score that selected the "
            "deployed payload, so its agreement on the sep15->oct03 transition is not independent",
            "Blind rows were never primary screened; the blind p_active values are extrapolation",
            "Nothing is refitted; all scores are frozen",
        ],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")

    primary_curve = [
        r
        for r in curve_rows
        if (r["metric"] == "O" and abs(r["parameter"] - prior) < 1e-9)
        or (r["metric"] == "M" and abs(r["parameter"] - phi) < 1e-9)
    ]
    report = [
        "# Gate-filtered D6 vector against the deployed recipe",
        "",
        "## Rate optima and the value at 49 calls",
        "",
        markdown_table(
            optimum_rows,
            [
                "split_id",
                "metric",
                "parameter",
                "model",
                "best_rate",
                "best_calls_of_750",
                "best_mcc",
                "mcc_at_49",
            ],
        ),
        "",
        "## Paired group bootstrap against the unfiltered control",
        "",
        markdown_table(
            boot_rows,
            [
                "split_id",
                "metric",
                "parameter",
                "model",
                "operating_point",
                "delta",
                "ci_low",
                "ci_high",
                "p_delta_gt_0",
            ],
        ),
        "",
        "## Nested rate selection",
        "",
        markdown_table(
            nested_rows,
            ["split_id", "metric", "parameter", "model", "n_calls", "mcc", "selected_rates"],
        ),
        "",
        "## Rate curve at the primary metric parameters",
        "",
        markdown_table(
            primary_curve,
            [
                "split_id",
                "metric",
                "parameter",
                "model",
                "rate",
                "blind_equivalent_calls",
                "mcc",
                "weighted_precision",
            ],
        ),
        "",
        "## Historical transition concordance",
        "",
        markdown_table(
            transition_rows,
            [
                "transition",
                "lb_delta",
                "column",
                "n_removed",
                "n_added",
                "mean_removed",
                "mean_added",
                "difference",
                "permutation_p",
                "sign_agreement",
            ],
        ),
        "",
        markdown_table(
            agreement,
            [
                "column",
                "agreeing_transitions",
                "transitions",
                "significant_transitions",
                "min_permutation_p",
                "agrees_on_sep15_to_oct03",
                "contains_deployed_components",
            ],
        ),
        "",
        "## Blind swap inventory at the deployed call count",
        "",
        markdown_table(
            swap_rows,
            [
                "model",
                "n_calls",
                "reference",
                "overlap",
                "added",
                "removed",
                "mean_p_active_retained",
                "mean_p_active_added",
                "mean_p_active_removed",
                "mean_neighbour_rate_retained",
                "mean_neighbour_rate_added",
                "mean_neighbour_rate_removed",
            ],
        ),
        "",
        "## Verdict against the rule declared in the config",
        "",
        markdown_table(
            verdict_rows,
            [
                "model",
                "metric_m_positive_both_splits",
                "metric_m_min_delta",
                "metric_m_min_ci_low",
                "metric_o_not_negative_either_split",
                "metric_o_min_delta",
                "metric_o_min_ci_high",
                "candidate_worthy",
            ],
        ),
        "",
    ]
    (out / "report.md").write_text("\n".join(report) + "\n")
    print(
        json.dumps(
            dict(runtime_seconds=summary["runtime_seconds"], verdict=verdict_rows), indent=2
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
