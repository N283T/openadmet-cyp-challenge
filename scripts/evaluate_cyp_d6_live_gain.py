"""Statistical post-mortem of the 2026-10-03 D6 Live MCC gain.

The public classification leaderboard scores an unknown subset of the 750 blind
molecules and reports bootstrap means and standard deviations of Accuracy,
Precision, Recall, F1 and MCC. This script reconstructs the scored-subset size,
the number of scored positives and each payload's confusion matrix from those
aggregates, sizes the September 15 -> October 3 change in units of individual
row outcomes, asks how often an unskilled call swap of the same size reproduces
a change that large, attributes the new call vector between its two score
components, and extends the historical concordance table by the new transition.

Part 1  Operating-point reconstruction (variance inversion, the prevalence
        identity, and a joint integer fit shared across all eight D6 states).
Part 2  Null resampling of an unskilled k-call swap at the estimated base rate.
Part 3  Attribution of the new vector between the deployed score and pool_lgbm.
Part 4  The eighth transition added to the frozen concordance table.

Nothing here fits a model, touches the challenge API, registers a candidate or
prepares a submission.
"""

import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl

ENDPOINT = "CYP2D6_is_TDI"
BLIND_ROWS = 750
DEPLOYED_BLIND_MODEL = "deployed_mean_block"
METRICS = ("accuracy", "precision", "recall", "f1", "mcc")


def fmt(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def markdown_table(rows, columns):
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(c)) for c in columns) + " |")
    return "\n".join(lines)


def metric_array(tp, fp, fn, tn):
    """Vectorised Accuracy, Precision, Recall, F1 and MCC from confusion counts."""
    tp, fp, fn, tn = (np.asarray(v, dtype=float) for v in (tp, fp, fn, tn))
    total = tp + fp + fn + tn
    accuracy = (tp + tn) / total
    called, positive, predicted = tp + fp, tp + fn, 2 * tp + fp + fn
    precision = np.divide(tp, called, out=np.zeros_like(tp), where=called > 0)
    recall = np.divide(tp, positive, out=np.zeros_like(tp), where=positive > 0)
    f1 = np.divide(2 * tp, predicted, out=np.zeros_like(tp), where=predicted > 0)
    denominator = np.sqrt(called * positive * (tn + fp) * (tn + fn))
    mcc = np.divide(tp * tn - fp * fn, denominator, out=np.zeros_like(tp), where=denominator > 0)
    return dict(accuracy=accuracy, precision=precision, recall=recall, f1=f1, mcc=mcc)


def bootstrap_metrics(tp, fp, fn, tn, samples, rng):
    """Mean and sd of each metric over `samples` multinomial resamples of the scored rows."""
    total = tp + fp + fn + tn
    draws = rng.multinomial(total, np.array([tp, fp, fn, tn], dtype=float) / total, size=samples)
    values = metric_array(draws[:, 0], draws[:, 1], draws[:, 2], draws[:, 3])
    return {k: (float(v.mean()), float(v.std(ddof=1))) for k, v in values.items()}


def read_live_row(path, username):
    """The participant's row of one saved leaderboard CSV."""
    table = pl.read_csv(path)
    row = table.filter(pl.col("Username") == username)
    assert row.height == 1, f"{path}: expected one {username} row, found {row.height}"
    record = row.row(0, named=True)
    return dict(
        submitted=record["Submitted"],
        rank=int(record["Rank"]),
        entrants=table.height,
        accuracy=float(record["Accuracy"]),
        accuracy_sd=float(record["Accuracy std"]),
        precision=float(record["Precision"]),
        precision_sd=float(record["Precision std"]),
        recall=float(record["Recall"]),
        recall_sd=float(record["Recall std"]),
        f1=float(record["F1 Score"]),
        f1_sd=float(record["F1 Score std"]),
        mcc=float(record["MCC"]),
        mcc_sd=float(record["MCC std"]),
    )


def variance_inversion(state):
    """P ~ r(1-r)/sd_r^2 and the analogous inversions on precision and accuracy."""
    return dict(
        positives_from_recall=state["recall"] * (1 - state["recall"]) / state["recall_sd"] ** 2,
        calls_from_precision=state["precision"]
        * (1 - state["precision"])
        / state["precision_sd"] ** 2,
        scored_from_accuracy=state["accuracy"]
        * (1 - state["accuracy"])
        / state["accuracy_sd"] ** 2,
    )


def prevalence_identity(state):
    """prevalence = (1-acc)/(1-2r+r/p); call fraction = prevalence*r/p."""
    ratio = state["recall"] / state["precision"]
    prevalence = (1 - state["accuracy"]) / (1 - 2 * state["recall"] + ratio)
    return dict(prevalence=prevalence, call_fraction=prevalence * ratio)


def best_cell(state, scored, positives, max_calls):
    """Best integer (TP, scored calls) for a given scored size and positive count."""
    true_positive = np.arange(0, positives + 1)[:, None]
    calls = np.arange(0, max_calls + 1)[None, :]
    false_positive = calls - true_positive
    false_negative = positives - true_positive
    true_negative = scored - true_positive - false_positive - false_negative
    valid = (false_positive >= 0) & (true_negative >= 0)
    values = metric_array(
        true_positive + 0 * calls,
        np.maximum(false_positive, 0),
        false_negative + 0 * calls,
        np.maximum(true_negative, 0),
    )
    loss = np.zeros_like(values["mcc"])
    for name in METRICS:
        loss += ((values[name] - state[name]) / state[f"{name}_sd"]) ** 2
    loss = np.where(valid, loss, np.inf)
    index = np.unravel_index(np.argmin(loss), loss.shape)
    return int(true_positive[index[0], 0]), int(calls[0, index[1]]), float(loss[index])


def joint_fit(states, scored_grid, positive_grid):
    """Shared (scored subset size, scored positives) across every saved state."""
    results = []
    for scored in scored_grid:
        for positives in positive_grid:
            if positives >= scored:
                continue
            total, per_state = 0.0, {}
            for name, state in states.items():
                true_positive, calls, loss = best_cell(
                    state, scored, positives, min(state["calls"], scored)
                )
                per_state[name] = dict(true_positive=true_positive, scored_calls=calls, loss=loss)
                total += loss
            results.append(dict(scored=scored, positives=positives, loss=total, states=per_state))
    results.sort(key=lambda r: r["loss"])
    return results


def null_swap(scored, positives, true_positive, scored_calls, total_calls, swaps, draws, rng):
    """Point MCC deltas of an unskilled swap of `swaps` calls, at the base rate."""
    false_positive = scored_calls - true_positive
    unscored_calls = total_calls - scored_calls
    uncalled = BLIND_ROWS - total_calls
    scored_uncalled = scored - scored_calls
    positive_uncalled = positives - true_positive
    negative_uncalled = scored_uncalled - positive_uncalled
    unscored_uncalled = uncalled - scored_uncalled
    assert min(unscored_calls, positive_uncalled, negative_uncalled, unscored_uncalled) >= 0

    removed = rng.multivariate_hypergeometric(
        np.array([true_positive, false_positive, unscored_calls]), swaps, size=draws
    )
    added = rng.multivariate_hypergeometric(
        np.array([positive_uncalled, negative_uncalled, unscored_uncalled]), swaps, size=draws
    )
    new_tp = true_positive - removed[:, 0] + added[:, 0]
    new_fp = false_positive - removed[:, 1] + added[:, 1]
    new_fn = positives - new_tp
    new_tn = scored - new_tp - new_fp - new_fn
    base = float(
        metric_array(
            true_positive,
            false_positive,
            positives - true_positive,
            scored - scored_calls - (positives - true_positive),
        )["mcc"]
    )
    return metric_array(new_tp, new_fp, new_fn, new_tn)["mcc"] - base, base, new_tp, new_fp


def exchangeable_scenarios(base, new, swaps):
    """Compositions of the changed calls consistent with both fitted confusion matrices.

    The retained calls contribute a fixed (true positive, scored call) pair; everything
    else must be carried by the `swaps` removed and the `swaps` added calls.
    """
    scenarios = []
    for retained_tp in range(0, min(base["true_positive"], new["true_positive"]) + 1):
        for retained_calls in range(
            retained_tp, min(base["scored_calls"], new["scored_calls"]) + 1
        ):
            removed_tp = base["true_positive"] - retained_tp
            added_tp = new["true_positive"] - retained_tp
            removed_calls = base["scored_calls"] - retained_calls
            added_calls = new["scored_calls"] - retained_calls
            removed_fp = removed_calls - removed_tp
            added_fp = added_calls - added_tp
            if min(removed_fp, added_fp) < 0:
                continue
            if max(removed_calls, added_calls) > swaps:
                continue
            scenarios.append(
                dict(
                    retained_true_positive=retained_tp,
                    retained_scored_calls=retained_calls,
                    changed_positives=removed_tp + added_tp,
                    changed_negatives=removed_fp + added_fp,
                    changed_unscored=2 * swaps - (removed_tp + added_tp) - (removed_fp + added_fp),
                    removed_true_positive=removed_tp,
                    added_true_positive=added_tp,
                    removed_false_positive=removed_fp,
                    added_false_positive=added_fp,
                )
            )
    return [s for s in scenarios if s["changed_unscored"] >= 0]


def exchangeable_null(base, scenario, swaps, draws, rng):
    """MCC deltas when the changed calls are split at random into added and removed halves."""
    colors = np.array(
        [scenario["changed_positives"], scenario["changed_negatives"], scenario["changed_unscored"]]
    )
    added = rng.multivariate_hypergeometric(colors, swaps, size=draws)
    new_tp = scenario["retained_true_positive"] + added[:, 0]
    new_calls = scenario["retained_scored_calls"] + added[:, 0] + added[:, 1]
    new_fp = new_calls - new_tp
    new_fn = base["positives"] - new_tp
    new_tn = base["scored"] - new_tp - new_fp - new_fn
    base_mcc = float(
        metric_array(
            base["true_positive"],
            base["false_positive"],
            base["false_negative"],
            base["true_negative"],
        )["mcc"]
    )
    return metric_array(new_tp, new_fp, new_fn, new_tn)["mcc"] - base_mcc


def percentile_rank(values):
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values))
    ranks[order] = np.arange(len(values))
    return (ranks + 0.5) / len(values)


def permutation_difference(removed, added, permutations, rng):
    """Two-sided permutation p for mean(added) - mean(removed)."""
    pooled = np.concatenate([removed, added])
    observed = added.mean() - removed.mean()
    n_removed = len(removed)
    count = 0
    for _ in range(permutations):
        shuffled = rng.permutation(pooled)
        if abs(shuffled[n_removed:].mean() - shuffled[:n_removed].mean()) >= abs(observed) - 1e-12:
            count += 1
    return float(observed), (count + 1) / (permutations + 1)


def part1(source, cfg, out):
    """Reconstruct the scored subset and every payload's confusion matrix."""
    states, rows = {}, []
    for name, spec in cfg["live_states"].items():
        state = read_live_row(source / cfg["inputs"][spec["leaderboard"]], cfg["username"])
        assert state["submitted"] == spec["submitted"], (name, state["submitted"])
        state["calls"] = spec["calls"]
        states[name] = state
    rng = np.random.default_rng(cfg["seed"])

    for name, state in states.items():
        inverted = variance_inversion(state)
        identity = prevalence_identity(state)
        rows.append(
            dict(
                state=name,
                submitted=state["submitted"],
                calls=state["calls"],
                mcc=state["mcc"],
                mcc_sd=state["mcc_sd"],
                positives_from_recall_sd=inverted["positives_from_recall"],
                calls_from_precision_sd=inverted["calls_from_precision"],
                scored_from_accuracy_sd=inverted["scored_from_accuracy"],
                prevalence_identity=identity["prevalence"],
                call_fraction_identity=identity["call_fraction"],
                scored_from_identity=inverted["calls_from_precision"] / identity["call_fraction"],
                positives_from_identity=identity["prevalence"]
                * inverted["calls_from_precision"]
                / identity["call_fraction"],
            )
        )

    ranked = joint_fit(
        states,
        range(cfg["scored_grid"][0], cfg["scored_grid"][1] + 1),
        range(cfg["positive_grid"][0], cfg["positive_grid"][1] + 1),
    )
    best = ranked[0]
    scored, positives = best["scored"], best["positives"]

    fit_rows, fitted = [], {}
    for name, state in states.items():
        cell = best["states"][name]
        true_positive, scored_calls = cell["true_positive"], cell["scored_calls"]
        false_positive = scored_calls - true_positive
        false_negative = positives - true_positive
        true_negative = scored - true_positive - false_positive - false_negative
        simulated = bootstrap_metrics(
            true_positive,
            false_positive,
            false_negative,
            true_negative,
            cfg["bootstrap_samples"],
            rng,
        )
        fitted[name] = dict(
            scored=scored,
            positives=positives,
            true_positive=true_positive,
            scored_calls=scored_calls,
            false_positive=false_positive,
            false_negative=false_negative,
            true_negative=true_negative,
            total_calls=state["calls"],
        )
        row = dict(
            state=name,
            submitted=state["submitted"],
            total_calls=state["calls"],
            scored_calls=scored_calls,
            true_positive=true_positive,
            false_positive=false_positive,
            false_negative=false_negative,
            fit_loss=cell["loss"],
        )
        for metric in METRICS:
            row[f"{metric}_reported"] = state[metric]
            row[f"{metric}_simulated"] = simulated[metric][0]
            row[f"{metric}_sd_reported"] = state[f"{metric}_sd"]
            row[f"{metric}_sd_simulated"] = simulated[metric][1]
        fit_rows.append(row)

    # Sensitivity: hold the scale-free identity fixed and vary the scored subset size.
    band = []
    for name in cfg["sensitivity_states"]:
        identity = prevalence_identity(states[name])
        for size in cfg["sensitivity_scored"]:
            scored_calls = identity["call_fraction"] * size
            scored_positives = identity["prevalence"] * size
            true_positive = states[name]["recall"] * scored_positives
            band.append(
                dict(
                    state=name,
                    scored=size,
                    scored_positives=scored_positives,
                    scored_calls=scored_calls,
                    true_positive=true_positive,
                    false_positive=scored_calls - true_positive,
                )
            )

    alternatives = [
        dict(
            rank=i + 1,
            scored=r["scored"],
            positives=r["positives"],
            loss=r["loss"],
            **{f"tp_{k}": v["true_positive"] for k, v in r["states"].items()},
            **{f"calls_{k}": v["scored_calls"] for k, v in r["states"].items()},
        )
        for i, r in enumerate(ranked[: cfg["report_alternatives"]])
    ]

    pl.DataFrame(rows).write_csv(out / "operating_point_estimators.csv")
    pl.DataFrame(fit_rows).write_csv(out / "operating_point_fit.csv")
    pl.DataFrame(band).write_csv(out / "sensitivity_band.csv")
    pl.DataFrame(alternatives).write_csv(out / "joint_fit_alternatives.csv")
    return states, fitted, rows, fit_rows, band, alternatives, scored, positives


def part2(cfg, states, fitted, out):
    """Null distribution of the MCC move produced by an unskilled k-call swap."""
    rng = np.random.default_rng(cfg["seed"] + 1)
    summary, quantile_rows = [], []
    for spec in cfg["null_transitions"]:
        base_name, new_name, swaps = spec["from"], spec["to"], spec["swaps"]
        base = fitted[base_name]
        observed = states[new_name]["mcc"] - states[base_name]["mcc"]
        deltas, base_mcc, new_tp, new_fp = null_swap(
            base["scored"],
            base["positives"],
            base["true_positive"],
            base["scored_calls"],
            base["total_calls"],
            swaps,
            cfg["null_draws"],
            rng,
        )
        greater = float(np.mean(deltas >= observed - 1e-12))
        extreme = float(np.mean(np.abs(deltas) >= abs(observed) - 1e-12))
        summary.append(
            dict(
                transition=f"{base_name}->{new_name}",
                swapped_calls=swaps,
                observed_delta=observed,
                reported_sd_new=states[new_name]["mcc_sd"],
                observed_in_sd=observed / states[new_name]["mcc_sd"],
                null_mean=float(deltas.mean()),
                null_sd=float(deltas.std(ddof=1)),
                null_p05=float(np.quantile(deltas, 0.05)),
                null_p50=float(np.quantile(deltas, 0.50)),
                null_p95=float(np.quantile(deltas, 0.95)),
                null_fraction_at_least=greater,
                null_fraction_as_extreme=extreme,
                base_mcc_point=base_mcc,
                mean_true_positive=float(new_tp.mean()),
                mean_false_positive=float(new_fp.mean()),
            )
        )
        for q in cfg["null_quantiles"]:
            quantile_rows.append(
                dict(
                    transition=f"{base_name}->{new_name}",
                    quantile=q,
                    delta=float(np.quantile(deltas, q)),
                )
            )

    # Does the bootstrap mean move with the point estimate? Check on the headline swap.
    head = cfg["null_transitions"][0]
    base = fitted[head["from"]]
    deltas, base_mcc, new_tp, new_fp = null_swap(
        base["scored"],
        base["positives"],
        base["true_positive"],
        base["scored_calls"],
        base["total_calls"],
        head["swaps"],
        cfg["bias_check_draws"],
        rng,
    )
    base_boot = bootstrap_metrics(
        base["true_positive"],
        base["false_positive"],
        base["false_negative"],
        base["true_negative"],
        cfg["bias_check_bootstrap"],
        rng,
    )["mcc"][0]
    boot_deltas = []
    for index in range(cfg["bias_check_draws"]):
        tp, fp = int(new_tp[index]), int(new_fp[index])
        fn = base["positives"] - tp
        tn = base["scored"] - tp - fp - fn
        boot_deltas.append(
            bootstrap_metrics(tp, fp, fn, tn, cfg["bias_check_bootstrap"], rng)["mcc"][0]
            - base_boot
        )
    boot_deltas = np.asarray(boot_deltas)
    bias_check = dict(
        draws=cfg["bias_check_draws"],
        point_mean=float(deltas.mean()),
        point_sd=float(deltas.std(ddof=1)),
        bootstrap_mean=float(boot_deltas.mean()),
        bootstrap_sd=float(boot_deltas.std(ddof=1)),
        mean_difference=float(boot_deltas.mean() - deltas.mean()),
    )
    # Sharper null: the recipe cannot order within the set of changed calls.
    exchange_rows, scenario_rows = [], []
    for spec in cfg["null_transitions"]:
        base_name, new_name, swaps = spec["from"], spec["to"], spec["swaps"]
        base, new = fitted[base_name], fitted[new_name]
        observed = states[new_name]["mcc"] - states[base_name]["mcc"]
        scenarios = exchangeable_scenarios(base, new, swaps)
        # The "central" scenario is the neutral composition: the retained, removed and
        # added groups share one scored fraction and one precision among their scored calls.
        retained_size = base["total_calls"] - swaps
        base_precision = base["true_positive"] / base["scored_calls"]
        distance = []
        for sc in scenarios:
            retained_calls = sc["retained_scored_calls"]
            removed_calls = sc["removed_true_positive"] + sc["removed_false_positive"]
            added_calls = sc["added_true_positive"] + sc["added_false_positive"]
            retained_rate = retained_calls / retained_size
            value = (retained_rate - removed_calls / swaps) ** 2 + (
                retained_rate - added_calls / swaps
            ) ** 2
            for calls, positives in (
                (retained_calls, sc["retained_true_positive"]),
                (removed_calls, sc["removed_true_positive"]),
            ):
                value += (positives / calls - base_precision) ** 2 if calls else 1.0
            distance.append(value)
        central = int(np.argmin(distance))
        pooled, p_values = [], []
        for index, scenario in enumerate(scenarios):
            scenario["central"] = index == central
            deltas = exchangeable_null(base, scenario, swaps, cfg["exchangeable_draws"], rng)
            p_value = float(np.mean(deltas >= observed - 1e-12))
            pooled.append(deltas)
            p_values.append(p_value)
            scenario_rows.append(
                dict(
                    transition=f"{base_name}->{new_name}",
                    **scenario,
                    fraction_at_least=p_value,
                    null_mean=float(deltas.mean()),
                )
            )
        pooled = np.concatenate(pooled)
        exchange_rows.append(
            dict(
                transition=f"{base_name}->{new_name}",
                swapped_calls=swaps,
                observed_delta=observed,
                scenarios=len(scenarios),
                pooled_mean=float(pooled.mean()),
                pooled_sd=float(pooled.std(ddof=1)),
                pooled_fraction_at_least=float(np.mean(pooled >= observed - 1e-12)),
                central_fraction_at_least=p_values[central],
                min_fraction_at_least=min(p_values),
                median_fraction_at_least=float(np.median(p_values)),
                max_fraction_at_least=max(p_values),
            )
        )

    pl.DataFrame(summary).write_csv(out / "null_distribution.csv")
    pl.DataFrame(quantile_rows).write_csv(out / "null_quantiles.csv")
    pl.DataFrame(exchange_rows).write_csv(out / "exchangeable_null.csv")
    pl.DataFrame(scenario_rows).write_csv(out / "exchangeable_scenarios.csv")
    return summary, quantile_rows, bias_check, exchange_rows, scenario_rows


def load_blind(source, cfg):
    """The 750 blind molecules with the deployed score, pool_lgbm and the gate columns."""
    inputs = cfg["inputs"]
    mapping = (
        pl.read_parquet(source / inputs["test_mapping"])
        .select("compound_id", "molecule_name")
        .unique()
    )
    assert mapping.height == BLIND_ROWS
    scores = pl.read_parquet(source / inputs["test_scores"])
    wide = mapping
    for key in cfg["score_columns"]:
        model, field = key.split(".")
        wide = wide.join(
            scores.filter(pl.col("model") == model).select("compound_id", pl.col(field).alias(key)),
            on="compound_id",
            how="inner",
            validate="1:1",
        )
    legacy = (
        pl.read_parquet(source / inputs["sep15_blind_scores"])
        .filter(pl.col("model") == DEPLOYED_BLIND_MODEL)
        .select("compound_id", legacy="score")
    )
    wide = wide.join(legacy, on="compound_id", how="inner", validate="1:1").sort("compound_id")
    assert wide.height == BLIND_ROWS
    return wide


def payload_calls(source, cfg, alias, mapping):
    payload = pl.read_parquet(source / cfg["inputs"][alias])
    calls = (
        payload.select(molecule_name="Molecule_Name", call=ENDPOINT)
        .join(mapping, on="molecule_name", how="inner", validate="1:1")
        .filter(pl.col("call"))["compound_id"]
        .to_list()
    )
    return set(calls)


def part3(source, cfg, blind, out):
    """Attribute the new call vector between the deployed score and pool_lgbm."""
    mapping = blind.select("compound_id", "molecule_name")
    sep15 = payload_calls(source, cfg, "payload_sep15", mapping)
    oct03 = payload_calls(source, cfg, "payload_oct03", mapping)
    assert len(sep15) == 49 and len(oct03) == 49

    compound = blind["compound_id"].to_numpy()
    legacy_rank = percentile_rank(blind["legacy"].to_numpy())
    pool_rank = percentile_rank(blind[cfg["pool_column"]].to_numpy())
    mean_rank = 0.5 * (legacy_rank + pool_rank)
    # The recorded diagnostic stores the percentile rank of the rank mean, not the mean itself.
    combined_rank = percentile_rank(mean_rank)

    def top(rank, k):
        return set(compound[np.argsort(-rank, kind="stable")[:k]])

    assert top(legacy_rank, 49) == sep15, "deployed score must reproduce the Sep 15 calls"
    assert top(mean_rank, 49) == oct03, "the rank mean must reproduce the Oct 3 calls"

    pool_top = top(pool_rank, 49)
    added, removed = oct03 - sep15, sep15 - oct03
    retained = oct03 & sep15

    group = np.where(
        np.isin(compound, list(added)),
        "added",
        np.where(
            np.isin(compound, list(removed)),
            "removed",
            np.where(np.isin(compound, list(retained)), "retained", "uncalled"),
        ),
    )
    table = blind.with_columns(
        group=pl.Series(group),
        legacy_rank=pl.Series(legacy_rank),
        pool_rank=pl.Series(pool_rank),
        mean_rank=pl.Series(mean_rank),
        combined_rank=pl.Series(combined_rank),
        pool_top49=pl.Series(np.isin(compound, list(pool_top))),
    ).filter(pl.col("group") != "uncalled")

    diagnostics = (
        pl.read_parquet(source / cfg["inputs"]["call_diagnostics"])
        .filter(pl.col("variant") == cfg["variant"])
        .select(
            "molecule_name",
            "group",
            "legacy_rank",
            "pool_rank",
            "combined_rank",
            "p_active_morgan",
            "top10_weighted_positive_rate",
        )
    )
    assert diagnostics.height == 66
    check = table.join(diagnostics, on="molecule_name", how="inner", suffix="_ref")
    assert check.height == 66, "changed-call identities must match the recorded diagnostic"
    identity = dict(
        group_mismatches=int((check["group"] != check["group_ref"]).sum()),
        max_legacy_rank_difference=float(
            (check["legacy_rank"] - check["legacy_rank_ref"]).abs().max()
        ),
        max_pool_rank_difference=float((check["pool_rank"] - check["pool_rank_ref"]).abs().max()),
        max_combined_rank_difference=float(
            (check["combined_rank"] - check["combined_rank_ref"]).abs().max()
        ),
        max_p_active_difference=float(
            (check["gate_x_tdi.p_active"] - check["p_active_morgan"]).abs().max()
        ),
    )
    assert identity["group_mismatches"] == 0
    assert max(v for k, v in identity.items() if k != "group_mismatches") < 1e-9

    # The CSV projection of the same diagnostic must list exactly the changed calls.
    changed = (
        pl.read_csv(source / cfg["inputs"]["changed_calls"])
        .filter(pl.col("variant") == cfg["variant"])
        .select("group", "molecule_name", "legacy_rank", "pool_rank", "combined_rank")
    )
    recorded = table.filter(pl.col("group") != "retained")
    paired = recorded.join(changed, on="molecule_name", how="inner", suffix="_csv")
    identity["changed_calls_rows"] = changed.height
    identity["changed_calls_matched"] = paired.height
    identity["changed_calls_group_mismatches"] = int((paired["group"] != paired["group_csv"]).sum())
    identity["changed_calls_max_rank_difference"] = max(
        float((paired[column] - paired[f"{column}_csv"]).abs().max())
        for column in ("legacy_rank", "pool_rank", "combined_rank")
    )
    assert identity["changed_calls_rows"] == 34
    assert identity["changed_calls_matched"] == 34
    assert identity["changed_calls_group_mismatches"] == 0
    assert identity["changed_calls_max_rank_difference"] < 1e-9

    table = check.drop(
        "group_ref", "legacy_rank_ref", "pool_rank_ref", "combined_rank_ref", "p_active_morgan"
    )
    group_rows = []
    for name in ("retained", "removed", "added"):
        part = table.filter(pl.col("group") == name)
        group_rows.append(
            dict(
                group=name,
                n=part.height,
                legacy_rank=float(part["legacy_rank"].mean()),
                pool_rank=float(part["pool_rank"].mean()),
                mean_rank=float(part["mean_rank"].mean()),
                p_active=float(part["gate_x_tdi.p_active"].mean()),
                neighbour_positive_rate=float(part["top10_weighted_positive_rate"].mean()),
                in_pool_top49=int(part["pool_top49"].sum()),
            )
        )
    attribution = dict(
        pool_top49_overlap_sep15=len(pool_top & sep15),
        pool_top49_overlap_oct03=len(pool_top & oct03),
        pool_top49_overlap_retained=len(pool_top & retained),
        additions_in_pool_top49=len(added & pool_top),
        removals_in_pool_top49=len(removed & pool_top),
        closer_vector="oct03" if len(pool_top & oct03) > len(pool_top & sep15) else "sep15",
        legacy_top49_overlap_oct03=len(sep15 & oct03),
    )
    table.sort("group", "mean_rank", descending=[False, True]).write_csv(
        out / "attribution_calls.csv"
    )
    pl.DataFrame(group_rows).write_csv(out / "attribution_groups.csv")
    return group_rows, attribution, identity, dict(sep15=sep15, oct03=oct03, pool_top=pool_top)


def part4(source, cfg, blind, sets, observed_delta, out):
    """Append the eighth transition to the frozen concordance table."""
    rng = np.random.default_rng(cfg["seed"] + 2)
    previous = pl.read_csv(source / cfg["inputs"]["transitions"])
    compound = blind["compound_id"].to_numpy()
    added_mask = np.isin(compound, list(sets["oct03"] - sets["sep15"]))
    removed_mask = np.isin(compound, list(sets["sep15"] - sets["oct03"]))

    rows = []
    for key in cfg["score_columns"]:
        values = blind[key].to_numpy()
        difference, p_value = permutation_difference(
            values[removed_mask], values[added_mask], cfg["permutations"], rng
        )
        rows.append(
            dict(
                transition="sep15->oct03",
                lb_delta=observed_delta,
                column=key,
                n_removed=int(removed_mask.sum()),
                n_added=int(added_mask.sum()),
                mean_removed=float(values[removed_mask].mean()),
                mean_added=float(values[added_mask].mean()),
                difference=difference,
                permutation_p=p_value,
                sign_agreement=bool(np.sign(difference) == np.sign(observed_delta)),
                selection_component=key in cfg["selection_components"],
            )
        )
    new = pl.DataFrame(rows)
    combined = pl.concat(
        [previous.select(new.columns[:-1]).with_columns(selection_component=pl.lit(False)), new]
    )
    tally = []
    for key in cfg["score_columns"]:
        part = combined.filter(pl.col("column") == key)
        tally.append(
            dict(
                column=key,
                transitions=part.height,
                sign_agreements=int(part["sign_agreement"].sum()),
                agreements_before=int(
                    part.filter(pl.col("transition") != "sep15->oct03")["sign_agreement"].sum()
                ),
                new_transition_agrees=bool(
                    part.filter(pl.col("transition") == "sep15->oct03")["sign_agreement"][0]
                ),
                significant=int((part["permutation_p"] < 0.05).sum()),
                selection_component=key in cfg["selection_components"],
            )
        )
    combined.write_csv(out / "transitions_updated.csv")
    pl.DataFrame(tally).write_csv(out / "concordance_tally.csv")
    return rows, tally


def build_report(
    scored,
    positives,
    estimators,
    fit_rows,
    band,
    alternatives,
    null_rows,
    bias_check,
    exchange_rows,
    scenario_rows,
    group_rows,
    attribution,
    concordance_rows,
    tally,
):
    return "\n".join(
        [
            "# Sizing the 2026-10-03 D6 Live gain",
            "",
            f"Scored subset estimate: {scored} of {BLIND_ROWS} rows, {positives} scored "
            f"positives (prevalence {positives / scored:.4f}).",
            "",
            "## Part 1 - operating-point estimators",
            "",
            markdown_table(
                estimators,
                [
                    "state",
                    "submitted",
                    "calls",
                    "mcc",
                    "positives_from_recall_sd",
                    "calls_from_precision_sd",
                    "scored_from_accuracy_sd",
                    "prevalence_identity",
                    "call_fraction_identity",
                    "scored_from_identity",
                    "positives_from_identity",
                ],
            ),
            "",
            "## Part 1 - joint integer fit",
            "",
            markdown_table(
                fit_rows,
                [
                    "state",
                    "total_calls",
                    "scored_calls",
                    "true_positive",
                    "false_positive",
                    "false_negative",
                    "fit_loss",
                    "mcc_reported",
                    "mcc_simulated",
                    "mcc_sd_reported",
                    "mcc_sd_simulated",
                    "precision_reported",
                    "precision_simulated",
                    "recall_reported",
                    "recall_simulated",
                ],
            ),
            "",
            "Alternative (scored, positives) pairs by total loss:",
            "",
            markdown_table(alternatives, list(alternatives[0].keys())),
            "",
            "## Part 1 - sensitivity band",
            "",
            markdown_table(band, list(band[0].keys())),
            "",
            "## Part 2 - unskilled swap null",
            "",
            markdown_table(
                null_rows,
                [
                    "transition",
                    "swapped_calls",
                    "observed_delta",
                    "observed_in_sd",
                    "null_mean",
                    "null_sd",
                    "null_p05",
                    "null_p50",
                    "null_p95",
                    "null_fraction_at_least",
                    "null_fraction_as_extreme",
                ],
            ),
            "",
            f"Bootstrap-mean versus point-estimate delta check: `{json.dumps(bias_check)}`",
            "",
            "## Part 2 - exchangeable split of the changed calls",
            "",
            markdown_table(
                exchange_rows,
                [
                    "transition",
                    "swapped_calls",
                    "observed_delta",
                    "scenarios",
                    "pooled_mean",
                    "pooled_sd",
                    "pooled_fraction_at_least",
                    "central_fraction_at_least",
                    "min_fraction_at_least",
                    "median_fraction_at_least",
                    "max_fraction_at_least",
                ],
            ),
            "",
            "Central scenario per transition: the neutral composition in which the retained,",
            "removed and added groups share one scored fraction and one precision among their",
            "scored calls. Every consistent scenario is in `exchangeable_scenarios.csv`.",
            "",
            markdown_table(
                [r for r in scenario_rows if r["central"]],
                [
                    "transition",
                    "changed_positives",
                    "changed_negatives",
                    "changed_unscored",
                    "removed_true_positive",
                    "added_true_positive",
                    "removed_false_positive",
                    "added_false_positive",
                    "fraction_at_least",
                ],
            ),
            "",
            "## Part 3 - attribution between the two score columns",
            "",
            markdown_table(group_rows, list(group_rows[0].keys())),
            "",
            f"Overlaps: `{json.dumps(attribution)}`",
            "",
            "## Part 4 - the eighth concordance transition",
            "",
            markdown_table(
                concordance_rows,
                [
                    "column",
                    "n_removed",
                    "n_added",
                    "mean_removed",
                    "mean_added",
                    "difference",
                    "permutation_p",
                    "sign_agreement",
                    "selection_component",
                ],
            ),
            "",
            markdown_table(tally, list(tally[0].keys())),
            "",
            "The deployed legacy score and `pool_lgbm` are the two components of the submitted",
            "rank mean, so their agreement on this transition is partly tautological.",
            "",
        ]
    )


def main():
    started = time.time()
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    states, fitted, estimators, fit_rows, band, alternatives, scored, positives = part1(
        source, cfg, out
    )
    null_rows, quantiles, bias_check, exchange_rows, scenario_rows = part2(cfg, states, fitted, out)
    blind = load_blind(source, cfg)
    group_rows, attribution, identity, sets = part3(source, cfg, blind, out)
    observed_delta = states["oct03"]["mcc"] - states["sep15"]["mcc"]
    concordance_rows, tally = part4(source, cfg, blind, sets, observed_delta, out)

    summary = dict(
        endpoint=ENDPOINT,
        runtime_seconds=round(time.time() - started, 2),
        scored_subset=scored,
        scored_positives=positives,
        scored_prevalence=positives / scored,
        sep15=fitted["sep15"],
        oct03=fitted["oct03"],
        true_positive_change=fitted["oct03"]["true_positive"] - fitted["sep15"]["true_positive"],
        false_positive_change=fitted["oct03"]["false_positive"] - fitted["sep15"]["false_positive"],
        scored_call_change=fitted["oct03"]["scored_calls"] - fitted["sep15"]["scored_calls"],
        observed_mcc_delta=observed_delta,
        null=null_rows,
        exchangeable_null=exchange_rows,
        bias_check=bias_check,
        attribution=attribution,
        identity_checks=identity,
        concordance=tally,
        status="Read-only post-mortem; no candidate, prepared submission or upload",
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    (out / "report.md").write_text(
        build_report(
            scored,
            positives,
            estimators,
            fit_rows,
            band,
            alternatives,
            null_rows,
            bias_check,
            exchange_rows,
            scenario_rows,
            group_rows,
            attribution,
            concordance_rows,
            tally,
        )
    )
    print(
        json.dumps(
            dict(
                scored=scored,
                positives=positives,
                true_positive_change=summary["true_positive_change"],
                false_positive_change=summary["false_positive_change"],
                runtime_seconds=summary["runtime_seconds"],
            )
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
