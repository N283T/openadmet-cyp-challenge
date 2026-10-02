"""Historical LB concordance of the D6 activity gate, legacy x gate OOF variants.

Part 1 scores the seven submitted D6 call sets and their chronological
transitions with blind score columns from run 7920b49a. For each transition it
compares the mean score of the removed calls with the mean of the added calls,
tests the difference by permutation, and records whether its sign agrees with
the observed Live MCC delta. The gate never saw blind labels, so this is a
frozen out-of-sample check of one score vector against the known LB directions.

Part 2 evaluates combinations of the legacy Sep 15 refit score with the gate on
the fixed umap/scaffold complete-pool OOF, with the legacy refit itself as the
control. Part 3 applies the same combinations to the deployed Sep 15 blind
score vector as a diagnostic swap inventory. Part 4 repeats Part 2 with the
pseudo-negative set restricted to primary log2fc >= 0.

Nothing here fits a model, registers a candidate or prepares a submission.
"""

import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl

ENDPOINT = "CYP2D6_is_TDI"
LEGACY = "legacy_sep15_refit"
CONTROL = "L0_legacy"
POOL = "pool_lgbm"
HISTORICAL_RATE = 49 / 750
DEPLOYED_BLIND_MODEL = "deployed_mean_block"


def percentile_rank(values):
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values))
    ranks[order] = np.arange(len(values))
    return (ranks + 0.5) / len(values)


def top_calls(score, fold, rate):
    """Call the top round(rate * n) molecules of each fold (stable compound order)."""
    calls = np.zeros(len(score), dtype=bool)
    for f in np.unique(fold):
        rows = np.flatnonzero(fold == f)
        k = int(round(rate * len(rows)))
        calls[rows[np.argsort(-score[rows], kind="stable")[:k]]] = True
    return calls


def mcc_from_counts(tp, fp, fn, tn):
    tp, fp, fn, tn = (np.asarray(v, dtype=float) for v in (tp, fp, fn, tn))
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return np.divide(tp * tn - fp * fn, denom, out=np.zeros_like(denom), where=denom > 0)


def confusion(y, calls):
    keep = np.isfinite(y)
    y, calls = y[keep].astype(bool), calls[keep]
    return (
        int((y & calls).sum()),
        int((~y & calls).sum()),
        int((y & ~calls).sum()),
        int((~y & ~calls).sum()),
    )


def mcc(y, calls):
    return float(mcc_from_counts(*confusion(y, calls)))


def group_bootstrap(y, groups, calls_model, calls_control, samples, seed):
    """Paired group bootstrap of pooled MCC(model) - MCC(control) on one label set."""
    keep = np.isfinite(y)
    y = y[keep].astype(bool)
    groups = groups[keep]
    _, index = np.unique(groups, return_inverse=True)
    count = index.max() + 1
    counts = []
    for calls in (calls_model[keep], calls_control[keep]):
        counts.append(
            [
                np.bincount(index, weights=(y & calls), minlength=count),
                np.bincount(index, weights=(~y & calls), minlength=count),
                np.bincount(index, weights=(y & ~calls), minlength=count),
                np.bincount(index, weights=(~y & ~calls), minlength=count),
            ]
        )
    draws = np.random.default_rng(seed).integers(0, count, (samples, count))
    weights = np.zeros((samples, count))
    for b in range(samples):
        weights[b] = np.bincount(draws[b], minlength=count)
    delta = mcc_from_counts(*(weights @ c for c in counts[0])) - mcc_from_counts(
        *(weights @ c for c in counts[1])
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


def fmt(value, digits=3):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "–"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def markdown_table(rows, columns):
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(c)) for c in columns) + " |")
    return "\n".join(lines)


def blind_columns(test_scores, keys):
    """Wide per-compound table of the requested `model.field` blind score columns."""
    out = None
    for key in keys:
        model, field = key.split(".")
        part = test_scores.filter(pl.col("model") == model).select(
            "compound_id", pl.col(field).alias(key)
        )
        assert part.height == 750 and part[key].null_count() == 0, key
        out = part if out is None else out.join(part, on="compound_id", validate="1:1")
    return out


def load_payloads(source, cfg, names):
    """Submitted D6 positive-call sets keyed by compound_id."""
    sets, rows = {}, []
    for date, alias in cfg["payload_inputs"].items():
        payload = pl.read_parquet(source / alias).join(
            names, left_on="Molecule_Name", right_on="molecule_name", how="left"
        )
        assert payload.height == 750 and payload["compound_id"].null_count() == 0, alias
        sets[date] = set(payload.filter(pl.col(ENDPOINT))["compound_id"].to_list())
        rows.append(dict(payload=date, alias=alias, positives=len(sets[date])))
    return sets, rows


def part1(source, cfg, names, out):
    """Transition concordance of blind score columns with the Live MCC deltas."""
    test_scores = pl.read_parquet(source / cfg["inputs"]["test_scores"])
    scores = blind_columns(test_scores, cfg["concordance_columns"])
    sets, payload_rows = load_payloads(source, cfg, names)
    keys = scores["compound_id"].to_list()
    values = {
        c: dict(zip(keys, scores[c].to_list(), strict=True)) for c in cfg["concordance_columns"]
    }

    gate = cfg["payload_gate_column"]
    for row in payload_rows:
        positives = sorted(sets[row["payload"]])
        row["live_mcc"] = cfg["live_mcc"][row["payload"]]
        for column in cfg["payload_summary_columns"]:
            series = np.array([values[column][c] for c in positives])
            row[f"mean_{column}"] = float(series.mean())
        row["n_below_0.3_" + gate] = int(
            (np.array([values[gate][c] for c in positives]) < 0.3).sum()
        )

    transition_rows = []
    for t_index, (source_date, target_date) in enumerate(cfg["transitions"]):
        removed = sorted(sets[source_date] - sets[target_date])
        added = sorted(sets[target_date] - sets[source_date])
        lb_delta = cfg["live_mcc"][target_date] - cfg["live_mcc"][source_date]
        for c_index, column in enumerate(cfg["concordance_columns"]):
            removed_values = np.array([values[column][c] for c in removed])
            added_values = np.array([values[column][c] for c in added])
            observed, p_value = permutation_difference(
                removed_values,
                added_values,
                cfg["permutations"],
                [cfg["seed"], t_index, c_index],
            )
            transition_rows.append(
                dict(
                    transition=f"{source_date}->{target_date}",
                    lb_mcc_from=cfg["live_mcc"][source_date],
                    lb_mcc_to=cfg["live_mcc"][target_date],
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
    payloads = pl.DataFrame(payload_rows, infer_schema_length=None)
    payloads.write_csv(out / "payload_summary.csv")
    agreement = (
        transitions.group_by("column")
        .agg(
            pl.col("sign_agreement").sum().alias("agreeing_transitions"),
            pl.len().alias("transitions"),
            pl.col("permutation_p").min().alias("min_permutation_p"),
            (pl.col("permutation_p") < 0.05).sum().alias("significant_transitions"),
        )
        .sort("agreeing_transitions", descending=True)
    )
    agreement.write_csv(out / "concordance_agreement.csv")
    return transition_rows, payload_rows, agreement.to_dicts(), sets, scores


def variant_scores(frame, taus):
    """Per-fold percentile-rank combinations of the legacy score with the gate."""
    fold = frame["fold"].to_numpy()
    legacy = frame["legacy"].to_numpy()
    pool = frame["pool_score"].to_numpy()
    pa_m = frame["p_active_m"].to_numpy()
    pa_c = frame["p_active_c"].to_numpy()
    legacy_rank, pool_rank = np.empty(len(legacy)), np.empty(len(legacy))
    for f in np.unique(fold):
        rows = np.flatnonzero(fold == f)
        legacy_rank[rows] = percentile_rank(legacy[rows])
        pool_rank[rows] = percentile_rank(pool[rows])
    scores = {
        CONTROL: legacy_rank,
        "L1_legacy_x_pa_morgan": legacy_rank * pa_m,
        "L2_legacy_x_sqrt_pa_morgan": legacy_rank * np.sqrt(pa_m),
        "L3_legacy_x_sqrt_pa_both": legacy_rank * np.sqrt(pa_m * pa_c),
        "L4_rankmean_legacy_pool": 0.5 * (legacy_rank + pool_rank),
    }
    for tau in taus:
        scores[f"L5_gate_filter_tau{tau:g}"] = legacy_rank + 10.0 * (pa_m >= tau)
    scores["B_pool_lgbm"] = pool_rank
    return scores, legacy_rank, pool_rank


def decision_row(y_official, y_testlike, calls, **keys):
    tp_o, fp_o, fn_o, tn_o = confusion(y_official, calls)
    tp_t, fp_t, fn_t, tn_t = confusion(y_testlike, calls)
    return dict(
        **keys,
        n_molecules=len(calls),
        n_calls=int(calls.sum()),
        mcc_official=float(mcc_from_counts(tp_o, fp_o, fn_o, tn_o)),
        precision_official=tp_o / (tp_o + fp_o) if tp_o + fp_o else None,
        tp_official=tp_o,
        fp_official=fp_o,
        mcc_testlike=float(mcc_from_counts(tp_t, fp_t, fn_t, tn_t)),
        precision_testlike=tp_t / (tp_t + fp_t) if tp_t + fp_t else None,
        tp_testlike=tp_t,
        fp_testlike=fp_t,
        recall_testlike=tp_t / (tp_t + fn_t) if tp_t + fn_t else None,
    )


def evaluate_variants(frame, scores, testlike_column, label_variant, split, cfg, rates):
    """Rate curve, nested rate selection and paired bootstrap against the legacy control."""
    fold = frame["fold"].to_numpy()
    groups = frame["group"].to_numpy().astype(str)
    y_official = frame["y_official"].to_numpy()
    y_testlike = frame[testlike_column].to_numpy()
    curve_rows, nested_rows, boot_rows = [], [], []
    calls_at = {}
    for model, score in scores.items():
        by_rate = {}
        for rate in rates:
            calls = top_calls(score, fold, rate)
            by_rate[rate] = calls
            keys = dict(
                split_id=split, label_variant=label_variant, model=model, rate=rate, fold=-1
            )
            curve_rows.append(decision_row(y_official, y_testlike, calls, **keys))
            for f in np.unique(fold):
                mask = fold == f
                curve_rows.append(
                    decision_row(
                        y_official[mask],
                        y_testlike[mask],
                        calls[mask],
                        **{**keys, "fold": int(f)},
                    )
                )
        selected = np.zeros(len(score), dtype=bool)
        for f in np.unique(fold):
            other = fold != f
            inner = {r: mcc(y_testlike[other], by_rate[r][other]) for r in rates}
            rate = max(rates, key=lambda r: (inner[r], -r))
            held = fold == f
            selected[held] = by_rate[rate][held]
            nested_rows.append(
                decision_row(
                    y_official[held],
                    y_testlike[held],
                    by_rate[rate][held],
                    split_id=split,
                    label_variant=label_variant,
                    model=model,
                    fold=int(f),
                    selected_rate=rate,
                    selection_mcc_testlike=inner[rate],
                )
            )
        nested_rows.append(
            decision_row(
                y_official,
                y_testlike,
                selected,
                split_id=split,
                label_variant=label_variant,
                model=model,
                fold=-1,
                selected_rate=None,
                selection_mcc_testlike=None,
            )
        )
        calls_at[model] = dict(nested=selected, historical=by_rate[HISTORICAL_RATE])
    for model in scores:
        if model == CONTROL:
            continue
        for operating in ["nested", "historical"]:
            for label, y in [("testlike", y_testlike), ("official", y_official)]:
                boot_rows.append(
                    dict(
                        split_id=split,
                        label_variant=label_variant,
                        model=model,
                        control=CONTROL,
                        operating_point=operating,
                        label_set=label,
                        **group_bootstrap(
                            y,
                            groups,
                            calls_at[model][operating],
                            calls_at[CONTROL][operating],
                            cfg["bootstrap_samples"],
                            cfg["seed"],
                        ),
                    )
                )
    return curve_rows, nested_rows, boot_rows


def part3(source, cfg, names, sep15_positive, blind, out):
    """Blind swap inventory of the legacy x gate variants (diagnostic only)."""
    deployed = pl.read_parquet(source / cfg["inputs"]["sep15_blind_scores"]).filter(
        pl.col("model") == DEPLOYED_BLIND_MODEL
    )
    assert deployed.height == 750
    ordered = deployed.sort(["score", "compound_id"], descending=[True, False])
    top49 = set(ordered["compound_id"].to_list()[:49])
    verification = dict(
        deployed_rows=deployed.height,
        deployed_calls=int(deployed["call"].sum()),
        top49_equals_sep15_payload=bool(top49 == sep15_positive),
        boundary_gap=float(ordered["score"].to_numpy()[48] - ordered["score"].to_numpy()[49]),
    )
    verification.update(
        blind_gate_ge_tau=int((blind[cfg["blind_gate_morgan_column"]] >= cfg["blind_tau"]).sum()),
        blind_tau=cfg["blind_tau"],
    )
    assert verification["top49_equals_sep15_payload"], verification

    table = (
        deployed.select("compound_id", pl.col("score").alias("legacy_score"))
        .join(blind, on="compound_id", validate="1:1")
        .join(names, on="compound_id", how="left")
        .sort("compound_id")
    )
    pool_score = table[cfg["blind_pool_column"]].to_numpy()
    legacy_rank = percentile_rank(table["legacy_score"].to_numpy())
    pool_rank = percentile_rank(pool_score)
    pa_m = table[cfg["blind_gate_morgan_column"]].to_numpy()
    pa_c = table[cfg["blind_gate_chemeleon_column"]].to_numpy()
    tau = cfg["blind_tau"]
    variants = {
        CONTROL: legacy_rank,
        "L1_legacy_x_pa_morgan": legacy_rank * pa_m,
        "L2_legacy_x_sqrt_pa_morgan": legacy_rank * np.sqrt(pa_m),
        "L3_legacy_x_sqrt_pa_both": legacy_rank * np.sqrt(pa_m * pa_c),
        "L4_rankmean_legacy_pool": 0.5 * (legacy_rank + pool_rank),
        f"L5_gate_filter_tau{tau:g}": legacy_rank + 10.0 * (pa_m >= tau),
    }
    ids = table["compound_id"].to_list()
    molecule = table["molecule_name"].to_list()
    is_sep15 = np.array([c in sep15_positive for c in ids])
    summary_rows, swap_rows, candidate_rows = [], [], []
    for variant, score in variants.items():
        for n_calls in cfg["blind_call_counts"]:
            order = np.argsort(-score, kind="stable")[:n_calls]
            calls = np.zeros(len(score), dtype=bool)
            calls[order] = True
            added = np.flatnonzero(calls & ~is_sep15)
            removed = np.flatnonzero(~calls & is_sep15)
            summary_rows.append(
                dict(
                    variant=variant,
                    n_calls=n_calls,
                    overlap_sep15=int((calls & is_sep15).sum()),
                    added=len(added),
                    removed=len(removed),
                    mean_pa_morgan_called=float(pa_m[calls].mean()),
                    mean_pa_morgan_added=float(pa_m[added].mean()) if len(added) else None,
                    mean_pa_morgan_removed=float(pa_m[removed].mean()) if len(removed) else None,
                )
            )
            candidate_rows.append(
                pl.DataFrame(
                    dict(
                        molecule_name=molecule,
                        variant=[variant] * len(ids),
                        n_calls=[n_calls] * len(ids),
                        call=calls,
                    )
                )
            )
            for action, rows in [("added", added), ("removed", removed)]:
                for i in rows:
                    swap_rows.append(
                        dict(
                            variant=variant,
                            n_calls=n_calls,
                            molecule_name=molecule[i],
                            action=action,
                            legacy_rank=float(legacy_rank[i]),
                            p_active_morgan=float(pa_m[i]),
                            p_active_chemeleon=float(pa_c[i]),
                            pool_score=float(pool_score[i]),
                            variant_score=float(score[i]),
                        )
                    )
    pl.concat(candidate_rows).write_parquet(out / "blind_candidates.parquet")
    swaps = pl.DataFrame(swap_rows)
    swaps.write_csv(out / "swaps.csv")
    pl.DataFrame(summary_rows, infer_schema_length=None).write_csv(out / "blind_swap_summary.csv")
    sep15_gate = dict(
        n_positives=int(is_sep15.sum()),
        mean_p_active_morgan=float(pa_m[is_sep15].mean()),
        mean_p_active_chemeleon=float(pa_c[is_sep15].mean()),
        **{f"below_{t}": int((pa_m[is_sep15] < t).sum()) for t in cfg["sep15_gate_thresholds"]},
    )
    return verification, summary_rows, sep15_gate, swaps


def main():
    started = time.time()
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    inputs = cfg["inputs"]
    rates = sorted({*cfg["rates"], HISTORICAL_RATE})
    names = (
        pl.read_parquet(source / inputs["test_mapping"])
        .select("compound_id", "molecule_name")
        .unique()
    )
    assert names.height == 750

    transition_rows, payload_rows, agreement_rows, payload_sets, blind = part1(
        source, cfg, names, out
    )

    # Development OOF: one wide row per compound and split family.
    oof = pl.read_parquet(source / inputs["oof_scores"])
    assert oof["endpoint"].unique().to_list() == [ENDPOINT]
    assignments = pl.read_parquet(source / inputs["assignments"])
    aliases = pl.read_parquet(source / inputs["aliases"])
    primary = (
        pl.read_csv(source / inputs["primary_screen"])
        .rename({"Molecule_Name": "molecule_name"})
        .filter(pl.col("enzyme") == cfg["primary_enzyme"])
        .join(
            aliases.select("molecule_name", "compound_id").unique(),
            on="molecule_name",
            how="left",
            validate="m:1",
        )
    )
    assert primary["compound_id"].null_count() == 0
    log2fc = primary.group_by("compound_id").agg(pl.col("log2fc_estimate").median().alias("log2fc"))

    def wide(split):
        rows = oof.filter(pl.col("split_id") == split)

        def column(model, field, alias):
            part = rows.filter(pl.col("model") == model).select(
                "compound_id", pl.col(field).alias(alias)
            )
            assert part.height == 6146 and part[alias].null_count() == 0, (model, field)
            return part

        frame = (
            rows.filter(pl.col("model") == LEGACY)
            .select(
                "compound_id",
                "fold",
                "pool_role",
                "y_official",
                "y_testlike",
                pl.col("score").alias("legacy"),
            )
            .join(column("gate_x_tdi", "p_active", "p_active_m"), on="compound_id", validate="1:1")
            .join(
                column("gate_x_tdi_chemeleon", "p_active", "p_active_c"),
                on="compound_id",
                validate="1:1",
            )
            .join(column(POOL, "score", "pool_score"), on="compound_id", validate="1:1")
            .join(
                assignments.filter(pl.col("strategy") == split).select(
                    "compound_id", "group", pl.col("fold").alias("split_fold")
                ),
                on="compound_id",
                how="left",
                validate="1:1",
            )
            .join(log2fc, on="compound_id", how="left")
            .sort("fold", "compound_id")
        )
        assert frame["group"].null_count() == 0
        assert (frame["fold"] == frame["split_fold"]).all()
        pseudo = frame["pool_role"] == "pseudo_negative"
        assert frame.filter(pseudo)["log2fc"].null_count() == 0
        strict = (
            pl.when(pseudo & (pl.col("log2fc") < cfg["sensitivity_min_log2fc"]))
            .then(None)
            .otherwise(pl.col("y_testlike"))
            .alias("y_testlike_strict")
        )
        return frame.with_columns(strict).with_columns(
            pl.col(["y_official", "y_testlike", "y_testlike_strict"]).fill_null(np.nan)
        )

    smoke = cfg.get("smoke")
    curve_rows, nested_rows, boot_rows, pool_rows = [], [], [], []
    for split in cfg["split_strategies"]:
        frame = wide(split)
        if smoke:
            rng = np.random.default_rng(cfg["seed"])
            keep = rng.choice(frame.height, smoke["development_sample"], replace=False)
            frame = frame[np.sort(keep)]
        scores, legacy_rank, _ = variant_scores(frame, cfg["gate_filter_taus"])
        pool_rows.append(
            dict(
                split_id=split,
                molecules=frame.height,
                official=int(np.isfinite(frame["y_official"].to_numpy()).sum()),
                testlike=int(np.isfinite(frame["y_testlike"].to_numpy()).sum()),
                testlike_strict=int(np.isfinite(frame["y_testlike_strict"].to_numpy()).sum()),
                legacy_gate_spearman=float(
                    frame.select(pl.corr("legacy", "p_active_m", method="spearman")).item()
                ),
                gate_ge_0_2=int((frame["p_active_m"] >= 0.2).sum()),
                gate_ge_0_3=int((frame["p_active_m"] >= 0.3).sum()),
                gate_ge_0_4=int((frame["p_active_m"] >= 0.4).sum()),
                legacy_rank_mean=float(legacy_rank.mean()),
            )
        )
        for label_variant, column in [
            ("testlike", "y_testlike"),
            ("testlike_strict", "y_testlike_strict"),
        ]:
            c, n, b = evaluate_variants(frame, scores, column, label_variant, split, cfg, rates)
            curve_rows += c
            nested_rows += n
            boot_rows += b

    curves = pl.DataFrame(curve_rows, infer_schema_length=None)
    nested = pl.DataFrame(nested_rows, infer_schema_length=None)
    bootstrap = pl.DataFrame(boot_rows, infer_schema_length=None)
    curves.write_csv(out / "oof_rate_curves.csv")
    nested.write_csv(out / "oof_nested_selection.csv")
    bootstrap.write_csv(out / "oof_bootstrap_deltas.csv")

    verification, swap_summary, sep15_gate, swaps = part3(
        source, cfg, names, payload_sets["sep15"], blind, out
    )

    # Sensitivity: does the variant ordering change with the stricter pseudo-negatives?
    def ordering(label_variant, split, metric="mcc_testlike"):
        rows = nested.filter(
            (pl.col("label_variant") == label_variant)
            & (pl.col("split_id") == split)
            & (pl.col("fold") == -1)
        ).sort(metric, descending=True)
        return rows["model"].to_list()

    sensitivity = []
    for split in cfg["split_strategies"]:
        base, strict = ordering("testlike", split), ordering("testlike_strict", split)
        sensitivity.append(
            dict(
                split_id=split,
                testlike_order=" > ".join(base),
                testlike_strict_order=" > ".join(strict),
                identical=base == strict,
                top_model_changed=base[0] != strict[0],
                spearman_of_ranks=float(
                    pl.DataFrame(
                        {"a": list(range(len(base))), "b": [strict.index(m) for m in base]}
                    )
                    .select(pl.corr("a", "b", method="spearman"))
                    .item()
                ),
            )
        )

    summary = dict(
        config=cfg,
        endpoint=ENDPOINT,
        smoke=bool(smoke),
        rates=rates,
        payloads=payload_rows,
        transitions=transition_rows,
        concordance_agreement=agreement_rows,
        oof_pool=pool_rows,
        oof_nested_pooled=[r for r in nested_rows if r["fold"] == -1],
        oof_nested_per_fold=[r for r in nested_rows if r["fold"] >= 0],
        oof_bootstrap=boot_rows,
        blind_verification=verification,
        blind_swap_summary=swap_summary,
        sep15_gate=sep15_gate,
        sensitivity=sensitivity,
        runtime_seconds=round(time.time() - started, 1),
        caveats=[
            "Seven transitions share one hidden blind set; they are not independent trials",
            "The public Live MCC is itself an estimate on an unknown scored subset",
            "Pseudo-negative labels are inferred from the primary screen, not measured",
            "Blind variant vectors are diagnostic; no candidate or submission was prepared",
        ],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")

    boot_columns = [
        "split_id",
        "model",
        "operating_point",
        "label_set",
        "delta",
        "ci_low",
        "ci_high",
        "p_delta_gt_0",
    ]
    boot_main = [r for r in boot_rows if r["label_variant"] == "testlike"]
    boot_strict = [r for r in boot_rows if r["label_variant"] == "testlike_strict"]
    report = [
        "# D6 activity-gate historical concordance and legacy x gate OOF variants",
        "",
        f"Smoke: {bool(smoke)}. Rate grid: " + ", ".join(f"{100 * r:.2f}%" for r in rates) + ".",
        "",
        "## Part 1 - submitted payloads",
        "",
        markdown_table(
            payload_rows,
            [
                "payload",
                "positives",
                "live_mcc",
                *[f"mean_{c}" for c in cfg["payload_summary_columns"]],
                f"n_below_0.3_{cfg['payload_gate_column']}",
            ],
        ),
        "",
        "## Part 1 - transition concordance",
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
        "## Part 1 - sign agreement by score column",
        "",
        markdown_table(agreement_rows, list(agreement_rows[0].keys())),
        "",
        "## Part 2 - OOF pools",
        "",
        markdown_table(pool_rows, list(pool_rows[0].keys())),
        "",
        "## Part 2 - nested selection (pooled held folds)",
        "",
        markdown_table(
            [r for r in nested_rows if r["fold"] == -1],
            [
                "split_id",
                "label_variant",
                "model",
                "n_calls",
                "mcc_testlike",
                "precision_testlike",
                "tp_testlike",
                "fp_testlike",
                "mcc_official",
                "precision_official",
                "tp_official",
                "fp_official",
            ],
        ),
        "",
        "## Part 2 - paired group bootstrap versus the legacy control",
        "",
        markdown_table(boot_main, boot_columns),
        "",
        "## Part 3 - blind swaps against the Sep 15 vector",
        "",
        f"Verification: `{json.dumps(verification)}`",
        "",
        markdown_table(swap_summary, list(swap_summary[0].keys())),
        "",
        f"Sep 15 gate profile: `{json.dumps(sep15_gate)}`",
        "",
        "## Part 4 - sensitivity to the stricter pseudo-negative set",
        "",
        markdown_table(sensitivity, list(sensitivity[0].keys())),
        "",
        "Paired group bootstrap under the stricter pseudo-negative set:",
        "",
        markdown_table(boot_strict, boot_columns),
        "",
        "Transitions are few, share one hidden label set and one bootstrap LB estimate;",
        "sign agreement here is calibration evidence about the score column, not proof",
        "that any variant would improve the Live MCC.",
        "",
    ]
    (out / "report.md").write_text("\n".join(report))
    print(
        json.dumps(
            dict(
                runtime_seconds=summary["runtime_seconds"],
                transitions=len(transition_rows),
                swaps=swaps.height,
                oof_rows=curves.height,
            )
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
