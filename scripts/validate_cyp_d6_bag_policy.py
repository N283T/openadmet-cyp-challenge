"""Fully nested candidate-policy validation; no outer labels enter a fit or cut."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from evaluate_cyp_d6_zero_start import bootstrap_mcc, fit_model, folds, predict
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, brier_score_loss, matthews_corrcoef


def select_cut(y, scores, samples):
    observed = np.isfinite(y)
    cuts = np.r_[-np.inf, np.unique(np.quantile(scores[observed], np.linspace(0, 1, 301))), np.inf]
    values, _ = bootstrap_mcc(y[observed], scores[observed], cuts, samples, 0)
    return float(cuts[int(np.argmax(values))])


def select_prior_cuts(y, scores, priors, samples, include_mean=False):
    """Train-only label-shift scenarios; each prior selects a cut, never a call budget."""
    observed = np.isfinite(y)
    truth, score = y[observed].astype(int), scores[observed]
    order = np.argsort(score, kind="stable")
    truth, score = truth[order], score[order]
    thresholds = np.r_[-np.inf, np.unique(np.quantile(score, np.linspace(0, 1, 301))), np.inf]
    cut_indices = np.searchsorted(score, thresholds, side="right")
    draws = np.random.default_rng(0).integers(0, len(truth), (samples, len(truth)))
    weights = np.zeros((samples, len(truth)), dtype=np.int32)
    np.add.at(weights, (np.arange(samples)[:, None], draws), 1)
    cp = np.pad(np.cumsum(weights * truth, axis=1), ((0, 0), (1, 0)))
    cn = np.pad(np.cumsum(weights * (1 - truth), axis=1), ((0, 0), (1, 0)))
    tpr = (cp[:, -1, None] - cp[:, cut_indices]) / np.maximum(cp[:, -1, None], 1)
    fpr = (cn[:, -1, None] - cn[:, cut_indices]) / np.maximum(cn[:, -1, None], 1)
    result = {}
    utilities = []
    for prior in priors:
        tp, fn, fp, tn = prior * tpr, prior * (1 - tpr), (1 - prior) * fpr, (1 - prior) * (1 - fpr)
        denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        value = np.divide(tp * tn - fp * fn, denom, out=np.zeros_like(denom), where=denom > 0)
        result[prior] = float(thresholds[np.argmax(value.mean(axis=0))])
        utilities.append(value.mean(axis=0))
    if include_mean:
        result["mean"] = float(thresholds[np.argmax(np.mean(utilities, axis=0))])
    return result


def metric(y, score, calls):
    observed = np.isfinite(y)
    truth = y[observed].astype(int)
    return dict(
        ap=float(average_precision_score(truth, score[observed])),
        mcc=float(matthews_corrcoef(truth, calls[observed])),
        full_pool_rate=float(calls.mean()),
        labelled_rate=float(calls[observed].mean()),
        precision=float(truth[calls[observed]].mean()) if calls[observed].any() else None,
        recall=float(calls[observed][truth == 1].mean()),
        n_labelled=int(observed.sum()),
    )


def evaluate(table, fp, x, y, tdi, shift, weights, subtype, tasks, development, blind, cfg):
    import torch

    torch.set_num_threads(cfg["cpu_threads"])
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    assignments = pl.read_parquet(inp / cfg["fixed_assignments"])
    rows, cuts_audit, leakage, fold_metrics = [], [], [], []
    selected = np.ones(len(y), dtype=bool)
    if cfg.get("selection_control"):
        primary = tasks[:, 16 + cfg["primary_enzyme_index"]]
        selected = np.isfinite(primary) & (primary < cfg["primary_activity_cut"])
        pl.DataFrame(
            {
                "compound_id": table["compound_id"],
                "endpoint": cfg["endpoint"],
                "label": y,
                "primary_log2fc": primary,
                "selected": selected,
            }
        ).write_parquet(out / "selection-control-source.parquet")
    partition_y = np.where(selected, y, np.nan)
    for strategy in cfg["split_strategies"]:
        ordered = table.select("compound_id").join(
            assignments.filter(pl.col("strategy") == strategy).select(
                "compound_id", "fold", "group"
            ),
            on="compound_id",
            how="left",
            maintain_order="left",
        )
        membership = ordered["fold"].to_numpy()
        groups = ordered["group"].fill_null("blind").to_numpy()
        assert np.isfinite(membership[development]).all()
        for seed in cfg["seeds"]:
            for outer in range(cfg["outer_folds"]):
                train = development[membership[development] != outer]
                test = development[membership[development] == outer]
                assert set(groups[train]).isdisjoint(groups[test])
                bag_splits = folds(
                    train, partition_y, groups, strategy, cfg["bag_folds"], seed + outer * 100
                )
                for kind in cfg["models"]:
                    fit_y = (
                        np.where(selected, y, np.nan)
                        if kind
                        in {
                            "activity_selected",
                            "gate_selected",
                            "chemeleon_selected_network",
                            "chemeleon_selected_primary_network",
                        }
                        else y
                    )
                    test_scores, test_margins, test_votes = [], [], []
                    rank_scores, rank_margins = [], []
                    bag_oof = np.full(len(table), np.nan)
                    inner_cuts = []
                    prior_margins = {p: [] for p in cfg.get("threshold_priors", [])}
                    if cfg.get("mean_prior_policy"):
                        prior_margins["mean"] = []
                    for bag, (bag_train, bag_valid) in enumerate(bag_splits):
                        assert set(groups[bag_train]).isdisjoint(groups[bag_valid])
                        inner_score = np.full(len(table), np.nan)
                        inner_splits = folds(
                            bag_train,
                            partition_y,
                            groups,
                            strategy,
                            cfg["inner_folds"],
                            seed + outer * 100 + bag + 1,
                        )
                        for inner, (fit, valid) in enumerate(inner_splits):
                            assert set(groups[fit]).isdisjoint(groups[valid])
                            assert not np.isin(np.r_[fit, valid], test).any()
                            fitted = fit_model(
                                kind,
                                fit,
                                fp,
                                x,
                                fit_y,
                                tdi,
                                shift,
                                seed + outer * 1000 + bag * 10 + inner,
                                weights,
                                subtype,
                                tasks,
                                cfg,
                            )
                            inner_score[valid] = predict(kind, fitted, valid, fp, x)
                        assert np.isfinite(inner_score[bag_train]).all()
                        cut = select_cut(
                            fit_y[bag_train], inner_score[bag_train], cfg["bootstrap_samples"]
                        )
                        fitted = fit_model(
                            kind,
                            bag_train,
                            fp,
                            x,
                            fit_y,
                            tdi,
                            shift,
                            seed + outer * 1000 + bag * 10 + 9,
                            weights,
                            subtype,
                            tasks,
                            cfg,
                        )
                        scores = predict(kind, fitted, test, fp, x)
                        bag_oof[bag_valid] = predict(kind, fitted, bag_valid, fp, x)
                        test_scores.append(scores)
                        test_margins.append(scores - cut)
                        test_votes.append(scores > cut)
                        prior_cuts = (
                            select_prior_cuts(
                                fit_y[bag_train],
                                inner_score[bag_train],
                                cfg.get("threshold_priors", []),
                                cfg["bootstrap_samples"],
                                cfg.get("mean_prior_policy", False),
                            )
                            if prior_margins
                            else {}
                        )
                        for prior, prior_cut in prior_cuts.items():
                            prior_margins[prior].append(scores - prior_cut)
                        if cfg.get("rank_transfer_ablation"):
                            inner_percentile = float((inner_score[bag_train] <= cut).mean())
                            percentile = rankdata(scores, method="average") / (len(scores) + 1)
                            rank_scores.append(percentile)
                            rank_margins.append(percentile - inner_percentile)
                        inner_cuts.append(cut)
                        labelled_train = bag_train[np.isfinite(fit_y[bag_train])]
                        labelled_valid = bag_valid[np.isfinite(fit_y[bag_valid])]
                        cuts_audit.append(
                            dict(
                                split_id=strategy,
                                seed=seed,
                                outer_fold=outer,
                                model=kind,
                                bag_fold=bag,
                                threshold=cut,
                                **{f"prior_{p}_cut": c for p, c in prior_cuts.items()},
                                train_labels=len(labelled_train),
                                inner_score_q10=float(
                                    np.quantile(inner_score[labelled_train], 0.1)
                                ),
                                inner_score_q50=float(
                                    np.quantile(inner_score[labelled_train], 0.5)
                                ),
                                inner_score_q90=float(
                                    np.quantile(inner_score[labelled_train], 0.9)
                                ),
                                refit_valid_q10=float(np.quantile(bag_oof[labelled_valid], 0.1)),
                                refit_valid_q50=float(np.quantile(bag_oof[labelled_valid], 0.5)),
                                refit_valid_q90=float(np.quantile(bag_oof[labelled_valid], 0.9)),
                            )
                        )
                    assert np.isfinite(bag_oof[train]).all()
                    scores_mean = np.mean(test_scores, axis=0)
                    margin_mean = np.mean(test_margins, axis=0)
                    pooled_cut = select_cut(fit_y[train], bag_oof[train], cfg["bootstrap_samples"])
                    # Independent full-refit control: tune a cut only on outer-train inner OOF.
                    single_inner = np.full(len(table), np.nan)
                    for number, (fit, valid) in enumerate(
                        folds(
                            train,
                            partition_y,
                            groups,
                            strategy,
                            cfg["inner_folds"],
                            seed + 100 + outer,
                        )
                    ):
                        assert set(groups[fit]).isdisjoint(groups[valid])
                        fitted = fit_model(
                            kind,
                            fit,
                            fp,
                            x,
                            fit_y,
                            tdi,
                            shift,
                            seed + outer * 10 + number,
                            weights,
                            subtype,
                            tasks,
                            cfg,
                        )
                        single_inner[valid] = predict(kind, fitted, valid, fp, x)
                    single_cut = select_cut(
                        fit_y[train], single_inner[train], cfg["bootstrap_samples"]
                    )
                    fitted = fit_model(
                        kind,
                        train,
                        fp,
                        x,
                        fit_y,
                        tdi,
                        shift,
                        seed + outer,
                        weights,
                        subtype,
                        tasks,
                        cfg,
                    )
                    single_score = predict(kind, fitted, test, fp, x)
                    policies = {
                        "centered_fold_mean": (scores_mean, margin_mean > 0),
                        "majority_vote": (scores_mean, np.sum(test_votes, axis=0) >= 3),
                        "pooled_oof_cut": (scores_mean, scores_mean > pooled_cut),
                        "single_refit": (single_score, single_score > single_cut),
                    }
                    if cfg.get("rank_transfer_ablation"):
                        policies["rank_margin_mean"] = (
                            np.mean(rank_scores, axis=0),
                            np.mean(rank_margins, axis=0) > 0,
                        )
                    for prior, margins in prior_margins.items():
                        policies[f"prior_{prior}_centered"] = (
                            scores_mean,
                            np.mean(margins, axis=0) > 0,
                        )
                    for policy, (score, calls) in policies.items():
                        entry = dict(
                            split_id=strategy, seed=seed, fold=outer, model=kind, policy=policy
                        )
                        fold_metrics.append({**entry, **metric(y[test], score, calls)})
                        for i, index in enumerate(test):
                            rows.append(
                                {
                                    **entry,
                                    "compound_id": table["compound_id"][int(index)],
                                    "molecule_name": table["molecule_name"][int(index)],
                                    "endpoint": cfg["endpoint"],
                                    "group": groups[index],
                                    "label": float(y[index]),
                                    "selected": bool(selected[index]),
                                    "score": float(score[i]),
                                    "call": bool(calls[i]),
                                }
                            )
                    leakage.append(
                        dict(
                            split_id=strategy,
                            seed=seed,
                            fold=outer,
                            model=kind,
                            outer_test_used_in_fit_or_cut=False,
                            mean_bag_threshold=float(np.mean(inner_cuts)),
                            pooled_cut=pooled_cut,
                            single_cut=single_cut,
                            outer_labels=int(np.isfinite(y[test]).sum()),
                        )
                    )
                    # Checkpoint keyed outputs after each completed block for recoverability.
                    pl.DataFrame(rows).write_parquet(out / "nested-oof.parquet")
                    pl.DataFrame(cuts_audit).write_csv(out / "threshold-scale-audit.csv")
                    pl.DataFrame(leakage).write_csv(out / "leakage-and-operating-points.csv")
                    print(
                        json.dumps(
                            {
                                **entry,
                                "policies": {
                                    p: metric(y[test], s, c) for p, (s, c) in policies.items()
                                },
                            }
                        ),
                        flush=True,
                    )
    frame = pl.DataFrame(rows)
    summaries = []
    for keys, sub in frame.group_by("split_id", "seed", "model", "policy", maintain_order=True):
        summary = dict(zip(["split_id", "seed", "model", "policy"], keys, strict=True))
        summary.update(
            metric(sub["label"].to_numpy(), sub["score"].to_numpy(), sub["call"].to_numpy())
        )
        if (
            summary["model"] != "historical_morgan_rbf_c1"
            and summary["policy"] != "rank_margin_mean"
        ):
            observed = sub.filter(pl.col("label").is_finite())
            summary["brier"] = float(brier_score_loss(observed["label"], observed["score"]))
        summaries.append(summary)
    pl.DataFrame(summaries).write_csv(out / "aggregate-metrics.csv")
    pl.DataFrame(fold_metrics).write_csv(out / "fold-metrics.csv")
    if not cfg.get("skip_historical_diagnostics"):
        historical_diagnostics(frame, inp, out, cfg)
    if cfg.get("selection_control"):
        selection_summaries = []
        for keys, sub in frame.filter(pl.col("selected")).group_by(
            "split_id", "seed", "model", "policy", maintain_order=True
        ):
            entry = dict(zip(["split_id", "seed", "model", "policy"], keys, strict=True))
            selection_summaries.append(
                {
                    **entry,
                    **metric(
                        sub["label"].to_numpy(), sub["score"].to_numpy(), sub["call"].to_numpy()
                    ),
                }
            )
        pl.DataFrame(selection_summaries).write_csv(out / "selected-only-metrics.csv")
    paired_intervals(frame, out, cfg)
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "aggregate": summaries,
                "no_upload": True,
                "blind_labels_used": False,
                "leakage_checks_passed": True,
                "limits": cfg["limitations"],
            },
            indent=2,
        )
        + "\n"
    )


def historical_diagnostics(frame, inp, out, cfg):
    from review_cyp_d6_causal_probe import calls_at_rate, projected

    curves = []
    rates = [0.02, 0.05, 0.06533333333333333, 0.1, 0.15, 0.2, 0.3]
    for split in cfg["split_strategies"]:
        source = frame.filter(
            (pl.col("split_id") == split) & (pl.col("policy") == "centered_fold_mean")
        )
        for seed in cfg["seeds"]:
            modern = source.filter(pl.col("seed") == seed)
            for alias in [f"control-{split}.parquet", f"ap-{split}.parquet"]:
                hist = pl.read_parquet(inp / alias)
                for kind in hist["model"].unique():
                    old = hist.filter(pl.col("model") == kind)
                    modern = pl.concat(
                        [
                            modern,
                            old.select(
                                "compound_id",
                                "fold",
                                "group",
                                "split_id",
                                "model",
                                pl.col("y_true").cast(pl.Float64).alias("label"),
                                "score",
                            )
                            .with_columns(
                                pl.lit(seed).alias("seed"),
                                pl.lit("centered_fold_mean").alias("policy"),
                                pl.lit(cfg["endpoint"]).alias("endpoint"),
                                pl.lit(None, dtype=pl.String).alias("molecule_name"),
                                pl.lit(False).alias("call"),
                            )
                            .select(modern.columns),
                        ],
                        how="vertical_relaxed",
                    )
            for kind in modern["model"].unique():
                sub = modern.filter(pl.col("model") == kind).sort("compound_id")
                reference = (
                    source.filter(pl.col("seed") == seed)
                    .filter(pl.col("model") == "fresh_d6_network")
                    .sort("compound_id")
                )
                assert sub["compound_id"].equals(reference["compound_id"])
                assert sub["fold"].equals(reference["fold"])
                np.testing.assert_allclose(sub["label"], reference["label"], equal_nan=True)
                y, scores, folds_ = (
                    sub["label"].to_numpy(),
                    sub["score"].to_numpy(),
                    sub["fold"].to_numpy(),
                )
                observed = np.isfinite(y)
                for rate in rates:
                    calls = calls_at_rate(scores, folds_, rate)
                    curves.append(
                        dict(
                            split_id=split,
                            seed=seed,
                            model=kind,
                            full_pool_rate=rate,
                            mcc=float(matthews_corrcoef(y[observed], calls[observed])),
                            projected_prior_03=projected(y[observed], calls[observed], 0.03),
                            projected_prior_07=projected(y[observed], calls[observed], 0.07),
                            projected_prior_12=projected(y[observed], calls[observed], 0.12),
                            usage="retrospective sensitivity; never threshold fitting",
                        )
                    )
    pl.DataFrame(curves).write_csv(out / "whole-ensemble-ranking-diagnostics.csv")


def paired_intervals(frame, out, cfg):
    """Resample molecular groups, keeping repeated seeds paired; conditional on fits."""
    comparisons = []
    baseline_model = cfg.get("baseline_model", "fresh_d6_network")
    for split in cfg["split_strategies"]:
        source = frame.filter((pl.col("split_id") == split) & pl.col("label").is_finite())
        targets = [(baseline_model, p) for p in ["single_refit", "majority_vote", "pooled_oof_cut"]]
        if cfg.get("rank_transfer_ablation"):
            targets.append((baseline_model, "rank_margin_mean"))
        targets += [(m, "centered_fold_mean") for m in cfg["models"] if m != baseline_model]
        baseline = source.filter(
            (pl.col("model") == baseline_model) & (pl.col("policy") == "centered_fold_mean")
        )
        rng = np.random.default_rng(20261002)
        for model, policy in targets:
            deltas = []
            joined = baseline.select(
                "compound_id", "seed", "group", "label", pl.col("call").alias("new_call")
            ).join(
                source.filter((pl.col("model") == model) & (pl.col("policy") == policy)).select(
                    "compound_id", "seed", pl.col("call").alias("other_call")
                ),
                on=["compound_id", "seed"],
                how="inner",
            )
            assert joined.height == baseline.height
            unique_groups = joined["group"].unique().sort().to_list()
            group_rows = {g: np.flatnonzero(joined["group"].to_numpy() == g) for g in unique_groups}
            y = joined["label"].to_numpy()
            a, b = joined["new_call"].to_numpy(), joined["other_call"].to_numpy()
            seeds = joined["seed"].to_numpy()
            for _ in range(cfg["bootstrap_samples"]):
                sample = rng.integers(0, len(unique_groups), len(unique_groups))
                idx = np.concatenate([group_rows[unique_groups[i]] for i in sample])
                delta = [
                    matthews_corrcoef(y[idx[seeds[idx] == s]], a[idx[seeds[idx] == s]])
                    - matthews_corrcoef(y[idx[seeds[idx] == s]], b[idx[seeds[idx] == s]])
                    for s in cfg["seeds"]
                ]
                deltas.append(np.mean(delta))
            actual = np.mean(
                [
                    matthews_corrcoef(y[seeds == s], a[seeds == s])
                    - matthews_corrcoef(y[seeds == s], b[seeds == s])
                    for s in cfg["seeds"]
                ]
            )
            comparisons.append(
                dict(
                    split_id=split,
                    baseline="network_centered_fold_mean",
                    comparator_model=model,
                    comparator_policy=policy,
                    delta=float(actual),
                    ci_low=float(np.quantile(deltas, 0.025)),
                    ci_high=float(np.quantile(deltas, 0.975)),
                    limits="conditional group bootstrap; seeds paired, "
                    "training/selection uncertainty omitted",
                )
            )
    pl.DataFrame(comparisons).write_csv(out / "paired-group-intervals.csv")
