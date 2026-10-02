"""Conditional saved-OOF TDI diagnostics; never supports deployment decisions."""

import itertools
import json
import os
import runpy
from functools import partial
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score


def rank_fold(score, folds):
    result = np.empty(len(score))
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = rankdata(score[mask]) / (mask.sum() + 1)
    return result


def calls_at(score, folds, rate):
    result = np.zeros(len(score), dtype=bool)
    for fold in np.unique(folds):
        idx = np.flatnonzero(folds == fold)
        k = max(1, int(round(len(idx) * rate)))
        result[idx[np.argsort(-score[idx], kind="stable")[:k]]] = True
    return result


def mcc(y, calls, prior=None):
    positive = y == 1
    if not positive.any() or positive.all():
        return float("nan")
    prior = positive.mean() if prior is None else prior
    tpr, fpr = calls[positive].mean(), calls[~positive].mean()
    tp, fn = prior * tpr, prior * (1 - tpr)
    fp, tn = (1 - prior) * fpr, (1 - prior) * (1 - fpr)
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denom) if denom else 0.0


def quality(y, score, folds, mask, objective, rate, prior):
    if objective == "AP":
        return float(average_precision_score(y[mask], score[mask]))
    return mcc(y[mask], calls_at(score, folds, rate)[mask], prior)


def compose(ns, target, additions=(), replacements=None):
    replacements = replacements or {}
    folds, y = ns[f"{target}_folds"], ns[f"{target}_y"]
    result = np.empty(len(y))
    for choice in ns[f"{target}_choices"]:
        members = choice["members"]
        if target == "a4":
            values = [replacements.get(name, ns["a4_atoms"][name]) for name in members]
            score = np.mean(values + list(additions), axis=0)
        else:
            values = []
            for name in members:
                value = ns["d6_atoms"][name]
                if name == "pfn::external_mixed_pca256":
                    value = (value + ns["d6_atoms"]["pfn::external_mixed_raw600"]) / 2
                values.append(replacements.get(name, value))
            core = rank_fold(np.mean(values + list(additions), axis=0), folds)
            tails = ns["d6_tail"]
            live = (2 * core + tails["logistic"]) / 3
            recipes = {
                "core": core,
                "core_logistic": live,
                "tail_rbf": (2 * live + tails["rbf"]) / 3,
                "tail_complement": (2 * live + tails["complement"]) / 3,
            }
            recipes["tail_average"] = (recipes["tail_rbf"] + recipes["tail_complement"]) / 2
            score = recipes[choice["recipe"]]
        held = folds == int(choice["held_fold"])
        result[held] = score[held]
    return result


def greedy_block(z, names, y, folds, mask, count, objective, rate, prior):
    chosen, total = [], np.zeros(len(y))
    for step in range(count):
        options = [i for i in range(len(names)) if i not in chosen]
        index = max(
            options,
            key=lambda i: (
                quality(y, (total + z[:, i]) / (step + 1), folds, mask, objective, rate, prior),
                -i,
            ),
        )
        chosen.append(index)
        total += z[:, index]
    return total / count, chosen


def greedy_replace(ns, target, z, names, y, folds, mask, cap, objective, rate, prior):
    members = list(ns["a4_atoms"]) if target == "a4" else ["pfn::external_mixed_pca256"]
    chosen, used, mapping = [], set(), {}
    score = compose(ns, target)
    for _ in range(min(cap, len(members))):
        incumbent = quality(y, score, folds, mask, objective, rate, prior)
        best = None
        for member in members:
            if member in mapping:
                continue
            for i in range(len(names)):
                if i in used:
                    continue
                proposal = compose(ns, target, replacements={**mapping, member: z[:, i]})
                q = quality(y, proposal, folds, mask, objective, rate, prior)
                if q > incumbent + 1e-12 and (best is None or q > best[0] + 1e-12):
                    best = (q, member, i, proposal)
        if best is None:
            break
        _, member, i, score = best
        mapping[member] = z[:, i]
        used.add(i)
        chosen.append({"member": member, "trunk": names[i]})
    return score, chosen


def bootstrap(y, selected, baseline, selected_calls, baseline_calls, groups, prior, draws, rng):
    indices = [np.flatnonzero(groups == group) for group in np.unique(groups)]
    deltas = []
    for _ in range(draws):
        sampled = rng.integers(0, len(indices), len(indices))
        rows = np.concatenate([indices[i] for i in sampled])
        if np.unique(y[rows]).size != 2:
            continue
        deltas.append(
            [
                average_precision_score(y[rows], selected[rows])
                - average_precision_score(y[rows], baseline[rows]),
                mcc(y[rows], selected_calls[rows], prior)
                - mcc(y[rows], baseline_calls[rows], prior),
            ]
        )
    return np.quantile(np.asarray(deltas), [0.025, 0.975], axis=0)


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    # Reuse the exact audited reconstruction code, with scratch outputs in this run.
    scratch = Path(os.environ["TMPDIR"]) / "reconstruction"
    scratch.mkdir()
    os.environ["CYP_OUTPUT_DIR"] = str(scratch)
    ns = runpy.run_path(
        str(Path(__file__).with_name("evaluate_classification_depth_whole_members.py"))
    )
    os.environ["CYP_OUTPUT_DIR"] = str(out)
    labels = pl.read_parquet(inp / "measurements.parquet")
    assignments = pl.read_parquet(inp / "assignments.parquet").filter(pl.col("strategy") == "umap")
    manifest = pl.read_parquet(inp / "tensor-manifest.parquet")
    old = pl.read_parquet(inp / "placement-oof.parquet")
    metric_rows, curve_rows, selection_rows, corr_rows, pair_rows = [], [], [], [], []
    call_rows, oof_rows, interval_rows, audit_rows = [], [], [], []
    rng = np.random.default_rng(cfg["seed"])
    for target in ["d6", "a4"]:
        setting = cfg[target]
        rate, prior = setting["call_rate"], setting["target_prior"]
        saved, endpoint = ns[f"{target}_saved"], ns[f"{target}_endpoint"]
        folds, y = ns[f"{target}_folds"], ns[f"{target}_y"]
        ids = saved.select("compound_id", "fold")
        ref = labels.filter(pl.col("endpoint") == endpoint).select("compound_id", "value")
        checked = (
            saved.select("compound_id", "fold", "y_true")
            .join(ref, on="compound_id", validate="1:1")
            .join(
                assignments.select("compound_id", pl.col("fold").alias("expected_fold")),
                on="compound_id",
                validate="1:1",
            )
        )
        assert checked.height == len(y)
        assert (
            checked.filter(
                (pl.col("y_true") != pl.col("value")) | (pl.col("fold") != pl.col("expected_fold"))
            ).height
            == 0
        )
        groups = ids.join(
            assignments.select("compound_id", "group"),
            on="compound_id",
            validate="1:1",
            maintain_order="left",
        )["group"].to_numpy()
        current = compose(ns, target)
        assert np.max(np.abs(current - ns[f"{target}_current"])) < 1e-12
        source = old.filter((pl.col("target") == target) & (pl.col("scenario") == "current"))
        assert (
            np.max(
                np.abs(
                    ids.join(
                        source.select("compound_id", "fold", "prediction"),
                        on=["compound_id", "fold"],
                        validate="1:1",
                        maintain_order="left",
                    )["prediction"].to_numpy()
                    - current
                )
            )
            < 1e-12
        )
        raw = pl.read_parquet(inp / cfg["trunk_predictions"][target]).filter(
            pl.col("role") == "train"
        )
        all_names = sorted(raw["model"].unique())
        names = [name for name in all_names if not name.startswith("molecular_only")]
        scores = {}
        raw_scores = {}
        for name in all_names:
            p = ids.join(
                raw.filter(pl.col("model") == name).select("compound_id", "fold", "prediction"),
                on=["compound_id", "fold"],
                validate="1:1",
                maintain_order="left",
            )["prediction"]
            assert p.is_finite().all()
            raw_scores[name] = p.to_numpy()
            scores[name] = rank_fold(p.to_numpy(), folds)
        z = np.column_stack([scores[name] for name in names])
        atoms = dict(ns[f"{target}_atoms"])
        if target == "d6":
            atoms = {name: rank_fold(value, folds) for name, value in atoms.items()}
            atoms.update({f"tail::{name}": value for name, value in ns["d6_tail"].items()})
        references = {"incumbent": rank_fold(current, folds), **atoms}
        base_rank = rank_fold(current, folds)
        for name in all_names:
            for reference, r in references.items():
                for fold in [-1, *sorted(np.unique(folds))]:
                    mask = np.ones(len(y), dtype=bool) if fold == -1 else folds == fold
                    corr_rows.append(
                        {
                            "target": target,
                            "model": name,
                            "reference": reference,
                            "fold": int(fold),
                            "prediction_spearman": float(
                                spearmanr(scores[name][mask], r[mask]).statistic
                            ),
                            "rank_residual_pearson": float(
                                np.corrcoef((scores[name] - y)[mask], (r - y)[mask])[0, 1]
                            ),
                            "rank_squared_loss_pearson": float(
                                np.corrcoef(((scores[name] - y) ** 2)[mask], ((r - y) ** 2)[mask])[
                                    0, 1
                                ]
                            ),
                        }
                    )
        for i, j in itertools.combinations(range(len(names)), 2):
            pair_rows.append(
                {
                    "target": target,
                    "left": names[i],
                    "right": names[j],
                    "prediction_spearman": float(spearmanr(z[:, i], z[:, j]).statistic),
                    "rank_residual_pearson": float(np.corrcoef(z[:, i] - y, z[:, j] - y)[0, 1]),
                }
            )
        scenarios = {
            name: value
            for name, value in ns[f"{target}_scenarios"].items()
            if name != "outer_selected_AP_placement"
        }
        # Ancillary molecular-only placements show how much a new readout alone contributes.
        for name in all_names:
            if name.startswith("molecular_only"):
                scenarios[f"molecular_add::{name}"] = compose(ns, target, [scores[name]])
                members = list(ns["a4_atoms"]) if target == "a4" else ["pfn::external_mixed_pca256"]
                for member in members:
                    scenarios[f"molecular_replace_{member}::{name}"] = compose(
                        ns, target, replacements={member: scores[name]}
                    )
        eligible = {
            name: score
            for name, score in scenarios.items()
            if not name.startswith(("standalone::", "molecular_"))
        }
        # Same fixed bank plus original members only: weight changes without any trunk.
        no_trunk = {"current": current}
        for name, score in atoms.items():
            no_trunk[f"original_add::{name}"] = compose(ns, target, [rank_fold(score, folds)])
        scenarios.update(no_trunk)
        baseline_calls = calls_at(current, folds, rate)
        policies = {}
        for objective in ["AP", "MCC"]:
            selected = np.empty(len(y))
            selected_rates = np.full(len(y), rate)
            joint = np.empty(len(y))
            joint_rates = np.empty(len(y))
            control = np.empty(len(y))
            block_predictions = np.empty(len(y))
            replacements = np.empty(len(y))
            for held in sorted(np.unique(folds)):
                train, valid = folds != held, folds == held
                q = partial(
                    quality, y, folds=folds, mask=train, objective=objective, rate=rate, prior=prior
                )
                chosen = max(
                    eligible, key=lambda name: (q(eligible[name]), name == "current", name)
                )
                selected[valid] = eligible[chosen][valid]
                selection_rows.append(
                    {
                        "target": target,
                        "objective": objective,
                        "family": "single",
                        "fold": int(held),
                        "recipe": chosen,
                        "rate": rate,
                    }
                )
                control_name = max(
                    no_trunk, key=lambda name: (q(no_trunk[name]), name == "current", name)
                )
                control[valid] = no_trunk[control_name][valid]
                selection_rows.append(
                    {
                        "target": target,
                        "objective": objective,
                        "family": "no_trunk",
                        "fold": int(held),
                        "recipe": control_name,
                        "rate": rate,
                    }
                )
                if objective == "MCC":
                    option = max(
                        ((name, r) for name in eligible for r in setting["rates"]),
                        key=lambda item: (
                            quality(y, eligible[item[0]], folds, train, "MCC", item[1], prior),
                            item[0] == "current",
                            -abs(item[1] - rate),
                            item[0],
                        ),
                    )
                    joint[valid] = eligible[option[0]][valid]
                    joint_rates[valid] = option[1]
                    selection_rows.append(
                        {
                            "target": target,
                            "objective": objective,
                            "family": "single_joint_rate",
                            "fold": int(held),
                            "recipe": option[0],
                            "rate": option[1],
                        }
                    )
                # Inner saved-fold validation chooses block size and placement; control is included.
                inner_options = {("current", 0, 0.0): current.copy()}
                for count in cfg["block_sizes"]:
                    for kind, weight in [
                        ("block_add", 0.0),
                        ("individual_add", 0.0),
                        *[("blend", w) for w in cfg["blend_weights"]],
                    ]:
                        inner_options[(kind, count, weight)] = np.empty(len(y))
                for inner in sorted(np.unique(folds[train])):
                    fit_mask, check_mask = train & (folds != inner), train & (folds == inner)
                    for count in cfg["block_sizes"]:
                        block, indices = greedy_block(
                            z, names, y, folds, fit_mask, count, objective, rate, prior
                        )
                        values = {
                            ("block_add", count, 0.0): compose(ns, target, [block]),
                            ("individual_add", count, 0.0): compose(
                                ns, target, [z[:, i] for i in indices]
                            ),
                        }
                        for w in cfg["blend_weights"]:
                            values[("blend", count, w)] = (1 - w) * base_rank + w * block
                        for key, score in values.items():
                            inner_options[key][check_mask] = score[check_mask]
                best = max(
                    inner_options,
                    key=lambda key: (q(inner_options[key]), key[0] == "current", -key[1], -key[2]),
                )
                kind, count, weight = best
                if kind == "current":
                    final, chosen_names = current, []
                else:
                    block, indices = greedy_block(
                        z, names, y, folds, train, count, objective, rate, prior
                    )
                    chosen_names = [names[i] for i in indices]
                    if kind == "block_add":
                        final = compose(ns, target, [block])
                    elif kind == "individual_add":
                        final = compose(ns, target, [z[:, i] for i in indices])
                    else:
                        final = (1 - weight) * base_rank + weight * block
                block_predictions[valid] = final[valid]
                selection_rows.append(
                    {
                        "target": target,
                        "objective": objective,
                        "family": "multi_block",
                        "fold": int(held),
                        "recipe": json.dumps(
                            {"kind": kind, "count": count, "weight": weight, "models": chosen_names}
                        ),
                        "rate": rate,
                    }
                )
                replacement, choices = greedy_replace(
                    ns,
                    target,
                    z,
                    names,
                    y,
                    folds,
                    train,
                    cfg["replacement_cap"],
                    objective,
                    rate,
                    prior,
                )
                replacements[valid] = replacement[valid]
                selection_rows.append(
                    {
                        "target": target,
                        "objective": objective,
                        "family": "multi_replace",
                        "fold": int(held),
                        "recipe": json.dumps(choices),
                        "rate": rate,
                    }
                )
                print(
                    f"{target} {objective} held={held} single={chosen} "
                    f"block={best} replacements={len(choices)}",
                    flush=True,
                )
            policies[f"single_selected_{objective}"] = (selected, selected_rates)
            policies[f"no_trunk_selected_{objective}"] = (control, selected_rates)
            policies[f"multi_block_selected_{objective}"] = (block_predictions, selected_rates)
            policies[f"multi_replace_selected_{objective}"] = (replacements, selected_rates)
            if objective == "MCC":
                policies["single_joint_rate_MCC"] = (joint, joint_rates)
        scenarios.update({name: values[0] for name, values in policies.items()})
        # Rate search and all fixed curves remain diagnostic; no full-pool inference.
        for name, score in scenarios.items():
            rates = policies[name][1] if name in policies else np.full(len(y), rate)
            call = np.zeros(len(y), dtype=bool)
            for held in np.unique(folds):
                mask = folds == held
                call[mask] = calls_at(score, folds, rates[mask][0])[mask]
            metric_rows.append(
                {
                    "target": target,
                    "scenario": name,
                    "n": len(y),
                    "AP": float(average_precision_score(y, score)),
                    "AUROC": float(roc_auc_score(y, score)),
                    "MCC_observed": mcc(y, call),
                    "MCC_target_prior": mcc(y, call, prior),
                    "added": int((call & ~baseline_calls).sum()),
                    "removed": int((~call & baseline_calls).sum()),
                    "positive_added": int((call & ~baseline_calls & (y == 1)).sum()),
                    "positive_removed": int((~call & baseline_calls & (y == 1)).sum()),
                }
            )
            for r in setting["rates"]:
                c = calls_at(score, folds, r)
                for p in setting["priors"]:
                    curve_rows.append(
                        {
                            "target": target,
                            "scenario": name,
                            "rate": r,
                            "prior": p,
                            "MCC_target_prior": mcc(y, c, p),
                            "MCC_observed": mcc(y, c),
                            "calls": int(c.sum()),
                        }
                    )
            if name == "current" or name in policies:
                oof_rows.append(
                    ids.with_columns(
                        pl.lit(endpoint).alias("endpoint"),
                        pl.lit("20260909T135130268493Z-edba6fc0:umap").alias("split_id"),
                        pl.Series("group", groups),
                        pl.Series("y_true", y),
                        pl.lit(name).alias("model"),
                        pl.Series("prediction", score),
                        pl.Series("call", call),
                        pl.Series("diagnostic_rate", rates),
                    )
                )
                for held in np.unique(folds):
                    mask = folds == held
                    metric_rows.append(
                        {
                            "target": target,
                            "scenario": f"{name}::fold{held}",
                            "n": int(mask.sum()),
                            "AP": float(average_precision_score(y[mask], score[mask])),
                            "AUROC": float(roc_auc_score(y[mask], score[mask])),
                            "MCC_observed": mcc(y[mask], call[mask]),
                            "MCC_target_prior": mcc(y[mask], call[mask], prior),
                        }
                    )
                changed = call != baseline_calls
                call_rows.append(
                    ids.filter(pl.Series(changed)).with_columns(
                        pl.lit(endpoint).alias("endpoint"),
                        pl.lit("20260909T135130268493Z-edba6fc0:umap").alias("split_id"),
                        pl.lit(name).alias("model"),
                        pl.Series("y_true", y[changed]),
                        pl.Series("baseline_call", baseline_calls[changed]),
                        pl.Series("call", call[changed]),
                        pl.Series("prediction", score[changed]),
                    )
                )
                if name != "current":
                    for comparison in ["current", f"no_trunk_selected_{name.rsplit('_', 1)[-1]}"]:
                        if comparison not in scenarios or name == comparison:
                            continue
                        reference = scenarios[comparison]
                        reference_call = calls_at(reference, folds, rate)
                        bounds = bootstrap(
                            y,
                            score,
                            reference,
                            call,
                            reference_call,
                            groups,
                            prior,
                            cfg["bootstrap_draws"],
                            rng,
                        )
                        interval_rows.append(
                            {
                                "target": target,
                                "scenario": name,
                                "reference": comparison,
                                "delta_AP": float(
                                    average_precision_score(y, score)
                                    - average_precision_score(y, reference)
                                ),
                                "AP_lower": float(bounds[0, 0]),
                                "AP_upper": float(bounds[1, 0]),
                                "delta_MCC": mcc(y, call, prior) - mcc(y, reference_call, prior),
                                "MCC_lower": float(bounds[0, 1]),
                                "MCC_upper": float(bounds[1, 1]),
                            }
                        )
        for name, p in raw_scores.items():
            c = calls_at(p, folds, rate)
            metric_rows.append(
                {
                    "target": target,
                    "scenario": f"raw_probability::{name}",
                    "n": len(y),
                    "AP": float(average_precision_score(y, p)),
                    "AUROC": float(roc_auc_score(y, p)),
                    "Brier": float(np.mean((p - y) ** 2)),
                    "MCC_observed": mcc(y, c),
                    "MCC_target_prior": mcc(y, c, prior),
                }
            )
        target_name = "cyp2d6" if target == "d6" else "cyp3a4"
        coverage = manifest.filter(
            (pl.col("target") == target_name) & (pl.col("role") == "train")
        ).height
        audit_rows.append(
            {
                "target": target,
                "labelled_rows": len(y),
                "total_pool_rows": assignments.height,
                "trunk_covered_rows": coverage,
                "missing_rows": assignments.height - coverage,
                "reconstruction_max_diff": float(np.max(np.abs(current - ns[f"{target}_current"]))),
                "deployment_decision_supported": False,
            }
        )
        print(f"{target} diagnostics complete", flush=True)
    for filename, rows in [
        ("metrics.parquet", metric_rows),
        ("rate-curves.parquet", curve_rows),
        ("selections.parquet", selection_rows),
        ("incumbent-correlations.parquet", corr_rows),
        ("trunk-correlations.parquet", pair_rows),
        ("intervals.parquet", interval_rows),
    ]:
        pl.DataFrame(rows).write_parquet(out / filename)
    pl.concat(oof_rows).write_parquet(out / "selected-oof.parquet")
    pl.concat(call_rows).write_parquet(out / "changed-labelled-calls.parquet")
    summary = {
        "config": cfg,
        "coverage": audit_rows,
        "new_model_fits": 0,
        "new_trunk_extractions": 0,
        "error_definition": (
            "Residual and squared loss of within-fold ranks versus binary labels; "
            "not calibrated-probability residuals. Raw Brier scores descriptive only."
        ),
        "uncertainty": (
            "Paired UMAP-group bootstrap conditional on saved scores and selected calls; "
            "no reranking or selector refit in bootstrap; no multiplicity correction."
        ),
        "limitations": (
            "Labelled-only operating pool; conditional saved-OOF selection, "
            "not fully nested base-model fitting. Multiple views reuse the same frozen trunk. "
            "No submission or deployment rate."
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
