"""Conditional saved-OOF trunk diversity, blocks, and multiple-member placement."""

import itertools
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from evaluate_cyp_boltz_incumbent import error, selected_weights


def greedy_block(z, low, high, mask, count):
    chosen, total = [], np.zeros(len(low))
    for step in range(count):
        scores = error(
            low[mask, None], high[mask, None], (total[mask, None] + z[mask]) / (step + 1)
        ).mean(axis=0)
        scores[chosen] = np.inf
        index = int(np.argmin(scores))
        chosen.append(index)
        total += z[:, index]
    return total / count, chosen


def replacements(z, x, weights, current, low, high, mask, count):
    p, choices, used_members, used_trunks = current.copy(), [], set(), set()
    for _ in range(count):
        best = (float(error(low[mask], high[mask], p[mask]).mean()), -1, -1)
        for j in range(x.shape[1]):
            if j in used_members or not np.any(weights[mask, j] > 0):
                continue
            trial = p[:, None] + weights[:, j, None] * (z - x[:, j, None])
            scores = error(low[mask, None], high[mask, None], trial[mask]).mean(axis=0)
            if used_trunks:
                scores[list(used_trunks)] = np.inf
            k = int(np.argmin(scores))
            best = min(best, (float(scores[k]), j, k))
        if best[1] < 0:
            break
        _, j, k = best
        p += weights[:, j] * (z[:, k] - x[:, j])
        used_members.add(j)
        used_trunks.add(k)
        choices.append((j, k))
    return p, choices


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    keys = ["compound_id", "endpoint", "fold"]
    metrics, fits, correlations, intervals, oof = [], [], [], [], []
    for target, setup in cfg["targets"].items():
        bank = pl.read_parquet(inp / setup["bank_oof"])
        for col, value in setup["filter"].items():
            bank = bank.filter(pl.col(col) == value)
        bank = bank.filter(pl.col("endpoint") == setup["endpoint"]).sort(keys)
        if bank.height != setup["expected_rows"] or bank.select(keys).is_duplicated().any():
            raise ValueError(f"Wrong cohort {target}")
        table = bank.select(*keys, "y_true", "conf_low", "conf_high")
        folds = table["fold"].to_numpy()
        current = bank["prediction"].to_numpy()
        if target == "cyp2d6":
            base = pl.read_parquet(inp / "d6-base.parquet").filter(
                (pl.col("endpoint") == setup["endpoint"]) & (pl.col("scenario") == "current_two")
            )
            p = (
                table.select(keys)
                .join(
                    base.select(*keys, "prediction"), on=keys, validate="1:1", maintain_order="left"
                )["prediction"]
                .to_numpy()
            )
            atom_names = ["submitted_two_member_core", "moment_matched_robust_branch"]
            x = np.column_stack([p, 2 * current - p])
            weights = np.full(x.shape, 0.5)
        else:
            specs = json.loads((inp / setup["bank_summary"]).read_text())["config"]["components"]
            atom_names = [s["name"] for s in specs]
            for name in atom_names:
                atom = pl.read_parquet(inp / f"{target}-{name}.parquet").filter(
                    pl.col("endpoint") == setup["endpoint"]
                )
                table = table.join(
                    atom.select(*keys, pl.col("prediction").alias(name)),
                    on=keys,
                    validate="1:1",
                    maintain_order="left",
                )
            x = table.select(atom_names).to_numpy()
            weights = np.zeros(x.shape)
            if setup["method"] == "equal_all":
                weights[:] = 1 / len(atom_names)
            else:
                saved = pl.read_csv(inp / setup["bank_fits"]).filter(
                    pl.col("method") == setup["method"]
                )
                for row in saved.iter_rows(named=True):
                    weights[folds == row["outer_fold"]] = json.loads(row["weights"])
        if np.max(np.abs((x * weights).sum(axis=1) - current)) > 1e-9:
            raise ValueError("Incumbent reconstruction mismatch")
        cohort = pl.read_parquet(inp / setup["cohort"])
        aligned = table.select("compound_id", "fold").join(
            cohort, on=["compound_id", "fold"], validate="1:1", maintain_order="left"
        )
        for a, b in [("y_true", "value"), ("conf_low", "conf_low"), ("conf_high", "conf_high")]:
            if not table[a].equals(aligned[b]):
                raise ValueError("Changed labels or bounds")
        trunk = pl.read_parquet(inp / setup["predictions"]).filter(pl.col("role") == "train")
        names = sorted(n for n in trunk["model"].unique() if not n.startswith("molecular_only"))
        z = np.column_stack(
            [
                table.select(keys)
                .join(
                    trunk.filter(pl.col("model") == name).select(*keys, "prediction"),
                    on=keys,
                    validate="1:1",
                    maintain_order="left",
                )["prediction"]
                .to_numpy()
                for name in names
            ]
        )
        if not np.isfinite(z).all() or not np.isfinite(x).all() or len(names) != 72:
            raise ValueError("Incomplete model bank")
        y, low, high = (table[c].to_numpy() for c in ["y_true", "conf_low", "conf_high"])
        for held in [-1, 0, 1, 2, 3, 4]:
            mask = np.ones(len(y), dtype=bool) if held == -1 else folds == held
            pred_corr = np.corrcoef(z[mask].T)
            err_corr = np.corrcoef((y[mask, None] - z[mask]).T)
            interval_corr = np.corrcoef(
                (z[mask] - np.clip(z[mask], low[mask, None], high[mask, None])).T
            )
            for i, j in itertools.combinations(range(len(names)), 2):
                correlations.append(
                    {
                        "target": target,
                        "fold": held,
                        "model_a": names[i],
                        "model_b": names[j],
                        "prediction_pearson": float(pred_corr[i, j]),
                        "residual_pearson": float(err_corr[i, j]),
                        "interval_residual_pearson": float(interval_corr[i, j]),
                    }
                )
        scenarios = {"current": current, "all_72_trunk_equal": z.mean(axis=1)}
        scenario_names = [
            "trunk_block",
            "trunk_block_added",
            "multiple_added",
            "multiple_replaced",
            "full_bank",
            "original_bank_control",
        ]
        scenarios.update({s: np.empty(len(y)) for s in scenario_names})
        for held in range(5):
            train = folds != held
            # All recipe/count/weight decisions exclude the reported outer fold.
            # Inner folds select block size and mixing strength using newly fitted
            # greedy block selections; base OOF remains conditionally reused.
            candidates = []
            for count in cfg["block_sizes"]:
                inner = np.empty(len(y))
                for inner_held in sorted(set(folds[train])):
                    fit = train & (folds != inner_held)
                    q, _ = greedy_block(z, low, high, fit, count)
                    inner[folds == inner_held] = q[folds == inner_held]
                for alpha in cfg["blend_weights"]:
                    q = (1 - alpha) * current[train] + alpha * inner[train]
                    candidates.append(
                        (float(error(low[train], high[train], q).mean()), count, alpha)
                    )
            _, count, alpha = min(candidates)
            block, selected = greedy_block(z, low, high, train, count)
            scenarios["trunk_block"][~train] = block[~train]
            scenarios["trunk_block_added"][~train] = (1 - alpha) * current[~train] + alpha * block[
                ~train
            ]
            fits.append(
                {
                    "target": target,
                    "fold": held,
                    "method": "block_added",
                    "choices": json.dumps([names[k] for k in selected]),
                    "parameter": float(alpha),
                }
            )
            # Multiple equal-strength additions, chosen by inner held-fold error.
            active = (weights > 0).sum(axis=1)
            add_choices = []
            for count in cfg["block_sizes"]:
                inner = np.empty(len(y))
                for h in sorted(set(folds[train])):
                    q, _ = greedy_block(z, low, high, train & (folds != h), count)
                    m = folds == h
                    inner[m] = (active[m] * current[m] + count * q[m]) / (active[m] + count)
                add_choices.append(
                    (float(error(low[train], high[train], inner[train]).mean()), count)
                )
            add_choices.append((float(error(low[train], high[train], current[train]).mean()), 0))
            _, count = min(add_choices)
            q, selected = greedy_block(z, low, high, train, count) if count else (current, [])
            scenarios["multiple_added"][~train] = (
                active[~train] * current[~train] + count * q[~train]
            ) / (active[~train] + count)
            fits.append(
                {
                    "target": target,
                    "fold": held,
                    "method": "multiple_added",
                    "choices": json.dumps([names[k] for k in selected]),
                    "parameter": float(count),
                }
            )
            replace_choices = []
            # D6 retains at least one original top-level branch; other targets can
            # replace up to three distinct atomic members at their saved weights.
            caps = [0, 1] if target == "cyp2d6" else [0, 1, 2, 3]
            for cap in caps:
                inner = np.empty(len(y))
                for h in sorted(set(folds[train])):
                    q, _ = replacements(
                        z, x, weights, current, low, high, train & (folds != h), cap
                    )
                    inner[folds == h] = q[folds == h]
                replace_choices.append(
                    (float(error(low[train], high[train], inner[train]).mean()), cap)
                )
            _, cap = min(replace_choices)
            q, selected = replacements(z, x, weights, current, low, high, train, cap)
            scenarios["multiple_replaced"][~train] = q[~train]
            fits.append(
                {
                    "target": target,
                    "fold": held,
                    "method": "multiple_replaced",
                    "choices": json.dumps([(atom_names[j], names[k]) for j, k in selected]),
                    "parameter": float(cap),
                }
            )
            for label, bank_x, labels in [
                ("full_bank", np.column_stack([x, z]), atom_names + names),
                ("original_bank_control", x, atom_names),
            ]:
                w = selected_weights(
                    bank_x[train], low[train], high[train], folds[train], "caruana"
                )
                scenarios[label][~train] = bank_x[~train] @ w
                fits.append(
                    {
                        "target": target,
                        "fold": held,
                        "method": label,
                        "choices": json.dumps(
                            {n: float(v) for n, v in zip(labels, w, strict=True) if v}
                        ),
                        "parameter": float(w[len(atom_names) :].sum()),
                    }
                )
        denominator = float(error(low, high, y.mean()).sum())
        groups = aligned["group"].to_numpy()
        unique = sorted(set(groups))
        lookup = {g: i for i, g in enumerate(unique)}
        index = np.array([lookup[g] for g in groups])
        rng = np.random.default_rng(20261001)
        draws = rng.integers(len(unique), size=(2000, len(unique)))
        counts = np.bincount(index, minlength=len(unique))
        for label, p in scenarios.items():
            for held in [-1, 0, 1, 2, 3, 4]:
                m = np.ones(len(y), dtype=bool) if held == -1 else folds == held
                metrics.append(
                    {
                        "target": target,
                        "scenario": label,
                        "fold": held,
                        "n": int(m.sum()),
                        "ST_RAE": float(
                            error(low[m], high[m], p[m]).sum()
                            / error(low[m], high[m], y[m].mean()).sum()
                        ),
                        "MAE": float(np.abs(y[m] - p[m]).mean()),
                    }
                )
            oof.append(
                table.select(keys).with_columns(
                    pl.lit(target).alias("target"),
                    pl.lit(label).alias("scenario"),
                    pl.Series("prediction", p),
                )
            )
            if label in ["current", "all_72_trunk_equal", "original_bank_control"]:
                continue
            for reference in (
                ["current", "original_bank_control"] if label == "full_bank" else ["current"]
            ):
                delta = error(low, high, p) - error(low, high, scenarios[reference])
                totals = np.bincount(index, weights=delta, minlength=len(unique))
                boot = totals[draws].sum(axis=1) / counts[draws].sum(axis=1) * len(y) / denominator
                intervals.append(
                    {
                        "target": target,
                        "scenario": label,
                        "reference": reference,
                        "delta_ST_RAE": float(delta.sum() / denominator),
                        "lower_95": float(np.quantile(boot, 0.025)),
                        "upper_95": float(np.quantile(boot, 0.975)),
                    }
                )
        print(json.dumps({"target": target, "n": len(y), "trunk_models": len(names)}), flush=True)
    pl.DataFrame(correlations).write_parquet(out / "trunk-pair-correlations.parquet")
    pl.DataFrame(metrics).write_parquet(out / "metrics.parquet")
    pl.DataFrame(fits).write_parquet(out / "selections.parquet")
    pl.DataFrame(intervals).write_parquet(out / "paired-uncertainty.parquet")
    pl.concat(oof).write_parquet(out / "oof.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "scope": (
                    "Conditional saved-OOF selection; held folds excluded from "
                    "recipe/count/weight selection; not fully nested base-model refits. "
                    "D6 keeps at least one incumbent branch during replacement. "
                    "All-72 equal mean is a diagnostic. Multiple additions preserve "
                    "relative incumbent weights."
                ),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
