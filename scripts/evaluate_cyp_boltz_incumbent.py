"""Whole-member trunk placement on exactly reproduced incumbent regression OOF."""

import itertools
import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def error(low, high, prediction):
    return np.maximum(low - prediction, 0) + np.maximum(prediction - high, 0)


def subset_weights(x, low, high):
    winner = (float("inf"), x.shape[1] + 1, ())
    for size in range(1, x.shape[1] + 1):
        for indices in itertools.combinations(range(x.shape[1]), size):
            result = (float(error(low, high, x[:, indices].mean(axis=1)).mean()), size, indices)
            if result < winner:
                winner = result
    weights = np.zeros(x.shape[1])
    weights[list(winner[2])] = 1 / len(winner[2])
    return weights


def caruana(x, low, high, steps):
    counts = np.zeros(x.shape[1])
    total = np.zeros(len(low))
    for step in range(steps):
        # Score columns separately so adding an unused column cannot alter the
        # floating-point reduction order or an incumbent tie-breaking decision.
        choices = [
            (float(error(low, high, (total + x[:, col]) / (step + 1)).sum()), col)
            for col in range(x.shape[1])
        ]
        _, chosen = min(choices)
        counts[chosen] += 1
        total += x[:, chosen]
    return counts / steps


def selected_weights(x, low, high, folds, method):
    if method == "equal_subset":
        return subset_weights(x, low, high)
    if method == "equal_all":
        return np.full(x.shape[1], 1 / x.shape[1])
    choices = []
    for steps in [5, 10, 20, 50]:
        prediction = np.empty(len(low))
        for held in sorted(set(folds)):
            train = folds != held
            weights = caruana(x[train], low[train], high[train], steps)
            prediction[~train] = x[~train] @ weights
        choices.append((float(error(low, high, prediction).mean()), steps))
    steps = min(choices)[1]
    return caruana(x, low, high, steps)


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    keys = ["compound_id", "endpoint", "fold"]
    rows, predictions, fits, intervals, checks = [], [], [], [], []
    for target, setup in cfg["targets"].items():
        bank = pl.read_parquet(inp / setup["bank_oof"])
        for column, value in setup["filter"].items():
            bank = bank.filter(pl.col(column) == value)
        bank = bank.filter(pl.col("endpoint") == setup["endpoint"]).sort(keys)
        if bank.height != setup["expected_rows"] or bank.select(keys).is_duplicated().any():
            raise ValueError(f"Wrong incumbent cohort: {target}")
        table = bank.select(*keys, "y_true", "conf_low", "conf_high")
        folds = table["fold"].to_numpy()
        if target == "cyp2d6":
            base = pl.read_parquet(inp / "d6-base.parquet").filter(
                (pl.col("endpoint") == setup["endpoint"]) & (pl.col("scenario") == "current_two")
            )
            table = table.join(
                base.select(*keys, pl.col("prediction").alias("base")), on=keys, validate="1:1"
            ).sort(keys)
            original_names = ["submitted_two_member_core", "moment_matched_robust_branch"]
            base_prediction = table["base"].to_numpy()
            x = np.column_stack(
                [base_prediction, 2 * bank["prediction"].to_numpy() - base_prediction]
            )
            weight_matrix = np.full(x.shape, 0.5)
        else:
            components = json.loads((inp / setup["bank_summary"]).read_text())["config"][
                "components"
            ]
            original_names = [item["name"] for item in components]
            for item in components:
                atom = pl.read_parquet(inp / f"{target}-{item['name']}.parquet").filter(
                    pl.col("endpoint") == setup["endpoint"]
                )
                table = table.join(
                    atom.select(*keys, pl.col("prediction").alias(item["name"])),
                    on=keys,
                    validate="1:1",
                ).sort(keys)
            x = table.select(original_names).to_numpy()
            weight_matrix = np.zeros(x.shape)
            if setup["method"] == "equal_all":
                weight_matrix[:] = 1 / len(original_names)
            else:
                saved_fits = pl.read_csv(inp / setup["bank_fits"]).filter(
                    pl.col("method") == setup["method"]
                )
                for row in saved_fits.iter_rows(named=True):
                    weight_matrix[folds == row["outer_fold"]] = json.loads(row["weights"])
        current = (x * weight_matrix).sum(axis=1)
        difference = float(np.max(np.abs(current - bank["prediction"].to_numpy())))
        if not np.isfinite(x).all() or difference > 1e-9:
            raise ValueError(f"Incumbent reconstruction failed: {target} {difference}")
        cohort = pl.read_parquet(inp / setup["cohort"]).select(
            "compound_id", "fold", "group", "value", "conf_low", "conf_high"
        )
        aligned = table.select(*keys, "y_true", "conf_low", "conf_high").join(
            cohort.rename(
                {"value": "cohort_value", "conf_low": "cohort_low", "conf_high": "cohort_high"}
            ),
            on=["compound_id", "fold"],
            validate="1:1",
            maintain_order="left",
        )
        if aligned.height != table.height or aligned["group"].null_count():
            raise ValueError("Missing matched cohort")
        for left, right in [
            ("y_true", "cohort_value"),
            ("conf_low", "cohort_low"),
            ("conf_high", "cohort_high"),
        ]:
            if not aligned[left].equals(aligned[right]):
                raise ValueError(f"Changed truth: {target} {left}")
        trunk = pl.read_parquet(inp / setup["predictions"]).filter(pl.col("role") == "train")
        scenarios = {"current": current}
        operations = {}
        for name in sorted(trunk["model"].unique()):
            part = trunk.filter(pl.col("model") == name)
            joined = table.select(*keys).join(
                part.select(*keys, "prediction"), on=keys, validate="1:1", maintain_order="left"
            )
            if joined.height != table.height or not joined["prediction"].is_finite().all():
                raise ValueError(f"Trunk alignment failed: {target} {name}")
            p = joined["prediction"].to_numpy()
            scenarios[f"standalone::{name}"] = p
            if name.startswith("molecular_only"):
                continue
            active = (weight_matrix > 0).sum(axis=1)
            if setup["method"] != "caruana":
                operations[f"add::{name}"] = (active * current + p) / (active + 1)
            for index, member in enumerate(original_names):
                operations[f"replace_{member}::{name}"] = current + weight_matrix[:, index] * (
                    p - x[:, index]
                )
        scenarios.update(operations)
        low, high = table["conf_low"].to_numpy(), table["conf_high"].to_numpy()
        y = table["y_true"].to_numpy()
        denominator = float(error(low, high, y.mean()).sum())
        conditional = np.empty(table.height)
        reselected = np.empty(table.height)
        reselected_control = np.empty(table.height)
        standalone_names = [
            name
            for name in scenarios
            if name.startswith("standalone::") and "molecular_only" not in name
        ]
        for held in range(5):
            train = folds != held
            choices = [
                (float(error(low[train], high[train], p[train]).mean()), name)
                for name, p in operations.items()
            ]
            choices.append(
                (float(error(low[train], high[train], current[train]).mean()), "current")
            )
            _, name = min(choices)
            conditional[~train] = scenarios[name][~train]
            fits.append(
                {
                    "target": target,
                    "held_fold": held,
                    "selection": "placement",
                    "chosen": name,
                    "weights": "",
                }
            )
            _, best = min(
                (float(error(low[train], high[train], scenarios[name][train]).mean()), name)
                for name in standalone_names
            )
            extension = np.column_stack([x, scenarios[best]])
            method = setup["method"] if target != "cyp2d6" else "caruana"
            original_weights = selected_weights(
                x[train], low[train], high[train], folds[train], method
            )
            reselected_control[~train] = x[~train] @ original_weights
            fits.append(
                {
                    "target": target,
                    "held_fold": held,
                    "selection": "original_bank_reselection",
                    "chosen": "no_trunk",
                    "weights": json.dumps(original_weights.tolist()),
                }
            )
            weights = selected_weights(
                extension[train], low[train], high[train], folds[train], method
            )
            reselected[~train] = extension[~train] @ weights
            fits.append(
                {
                    "target": target,
                    "held_fold": held,
                    "selection": "bank_reselection",
                    "chosen": best,
                    "weights": json.dumps(weights.tolist()),
                }
            )
        scenarios["outer_selected_placement"] = conditional
        scenarios["original_bank_reselection"] = reselected_control
        scenarios["outer_selected_trunk_bank_reselection"] = reselected
        for name, p in scenarios.items():
            for held in [-1, 0, 1, 2, 3, 4]:
                mask = np.ones(len(y), dtype=bool) if held == -1 else folds == held
                d = float(error(low[mask], high[mask], y[mask].mean()).sum())
                rows.append(
                    {
                        "target": target,
                        "scenario": name,
                        "fold": held,
                        "n": int(mask.sum()),
                        "ST_RAE": float(error(low[mask], high[mask], p[mask]).sum() / d),
                        "MAE": float(np.abs(y[mask] - p[mask]).mean()),
                    }
                )
            predictions.append(
                table.select(*keys).with_columns(
                    pl.lit(target).alias("target"),
                    pl.lit(name).alias("scenario"),
                    pl.Series("prediction", p),
                )
            )
        groups = aligned["group"].to_numpy()
        unique = sorted(set(groups))
        group_index = {group: i for i, group in enumerate(unique)}
        index = np.array([group_index[g] for g in groups])
        rng = np.random.default_rng(20261001)
        draws = rng.integers(0, len(unique), size=(2000, len(unique)))
        counts = np.bincount(index, minlength=len(unique))
        for name, reference in [
            ("outer_selected_placement", "current"),
            ("outer_selected_trunk_bank_reselection", "current"),
            ("bank_extension_vs_reselected_control", "original_bank_reselection"),
        ]:
            candidate = (
                reselected if name == "bank_extension_vs_reselected_control" else scenarios[name]
            )
            delta = error(low, high, candidate) - error(low, high, scenarios[reference])
            totals = np.bincount(index, weights=delta, minlength=len(unique))
            bootstrap = (
                totals[draws].sum(axis=1) / counts[draws].sum(axis=1) / (denominator / len(y))
            )
            intervals.append(
                {
                    "target": target,
                    "scenario": name,
                    "reference": reference,
                    "delta_ST_RAE": float(delta.sum() / denominator),
                    "lower_95": float(np.quantile(bootstrap, 0.025)),
                    "upper_95": float(np.quantile(bootstrap, 0.975)),
                    "groups": len(unique),
                    "conditional_on_saved_oof": True,
                }
            )
        checks.append(
            {"target": target, "n": table.height, "incumbent_max_absolute_difference": difference}
        )
        pl.DataFrame(rows).write_parquet(out / "metrics.parquet")
        print(json.dumps(checks[-1]), flush=True)
    pl.concat(predictions).write_parquet(out / "oof.parquet")
    pl.DataFrame(fits).write_parquet(out / "selections.parquet")
    pl.DataFrame(intervals).write_parquet(out / "paired-uncertainty.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "checks": checks,
                "scope": (
                    "Conditional saved-OOF selection, not fully nested base-model refits; "
                    "no candidate or submission"
                ),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
