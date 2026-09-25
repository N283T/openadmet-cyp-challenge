"""Reselect submitted 2C9 and 3A4 banks after whole-member depth additions/swaps."""

import importlib.util
import itertools
import json
import os
from pathlib import Path

import numpy as np
import polars as pl


def load_scorer(path):
    spec = importlib.util.spec_from_file_location("official_scoring", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    scorer = load_scorer(inp / cfg["scorer"])
    keys = ["compound_id", "endpoint", "fold"]

    def score(y, low, high, prediction):
        return float(scorer.rae_soft_threshold_absolute_error(
            y, prediction, y_true_lower=low, y_true_upper=high
        ))

    def choose_subset(x, y, low, high):
        count = x.shape[1]
        winner = (float("inf"), count + 1, ())
        for size in range(1, count + 1):
            for subset in itertools.combinations(range(count), size):
                candidate = x[:, subset].mean(axis=1)
                result = (score(y, low, high, candidate), size, subset)
                if result < winner:
                    winner = result
        return winner[2]

    def caruana_weights(x, y, low, high, steps):
        counts = np.zeros(x.shape[1], dtype=int)
        prediction_sum = np.zeros(len(y), dtype=float)
        for step in range(steps):
            choices = [
                (score(y, low, high, (prediction_sum + x[:, col]) / (step + 1)), col)
                for col in range(x.shape[1])
            ]
            _, selected = min(choices)
            counts[selected] += 1
            prediction_sum += x[:, selected]
        return counts / counts.sum()

    def choose_caruana_steps(x, y, low, high, folds):
        choices = []
        for steps in cfg["caruana_steps"]:
            pred = np.empty(len(y), dtype=float)
            for inner in sorted(set(folds)):
                train = folds != inner
                weights = caruana_weights(x[train], y[train], low[train], high[train], steps)
                pred[~train] = x[~train] @ weights
            choices.append((score(y, low, high, pred), steps))
        return min(choices, key=lambda item: (item[0], item[1]))[1]

    rows, fits, predictions = [], [], []
    for target, setup in cfg["targets"].items():
        endpoint = setup["endpoint"]
        summary = json.loads((inp / setup["bank_summary"]).read_text())
        components = summary["config"]["components"]
        bank = pl.read_parquet(inp / setup["bank_oof"]).filter(
            pl.col("method") == setup["method"]
        ).sort(keys)
        table = bank.select(*keys, "y_true", "conf_low", "conf_high")
        original_names = [item["name"] for item in components]
        for item in components:
            atom = pl.read_parquet(inp / f"{target}-{item['name']}.parquet").filter(
                pl.col("endpoint") == endpoint
            )
            table = table.join(
                atom.select(*keys, pl.col("prediction").alias(item["name"])),
                on=keys, validate="1:1"
            ).sort(keys)
        assert table.height == setup["expected_rows"]

        def aligned(filename, name, depth=None, endpoint=endpoint, table=table):
            frame = pl.read_parquet(inp / filename).filter(pl.col("endpoint") == endpoint)
            if depth is not None:
                frame = frame.filter(pl.col("depth") == depth)
            assert frame.height == table.height
            return table.select(*keys).join(
                frame.select(*keys, pl.col("prediction").alias(name)),
                on=keys, validate="1:1"
            ).sort(keys)[name].to_numpy()

        candidates = {
            "chemeleon_d1_tabicl": aligned(cfg["chemeleon_tabicl"], "d1_tabicl"),
            "chemeleon_1_3_6_lgbm": np.mean([
                aligned(cfg["chemeleon_depths"], f"d{depth}", depth)
                for depth in (1, 3, 6)
            ], axis=0),
            "jeremy_1_2_4_lgbm": np.mean([
                aligned(setup["jeremy_oof"], f"j{depth}", depth)
                for depth in (1, 2, 4)
            ], axis=0),
        }
        y = table["y_true"].to_numpy()
        low = table["conf_low"].to_numpy()
        high = table["conf_high"].to_numpy()
        folds = table["fold"].to_numpy()
        scenarios = {"original_bank": (original_names, table.select(original_names).to_numpy())}
        for name in cfg["candidate_names"]:
            candidate = candidates[name]
            scenarios[f"add_{name}"] = (
                [*original_names, name],
                np.column_stack([table.select(original_names).to_numpy(), candidate]),
            )
            old = (
                setup["ft_member"] if name.startswith("jeremy") else setup["chemeleon_member"]
            )
            swapped = table.select(original_names).to_numpy().copy()
            swapped[:, original_names.index(old)] = candidate
            scenarios[f"swap_{old}__{name}"] = (original_names, swapped)
        for scenario, (names, x) in scenarios.items():
            prediction = np.empty(len(y), dtype=float)
            for held in range(5):
                train = folds != held
                if setup["method"] == "equal_subset":
                    subset = choose_subset(x[train], y[train], low[train], high[train])
                    prediction[~train] = x[~train][:, subset].mean(axis=1)
                    selected_names = [names[index] for index in subset]
                    hyperparameter = ""
                else:
                    steps = choose_caruana_steps(
                        x[train], y[train], low[train], high[train], folds[train]
                    )
                    weights = caruana_weights(
                        x[train], y[train], low[train], high[train], steps
                    )
                    prediction[~train] = x[~train] @ weights
                    selected_names = [names[index] for index in np.flatnonzero(weights > 0)]
                    hyperparameter = str(steps)
                candidate_selected = False
                if scenario.startswith("add_"):
                    candidate_selected = scenario.removeprefix("add_") in selected_names
                elif scenario.startswith("swap_"):
                    replaced_name = scenario.removeprefix("swap_").split("__", 1)[0]
                    candidate_selected = replaced_name in selected_names
                fits.append({
                    "target": target, "scenario": scenario, "held_fold": held,
                    "hyperparameter": hyperparameter,
                    "members": "+".join(selected_names),
                    "candidate_selected": candidate_selected,
                })
                print(f"target={target} scenario={scenario} fold={held}", flush=True)
            if scenario == "original_bank":
                difference = float(np.max(np.abs(prediction - bank["prediction"].to_numpy())))
                assert difference < 1e-9, (target, difference)
            rows.append({
                "target": target,
                "scenario": scenario,
                "n": len(y),
                "ST_RAE": score(y, low, high, prediction),
                "MAE": float(np.abs(y - prediction).mean()),
            })
            predictions.append(table.select(*keys).with_columns(
                pl.lit(target).alias("target"),
                pl.lit(scenario).alias("scenario"),
                pl.Series("prediction", prediction),
            ))
    pl.DataFrame(rows).write_csv(out / "metrics.csv")
    pl.DataFrame(fits).write_csv(out / "fits.csv")
    pl.concat(predictions).write_parquet(out / "oof.parquet")
    (out / "summary.json").write_text(json.dumps({"config": cfg}, indent=2) + "\n")


if __name__ == "__main__":
    main()
