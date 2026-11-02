"""Test whether analog-proxy OOF rows explain why nested blends failed on the Live board.

The official test set is an analog expansion around top training hits, so OOF rows that are
analogs of training anchors are a closer stand-in for it than the full cohort. Per endpoint:
(1) score the sent blends on all, analog and other rows next to their observed Live deltas,
(2) repeat the nested anchored selection with the analog rows as the objective, and
(3) scan global shifts of the incumbent on analog and other rows.
"""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import regression_bank_selection as sel
import regression_grand_bank as grand


def delta(y, low, high, blend, incumbent, rows):
    return sel.st_rae(y[rows], low[rows], high[rows], blend[rows]) - sel.st_rae(
        y[rows], low[rows], high[rows], incumbent[rows]
    )


def bootstrap(rng, y, low, high, blend, incumbent, rows, draws):
    index = np.flatnonzero(rows)
    values = []
    for _ in range(draws):
        pick = rng.choice(index, len(index))
        values.append(
            sel.st_rae(y[pick], low[pick], high[pick], blend[pick])
            - sel.st_rae(y[pick], low[pick], high[pick], incumbent[pick])
        )
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def nested(x, y, low, high, folds, incumbent, rows_for_fit, cfg):
    blend = np.empty(len(y))
    for outer in sorted(set(folds)):
        fit, held = folds != outer, folds == outer
        rows = np.flatnonzero(fit)
        rows = rows_for_fit(rows)
        weights, _ = sel.procedure_weights(
            x[rows], y[rows], low[rows], high[rows], folds[rows], "greedy_anchored", cfg
        )
        blend[held] = sel.moment_match(x[held] @ weights, incumbent[held])
    return blend


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    rng = np.random.default_rng(cfg["seed"])
    root = Path(cfg["runs_root"])
    analog = set(pl.read_parquet(inp / "analog-pairs.parquet")["analog_id"].to_list())
    test_ids = (
        pl.read_parquet(inp / "test-mapping.parquet").sort("compound_id")["compound_id"].to_list()
    )
    cohorts, truth = {}, {}
    for endpoint, control in cfg["controls"].items():
        frame = pl.read_parquet(inp / control["file"])
        for column, value in control["filter"].items():
            frame = frame.filter(pl.col(column) == value)
        frame = frame.filter(pl.col("endpoint") == endpoint).sort("compound_id")
        assert frame.height == control["rows"]
        cohorts[endpoint], truth[endpoint] = frame.select("compound_id", "fold"), frame
    bank = grand.collect_oof(root, cohorts, grand.collect_blind(root, test_ids))
    cypmol = pl.read_parquet(inp / "cypmol-oof.parquet")

    yardstick, selection, shifts = [], [], []
    for endpoint, frame in truth.items():
        ids = frame["compound_id"].to_list()
        y, low, high, folds = (
            frame[c].to_numpy() for c in ("y_true", "conf_low", "conf_high", "fold")
        )
        incumbent = frame["prediction"].to_numpy()
        is_analog = np.array([i in analog for i in ids])
        x = np.column_stack([incumbent] + [m["oof"] for m in bank[endpoint]])
        everything = np.ones(len(y), bool)

        blends = {"286a7b27": nested(x, y, low, high, folds, incumbent, lambda r: r, cfg)}
        weight = cfg["cypmol_weights"].get(endpoint)
        if weight:
            wide = (
                cypmol.filter(pl.col("endpoint") == endpoint)
                .pivot(on="model", index="compound_id", values="prediction")
                .sort("compound_id")
            )
            member = (
                frame.select("compound_id")
                .join(wide, on="compound_id", how="left")
                .drop("compound_id")
                .to_numpy()
                .mean(axis=1)
            )
            mixed = np.empty(len(y))
            for fold in sorted(set(folds)):
                held = folds == fold
                mixed[held] = sel.moment_match(
                    (1 - weight) * incumbent[held] + weight * member[held], incumbent[held]
                )
            blends["3dc232c5"] = mixed
        for name, blend in blends.items():
            record = {"endpoint": endpoint, "blend": name, "analog_rows": int(is_analog.sum())}
            for tag, rows in (("all", everything), ("analog", is_analog), ("other", ~is_analog)):
                record[f"delta_{tag}"] = delta(y, low, high, blend, incumbent, rows)
            lo95, hi95 = bootstrap(rng, y, low, high, blend, incumbent, is_analog, cfg["bootstrap"])
            record["analog_ci"] = [lo95, hi95]
            record["live_delta"] = cfg["live_delta"].get(endpoint, {}).get(name)
            yardstick.append(record)

        for duplication, label in ((1, "analog_only"), (3, "analog_x3")):

            def rows_for_fit(rows, duplication=duplication, is_analog=is_analog):
                if duplication == 1:
                    return rows[is_analog[rows]]
                return np.concatenate([rows] + [rows[is_analog[rows]]] * (duplication - 1))

            blend = nested(x, y, low, high, folds, incumbent, rows_for_fit, cfg)
            lo95, hi95 = bootstrap(rng, y, low, high, blend, incumbent, is_analog, cfg["bootstrap"])
            selection.append(
                {
                    "endpoint": endpoint,
                    "objective": label,
                    "delta_analog": delta(y, low, high, blend, incumbent, is_analog),
                    "analog_ci": [lo95, hi95],
                    "delta_other": delta(y, low, high, blend, incumbent, ~is_analog),
                }
            )
        for shift in cfg["shifts"]:
            shifted = incumbent + shift
            shifts.append(
                {
                    "endpoint": endpoint,
                    "shift": shift,
                    "analog": sel.st_rae(
                        y[is_analog], low[is_analog], high[is_analog], shifted[is_analog]
                    ),
                    "other": sel.st_rae(
                        y[~is_analog], low[~is_analog], high[~is_analog], shifted[~is_analog]
                    ),
                }
            )
        print(json.dumps(yardstick[-1]), flush=True)
    (out / "yardstick.json").write_text(json.dumps(yardstick, indent=2) + "\n")
    (out / "analog-selection.json").write_text(json.dumps(selection, indent=2) + "\n")
    pl.DataFrame(shifts).write_csv(out / "shift-scan.csv")
    print(pl.DataFrame(selection))


if __name__ == "__main__":
    main()
