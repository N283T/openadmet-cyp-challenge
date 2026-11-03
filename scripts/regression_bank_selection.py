"""Nested reselection over the incumbent and every frozen-embedding member on the shared OOF cohort.

Every screen reuses the incumbent's own OOF rows and series folds, so the incumbent and all member
columns of one endpoint are directly comparable. Each procedure is chosen on four folds and scored
on the fifth: the incumbent alone, the best single column, and two greedy with-replacement
selections (one anchored on the incumbent, one free) whose step count is picked by an inner
leave-one-fold-out loop. Scores are ST-RAE pooled over the five outer folds.
"""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

KEYS = ["compound_id", "endpoint", "fold"]
INCUMBENT = "incumbent"


def st_rae(y, low, high, p):
    def error(q):
        return np.maximum(low - q, 0) + np.maximum(q - high, 0)

    return float(error(p).sum() / error(y.mean()).sum())


def moment_match(p, reference):
    return reference.mean() + (p - p.mean()) / p.std(ddof=1) * reference.std(ddof=1)


def greedy_path(x, y, low, high, k_max, anchor, anchor_count):
    """Counts after 0..k_max greedy additions; the anchor column starts with `anchor_count`."""
    counts = np.zeros(x.shape[1])
    if anchor_count:
        counts[anchor] = anchor_count
    total, running = counts.sum(), x @ counts
    path = [counts.copy()]
    for _ in range(k_max):
        scores = [
            st_rae(y, low, high, (running + x[:, c]) / (total + 1)) for c in range(x.shape[1])
        ]
        pick = int(np.argmin(scores))
        counts[pick] += 1
        running += x[:, pick]
        total += 1
        path.append(counts.copy())
    return path


def choose_steps(x, y, low, high, folds, k_grid, anchor, anchor_count):
    k_max = max(k_grid)
    pooled = {k: np.empty(len(y)) for k in k_grid}
    for inner in sorted(set(folds)):
        fit = folds != inner
        path = greedy_path(x[fit], y[fit], low[fit], high[fit], k_max, anchor, anchor_count)
        for k in k_grid:
            weights = path[k] / path[k].sum()
            pooled[k][~fit] = x[~fit] @ weights
    scored = [(st_rae(y, low, high, pooled[k]), k) for k in k_grid]
    return min(scored)[1]


def procedure_weights(x, y, low, high, folds, name, cfg):
    columns = x.shape[1]
    one_hot = np.zeros(columns)
    if name == "incumbent":
        one_hot[0] = 1.0
        return one_hot, 0
    if name == "best_single":
        scores = [st_rae(y, low, high, x[:, c]) for c in range(columns)]
        one_hot[int(np.argmin(scores))] = 1.0
        return one_hot, 1
    anchored = name == "greedy_anchored"
    anchor_count = cfg["anchor_count"] if anchored else 0
    k_grid = cfg["k_grid"]
    k = choose_steps(x, y, low, high, folds, k_grid, 0, anchor_count)
    path = greedy_path(x, y, low, high, k, 0, anchor_count)
    return path[k] / path[k].sum(), k


def build_table(inp, cfg, endpoint):
    control = cfg["controls"][endpoint]
    frame = pl.read_parquet(inp / control["file"])
    for column, value in control["filter"].items():
        frame = frame.filter(pl.col(column) == value)
    frame = frame.filter(pl.col("endpoint") == endpoint)
    table = frame.select(
        *KEYS, "y_true", "conf_low", "conf_high", pl.col("prediction").alias(INCUMBENT)
    )
    table = table.sort(KEYS)
    assert table.height == control["rows"] == table.select(KEYS).unique().height
    for screen, file in cfg["screens"].items():
        wide = (
            pl.read_parquet(inp / file)
            .filter(pl.col("endpoint") == endpoint)
            .pivot(on="model", index=KEYS, values="prediction")
        )
        wide = wide.rename({c: f"{screen}:{c}" for c in wide.columns if c not in KEYS})
        table = table.join(wide, on=KEYS, how="left", validate="1:1")
    assert table.null_count().sum_horizontal().item() == 0
    return table.sort(KEYS)


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    rng = np.random.default_rng(cfg["seed"])
    procedures = ["incumbent", "best_single", "greedy_anchored", "greedy_free"]
    standalone, nested, boot, finals = [], [], [], []
    for endpoint in cfg["controls"]:
        table = build_table(inp, cfg, endpoint)
        names = [c for c in table.columns if c not in KEYS + ["y_true", "conf_low", "conf_high"]]
        x = table.select(names).to_numpy()
        y, low, high = (table[c].to_numpy() for c in ("y_true", "conf_low", "conf_high"))
        folds = table["fold"].to_numpy()
        base = st_rae(y, low, high, x[:, 0])
        for column, name in enumerate(names):
            standalone.append(
                {
                    "endpoint": endpoint,
                    "column": name,
                    "st_rae": st_rae(y, low, high, x[:, column]),
                    "incumbent": base,
                }
            )
        raw = {p: np.empty(len(y)) for p in procedures}
        matched = {p: np.empty(len(y)) for p in procedures}
        for outer in sorted(set(folds)):
            fit, held = folds != outer, folds == outer
            for procedure in procedures:
                weights, k = procedure_weights(
                    x[fit], y[fit], low[fit], high[fit], folds[fit], procedure, cfg
                )
                raw[procedure][held] = x[held] @ weights
                matched[procedure][held] = moment_match(raw[procedure][held], x[held, 0])
                nested.append(
                    {
                        "endpoint": endpoint,
                        "outer_fold": int(outer),
                        "procedure": procedure,
                        "steps": k,
                        "members": json.dumps(
                            {names[c]: round(float(w), 3) for c, w in enumerate(weights) if w > 0}
                        ),
                    }
                )
        for procedure in procedures:
            wins = sum(
                st_rae(
                    y[folds == f], low[folds == f], high[folds == f], matched[procedure][folds == f]
                )
                < st_rae(y[folds == f], low[folds == f], high[folds == f], x[folds == f, 0])
                for f in sorted(set(folds))
            )
            deltas = []
            for _ in range(cfg["bootstrap"]):
                pick = rng.integers(0, len(y), len(y))
                deltas.append(
                    st_rae(y[pick], low[pick], high[pick], matched[procedure][pick])
                    - st_rae(y[pick], low[pick], high[pick], x[pick, 0])
                )
            boot.append(
                {
                    "endpoint": endpoint,
                    "procedure": procedure,
                    "incumbent": base,
                    "raw": st_rae(y, low, high, raw[procedure]),
                    "moment_matched": st_rae(y, low, high, matched[procedure]),
                    "delta_matched": st_rae(y, low, high, matched[procedure]) - base,
                    "delta_lo": float(np.quantile(deltas, 0.025)),
                    "delta_hi": float(np.quantile(deltas, 0.975)),
                    "folds_won": int(wins),
                }
            )
            weights, k = procedure_weights(x, y, low, high, folds, procedure, cfg)
            finals.append(
                {
                    "endpoint": endpoint,
                    "procedure": procedure,
                    "steps": k,
                    "weights": json.dumps(
                        {names[c]: round(float(w), 4) for c, w in enumerate(weights) if w > 0}
                    ),
                }
            )
        print(endpoint, flush=True)
    pl.DataFrame(standalone).write_csv(out / "standalone.csv")
    pl.DataFrame(nested).write_csv(out / "nested-choices.csv")
    pl.DataFrame(boot).write_csv(out / "nested-scores.csv")
    pl.DataFrame(finals).write_csv(out / "final-weights.csv")
    print(pl.DataFrame(boot))


if __name__ == "__main__":
    main()
