"""Conditional parent-group bootstrap for selected fixed depth scenarios."""

import importlib.util
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
    labels = (
        pl.read_parquet(inp / cfg["labels"])
        .filter(pl.col("endpoint") == cfg.get("endpoint", "CYP2D6_pIC50_direct_inhibition"))
        .select(*keys, "y_true", "conf_low", "conf_high")
    )
    groups = (
        pl.read_parquet(inp / cfg["assignments"])
        .filter(pl.col("strategy") == "series")
        .select("compound_id", "group")
    )
    table = labels.join(groups, on="compound_id", validate="m:1").sort(keys)
    scenarios = pl.read_parquet(inp / cfg["scenario_oof"])
    names = [cfg["control"], *cfg["candidates"]]
    for name in names:
        part = scenarios.filter(pl.col("scenario") == name).select(
            *keys, pl.col("prediction").alias(name)
        )
        table = table.join(part, on=keys, validate="1:1")
        assert table.height == cfg.get("expected_rows", 1457)
    y = table["y_true"].to_numpy()
    low = table["conf_low"].to_numpy()
    high = table["conf_high"].to_numpy()
    group = table["group"].to_numpy()
    unique_groups = np.unique(group)
    indices = [np.flatnonzero(group == value) for value in unique_groups]
    values = {name: table[name].to_numpy() for name in names}

    def score(name, selected):
        return float(
            scorer.rae_soft_threshold_absolute_error(
                y[selected],
                values[name][selected],
                y_true_lower=low[selected],
                y_true_upper=high[selected],
            )
        )

    rng = np.random.default_rng(cfg["seed"])
    differences = {name: [] for name in cfg["candidates"]}
    for _ in range(cfg["bootstrap_replicates"]):
        sample = rng.integers(len(unique_groups), size=len(unique_groups))
        selected = np.concatenate([indices[index] for index in sample])
        base = score(cfg["control"], selected)
        for name in cfg["candidates"]:
            differences[name].append(score(name, selected) - base)
    rows = []
    whole = np.arange(len(y))
    base_score = score(cfg["control"], whole)
    for name in cfg["candidates"]:
        delta = np.asarray(differences[name])
        rows.append(
            {
                "scenario": name,
                "control_ST_RAE": base_score,
                "candidate_ST_RAE": score(name, whole),
                "delta_ST_RAE": score(name, whole) - base_score,
                "bootstrap_025": float(np.quantile(delta, 0.025)),
                "bootstrap_975": float(np.quantile(delta, 0.975)),
                "fraction_better": float(np.mean(delta < 0)),
                "groups": len(unique_groups),
                "replicates": cfg["bootstrap_replicates"],
            }
        )
    pl.DataFrame(rows).write_csv(out / "bootstrap.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "limitation": (
                    "Conditional on fitted OOF models; ignores post-screen selection, "
                    "refitting, and distribution shift"
                ),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
