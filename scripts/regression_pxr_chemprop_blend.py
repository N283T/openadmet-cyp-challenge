"""Evaluate light endpoint-specific blends of strict PXR-Chemprop regression readouts."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


keys = ["compound_id", "endpoint", "fold"]
truth = ["y_true", "conf_low", "conf_high", "is_neighbor"]
components = {name: read(name).sort(keys) for name in cfg["component_names"]}
pxr = read("pxr").sort(keys)
base = components[cfg["component_names"][0]]
assert base.height == pxr.height == 6384
for frame in [*components.values(), pxr]:
    assert not frame.select(keys).is_duplicated().any()
    for column in keys + truth:
        assert frame[column].equals(base[column]), (column, frame.shape)

top2 = (
    read("caruana_top2")
    .filter(
        (pl.col("strategy") == "crossfold_topk")
        & (pl.col("weighting") == "weighted")
        & (pl.col("k") == 2)
    )
    .sort(keys)
)
assert top2.height == base.filter(pl.col("endpoint") == cfg["a4_endpoint"]).height

spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
assert spec is not None and spec.loader is not None
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)


def score(frame: pl.DataFrame, prediction: np.ndarray) -> float:
    return float(
        scorer.rae_soft_threshold_absolute_error(
            frame["y_true"].to_numpy(),
            prediction,
            y_true_lower=frame["conf_low"].to_numpy(),
            y_true_upper=frame["conf_high"].to_numpy(),
        )
    )


baseline_parts = {
    endpoint: names for endpoint, names in cfg["endpoint_components"].items()
}
baseline = np.empty(base.height, dtype=float)
for endpoint, names in baseline_parts.items():
    mask = base["endpoint"].to_numpy() == endpoint
    values = [components[name]["prediction"].to_numpy()[mask] for name in names]
    if endpoint == cfg["a4_endpoint"]:
        values.insert(0, top2["caruana_prediction"].to_numpy())
    baseline[mask] = np.mean(values, axis=0)

pxr_prediction = pxr["prediction"].to_numpy()
weights = cfg["weights"]
metric_rows = []
fold_rows = []
crossfit_rows = []
bootstrap_rows = []
rng = np.random.default_rng(cfg["bootstrap_seed"])

for (endpoint,), group in base.with_row_index("row_index").group_by(
    "endpoint", maintain_order=True
):
    ix = group["row_index"].to_numpy()
    folds = group["fold"].to_numpy()
    old = baseline[ix]
    new = pxr_prediction[ix]
    endpoint_scores = {}
    for weight in weights:
        prediction = (1 - weight) * old + weight * new
        endpoint_scores[weight] = score(group, prediction)
        metric_rows.append(
            {
                "endpoint": endpoint,
                "weight": weight,
                "ST_RAE": endpoint_scores[weight],
                "MAE": float(np.abs(group["y_true"].to_numpy() - prediction).mean()),
                "Spearman": float(spearmanr(group["y_true"].to_numpy(), prediction).statistic),
            }
        )
        for fold in sorted(set(folds)):
            valid = folds == fold
            fold_rows.append(
                {
                    "endpoint": endpoint,
                    "weight": weight,
                    "fold": int(fold),
                    "ST_RAE": score(group.filter(pl.Series(valid)), prediction[valid]),
                }
            )

    crossfit_prediction = np.empty(group.height, dtype=float)
    selections = []
    for fold in sorted(set(folds)):
        train = folds != fold
        valid = ~train
        train_frame = group.filter(pl.Series(train))
        scored = [
            (score(train_frame, ((1 - w) * old + w * new)[train]), w) for w in weights
        ]
        _, selected = min(scored, key=lambda item: (item[0], item[1]))
        crossfit_prediction[valid] = ((1 - selected) * old + selected * new)[valid]
        selections.append(float(selected))
    crossfit_score = score(group, crossfit_prediction)
    best_weight = min(weights, key=lambda w: (endpoint_scores[w], w))
    best_prediction = (1 - best_weight) * old + best_weight * new
    delta = endpoint_scores[best_weight] - endpoint_scores[0.0]
    bootstrap = []
    for _ in range(cfg["bootstrap_replicates"]):
        sample = rng.integers(0, group.height, group.height)
        sampled = group[sample]
        bootstrap.append(
            score(sampled, best_prediction[sample]) - score(sampled, old[sample])
        )
    bootstrap_rows.append(
        {
            "endpoint": endpoint,
            "best_full_oof_weight": float(best_weight),
            "baseline_ST_RAE": endpoint_scores[0.0],
            "best_full_oof_ST_RAE": endpoint_scores[best_weight],
            "delta": delta,
            "delta_ci_low": float(np.quantile(bootstrap, 0.025)),
            "delta_ci_high": float(np.quantile(bootstrap, 0.975)),
            "crossfit_ST_RAE": crossfit_score,
            "fold_selected_weights": ";".join(f"{value:g}" for value in selections),
            "prediction_pearson": float(np.corrcoef(old, new)[0, 1]),
            "residual_pearson": float(
                np.corrcoef(
                    old - group["y_true"].to_numpy(),
                    new - group["y_true"].to_numpy(),
                )[0, 1]
            ),
        }
    )
    crossfit_rows.append(
        group.select(keys + truth).with_columns(
            pl.Series("baseline_prediction", old),
            pl.Series("pxr_prediction", new),
            pl.Series("crossfit_prediction", crossfit_prediction),
        )
    )

metrics = pl.DataFrame(metric_rows).sort("endpoint", "weight")
summary = pl.DataFrame(bootstrap_rows).sort("endpoint")
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(fold_rows).sort("endpoint", "weight", "fold").write_csv(out / "fold-metrics.csv")
summary.write_csv(out / "summary.csv")
pl.concat(crossfit_rows).sort(keys).write_parquet(out / "predictions.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "interpretation": (
                "Full-OOF weights are deployment diagnostics; fold-selected crossfit scores "
                "evaluate weight selection using only the other four folds."
            ),
            "rows": summary.to_dicts(),
        },
        indent=2,
    )
    + "\n"
)
print(summary, flush=True)
