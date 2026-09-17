"""Outer-cross-fit ensemble methods on an atomic regression prediction bank."""

import importlib.util
import itertools
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.optimize import minimize
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]
truth = ["y_true", "conf_low", "conf_high", "is_neighbor"]

spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)


def st_rae(y: np.ndarray, low: np.ndarray, high: np.ndarray, prediction: np.ndarray) -> float:
    return float(
        scorer.rae_soft_threshold_absolute_error(y, prediction, y_true_lower=low, y_true_upper=high)
    )


def interval_error(low: np.ndarray, high: np.ndarray, prediction: np.ndarray) -> float:
    return float((np.maximum(low - prediction, 0) + np.maximum(prediction - high, 0)).mean())


frames = []
for item in cfg["components"]:
    frame = pl.read_parquet(inp / item["file"]).sort(keys)
    if cfg.get("endpoints"):
        frame = frame.filter(pl.col("endpoint").is_in(cfg["endpoints"]))
    assert frame.height > 0 and frame.select(keys).n_unique() == frame.height
    frames.append(frame)
base = frames[0]
for frame in frames[1:]:
    for column in keys + truth:
        assert frame[column].equals(base[column]), column

names = [item["name"] for item in cfg["components"]]
x_all = np.column_stack([frame["prediction"].to_numpy() for frame in frames])
assert np.isfinite(x_all).all()
control = (
    pl.read_parquet(inp / cfg["current_control"]).filter(
        pl.col("endpoint").is_in(cfg.get("endpoints", []))
    )
    if cfg.get("endpoints")
    else pl.read_parquet(inp / cfg["current_control"])
)
control = control.sort(keys).select(*keys, pl.col("current").alias("current_control"))
base = base.join(control, on=keys, validate="1:1")
assert base["current_control"].is_finite().all()


def simplex_fit(x, y, low, high, loss_name, regularization):
    count = x.shape[1]

    def objective(weights):
        prediction = x @ weights
        if loss_name == "mse":
            loss = float(np.mean((prediction - y) ** 2))
        else:
            loss = interval_error(low, high, prediction)
        return loss + regularization * float(np.sum(weights**2))

    result = minimize(
        objective,
        np.full(count, 1 / count),
        method="SLSQP",
        bounds=[(0.0, 1.0)] * count,
        constraints={"type": "eq", "fun": lambda weights: weights.sum() - 1},
        options={"maxiter": 1000, "ftol": 1e-11},
    )
    if not result.success:
        raise RuntimeError(result.message)
    weights = np.maximum(result.x, 0)
    return weights / weights.sum()


def caruana_fit(x, y, low, high, steps):
    counts = np.zeros(x.shape[1], dtype=int)
    prediction_sum = np.zeros(x.shape[0], dtype=float)
    for step in range(steps):
        choices = []
        for column in range(x.shape[1]):
            prediction = (prediction_sum + x[:, column]) / (step + 1)
            choices.append((st_rae(y, low, high, prediction), column))
        _, selected = min(choices)
        counts[selected] += 1
        prediction_sum += x[:, selected]
    return counts / counts.sum()


def inner_select(x, y, low, high, folds, method, grid):
    rows = []
    for value in grid:
        prediction = np.empty(len(y), dtype=float)
        for held in sorted(set(folds)):
            train, valid = folds != held, folds == held
            if method == "ridge":
                model = make_pipeline(StandardScaler(), Ridge(alpha=value)).fit(x[train], y[train])
                prediction[valid] = model.predict(x[valid])
            elif method.startswith("simplex_"):
                weights = simplex_fit(
                    x[train],
                    y[train],
                    low[train],
                    high[train],
                    method.removeprefix("simplex_"),
                    value,
                )
                prediction[valid] = x[valid] @ weights
            else:
                weights = caruana_fit(x[train], y[train], low[train], high[train], value)
                prediction[valid] = x[valid] @ weights
        rows.append((st_rae(y, low, high, prediction), value))
    return min(rows, key=lambda item: (item[0], item[1]))[1]


prediction_rows = []
metric_rows = []
fold_metric_rows = []
fit_rows = []
indexed = base.with_row_index("row_index")

for (endpoint,), group in indexed.group_by("endpoint", maintain_order=True):
    if cfg.get("endpoints") and endpoint not in cfg["endpoints"]:
        continue
    group = group.sort(keys)
    ix = group["row_index"].to_numpy()
    x = x_all[ix]
    y = group["y_true"].to_numpy()
    low = group["conf_low"].to_numpy()
    high = group["conf_high"].to_numpy()
    folds = group["fold"].to_numpy()
    methods = {
        "current_control": group["current_control"].to_numpy(),
        "equal_all": x.mean(axis=1),
        "median_all": np.median(x, axis=1),
    }
    fitted_predictions = {
        name: np.empty(group.height, dtype=float)
        for name in [
            "best_single",
            "equal_subset",
            "ridge",
            "positive_linear",
            "simplex_mse",
            "simplex_interval",
            "caruana",
        ]
    }

    subsets = [
        subset
        for size in range(1, len(names) + 1)
        for subset in itertools.combinations(range(len(names)), size)
    ]
    for held in sorted(set(folds)):
        train, valid = folds != held, folds == held
        train_args = (y[train], low[train], high[train])

        singles = [(st_rae(*train_args, x[train, column]), column) for column in range(len(names))]
        _, single = min(singles)
        weights = np.zeros(len(names))
        weights[single] = 1
        fitted_predictions["best_single"][valid] = x[valid] @ weights
        fit_rows.append(
            {
                "endpoint": endpoint,
                "outer_fold": int(held),
                "method": "best_single",
                "hyperparameter": "",
                "members": names[single],
                "weights": json.dumps(weights.tolist()),
            }
        )

        choices = [
            (st_rae(*train_args, x[train][:, subset].mean(axis=1)), subset) for subset in subsets
        ]
        _, subset = min(choices, key=lambda item: (item[0], len(item[1]), item[1]))
        weights = np.zeros(len(names))
        weights[list(subset)] = 1 / len(subset)
        fitted_predictions["equal_subset"][valid] = x[valid] @ weights
        fit_rows.append(
            {
                "endpoint": endpoint,
                "outer_fold": int(held),
                "method": "equal_subset",
                "hyperparameter": "",
                "members": "+".join(names[i] for i in subset),
                "weights": json.dumps(weights.tolist()),
            }
        )

        alpha = inner_select(
            x[train], y[train], low[train], high[train], folds[train], "ridge", cfg["ridge_alpha"]
        )
        ridge = make_pipeline(StandardScaler(), Ridge(alpha=alpha)).fit(x[train], y[train])
        fitted_predictions["ridge"][valid] = ridge.predict(x[valid])
        fit_rows.append(
            {
                "endpoint": endpoint,
                "outer_fold": int(held),
                "method": "ridge",
                "hyperparameter": str(alpha),
                "members": "+".join(names),
                "weights": json.dumps(ridge[-1].coef_.tolist()),
            }
        )

        positive = LinearRegression(positive=True).fit(x[train], y[train])
        fitted_predictions["positive_linear"][valid] = positive.predict(x[valid])
        fit_rows.append(
            {
                "endpoint": endpoint,
                "outer_fold": int(held),
                "method": "positive_linear",
                "hyperparameter": "",
                "members": "+".join(names),
                "weights": json.dumps(positive.coef_.tolist()),
            }
        )

        for method in ["simplex_mse", "simplex_interval"]:
            regularization = inner_select(
                x[train],
                y[train],
                low[train],
                high[train],
                folds[train],
                method,
                cfg["simplex_lambda"],
            )
            weights = simplex_fit(
                x[train],
                y[train],
                low[train],
                high[train],
                method.removeprefix("simplex_"),
                regularization,
            )
            fitted_predictions[method][valid] = x[valid] @ weights
            fit_rows.append(
                {
                    "endpoint": endpoint,
                    "outer_fold": int(held),
                    "method": method,
                    "hyperparameter": str(regularization),
                    "members": "+".join(
                        name for name, weight in zip(names, weights, strict=True) if weight > 1e-8
                    ),
                    "weights": json.dumps(weights.tolist()),
                }
            )

        steps = inner_select(
            x[train],
            y[train],
            low[train],
            high[train],
            folds[train],
            "caruana",
            cfg["caruana_steps"],
        )
        weights = caruana_fit(x[train], y[train], low[train], high[train], steps)
        fitted_predictions["caruana"][valid] = x[valid] @ weights
        fit_rows.append(
            {
                "endpoint": endpoint,
                "outer_fold": int(held),
                "method": "caruana",
                "hyperparameter": str(steps),
                "members": "+".join(
                    name for name, weight in zip(names, weights, strict=True) if weight > 0
                ),
                "weights": json.dumps(weights.tolist()),
            }
        )

    methods.update(fitted_predictions)
    for method, prediction in methods.items():
        assert np.isfinite(prediction).all()
        prediction_rows.append(
            group.select(*keys, *truth).with_columns(
                pl.lit(method).alias("method"), pl.Series("prediction", prediction)
            )
        )
        metric_rows.append(
            {
                "endpoint": endpoint,
                "method": method,
                "ST_RAE": st_rae(y, low, high, prediction),
                "MAE": float(np.abs(y - prediction).mean()),
                "Spearman": float(spearmanr(y, prediction).statistic),
            }
        )
        for held in sorted(set(folds)):
            mask = folds == held
            fold_metric_rows.append(
                {
                    "endpoint": endpoint,
                    "fold": int(held),
                    "method": method,
                    "ST_RAE": st_rae(y[mask], low[mask], high[mask], prediction[mask]),
                }
            )

metrics = pl.DataFrame(metric_rows).sort("endpoint", "ST_RAE")
macro = metrics.group_by("method").agg(pl.col("ST_RAE", "MAE", "Spearman").mean()).sort("ST_RAE")
pl.concat(prediction_rows).sort("endpoint", "method", "fold", "compound_id").write_parquet(
    out / "oof.parquet"
)
metrics.write_csv(out / "metrics.csv")
macro.write_csv(out / "macro.csv")
pl.DataFrame(fold_metric_rows).sort("endpoint", "method", "fold").write_csv(
    out / "fold-metrics.csv"
)
pl.DataFrame(fit_rows).sort("endpoint", "method", "outer_fold").write_csv(out / "fits.csv")
(out / "bank.json").write_text(
    json.dumps({"atomic_members": names, "control_is_member": False}, indent=2)
)
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "rows": base.height,
            "primary_selection_scope": "endpoint",
            "macro_use": "descriptive_only",
            "endpoint_metrics": metrics.to_dicts(),
            "macro": macro.to_dicts(),
        },
        indent=2,
    )
)
print(macro)
print(metrics)
