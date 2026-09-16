"""Screen uncertainty-weighted ranking, molecular kernels, and classical D6 readouts."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from scipy.stats import spearmanr
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import HuberRegressor, LogisticRegression, RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import SplineTransformer, StandardScaler
from threadpoolctl import threadpool_limits

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoint = cfg["endpoint"]
keys = ["compound_id", "endpoint", "fold"]


def read(index: int) -> pl.DataFrame:
    paths = list(inp.glob(f"{index:03d}-*"))
    if len(paths) != 1:
        raise ValueError(f"Expected one input at index {index}, found {paths}")
    return pl.read_parquet(paths[0])


def st_rae(frame: pl.DataFrame, prediction: np.ndarray) -> float:
    lower = frame["conf_low"].to_numpy()
    upper = frame["conf_high"].to_numpy()
    y = frame["y_true"].to_numpy()
    error = np.maximum(lower - prediction, 0) + np.maximum(prediction - upper, 0)
    center = np.full_like(y, y.mean())
    denominator = (np.maximum(lower - center, 0) + np.maximum(center - upper, 0)).sum()
    if denominator <= 0:
        raise ValueError("Non-positive ST-RAE denominator")
    return float(error.sum() / denominator)


def confident_pairs(
    lower: np.ndarray, upper: np.ndarray, maximum: int, seed: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    first, second = np.triu_indices(len(lower), k=1)
    first_higher = lower[first] > upper[second]
    second_higher = lower[second] > upper[first]
    high = np.concatenate([first[first_higher], second[second_higher]])
    low = np.concatenate([second[first_higher], first[second_higher]])
    gap = lower[high] - upper[low]
    total = len(high)
    if total > maximum:
        probability = gap / gap.sum()
        chosen = np.random.default_rng(seed).choice(
            total, size=maximum, replace=False, p=probability
        )
        high, low, gap = high[chosen], low[chosen], gap[chosen]
    weight = np.clip(gap / np.median(gap), 0.25, 4.0)
    return high, low, weight, total


def positive_affine(
    train_score: np.ndarray, train_y: np.ndarray, score: np.ndarray
) -> tuple[np.ndarray, float, float]:
    variance = float(np.var(train_score))
    if variance <= 0:
        raise ValueError("Constant utility")
    slope = float(np.cov(train_score, train_y, ddof=0)[0, 1] / variance)
    if slope <= 0:
        raise ValueError(f"Non-positive calibration slope: {slope}")
    intercept = float(train_y.mean() - slope * train_score.mean())
    return intercept + slope * score, slope, intercept


def tanimoto(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    intersection = left @ right.T
    denominator = left.sum(axis=1)[:, None] + right.sum(axis=1)[None, :] - intersection
    return np.divide(
        intersection,
        denominator,
        out=np.zeros_like(intersection, dtype=np.float32),
        where=denominator > 0,
    )


def fit_interval_pair(
    x_train: np.ndarray,
    x_all: np.ndarray,
    y_train: np.ndarray,
    lower_train: np.ndarray,
    upper_train: np.ndarray,
    high: np.ndarray,
    low: np.ndarray,
    pair_weight: np.ndarray,
    pair_strength: float,
    seed: int,
    path: Path,
) -> np.ndarray:
    torch.manual_seed(seed)
    model = torch.nn.Linear(x_train.shape[1], 1)
    optimizer = torch.optim.AdamW(model.parameters(), lr=0.02, weight_decay=1e-3)
    tx = torch.from_numpy(x_train.astype(np.float32))
    ty = torch.from_numpy(y_train.astype(np.float32))
    tlower = torch.from_numpy(lower_train.astype(np.float32))
    tupper = torch.from_numpy(upper_train.astype(np.float32))
    thigh = torch.from_numpy(high)
    tlow = torch.from_numpy(low)
    tweight = torch.from_numpy(pair_weight.astype(np.float32))
    for _ in range(cfg["interval_epochs"]):
        optimizer.zero_grad(set_to_none=True)
        score = model(tx).squeeze(-1)
        interval_loss = (torch.relu(tlower - score) + torch.relu(score - tupper)).mean()
        point_loss = torch.nn.functional.smooth_l1_loss(score, ty)
        rank_loss = (tweight * torch.nn.functional.softplus(-(score[thigh] - score[tlow]))).mean()
        loss = interval_loss + 0.1 * point_loss + pair_strength * rank_loss
        loss.backward()
        optimizer.step()
    torch.save(model.state_dict(), path)
    with torch.no_grad():
        return model(torch.from_numpy(x_all.astype(np.float32))).squeeze(-1).numpy().astype(float)


assignments, labels, scores, classic, nested, latest, current = [read(index) for index in range(7)]
labels = labels.filter(pl.col("endpoint") == endpoint)
score_columns = sorted(column for column in scores.columns if column.startswith("score_"))
fingerprint_columns = sorted(column for column in classic.columns if column.startswith("morgan_"))
if not fingerprint_columns:
    fingerprint_columns = sorted(column for column in classic.columns if column.startswith("fp_"))
if (len(score_columns), len(fingerprint_columns)) != (47, 2048):
    raise ValueError("Unexpected feature dimensions")
static = scores.select("compound_id", *score_columns).join(
    classic.select("compound_id", *fingerprint_columns), on="compound_id", validate="1:1"
)
base = labels.join(assignments, on="compound_id", validate="1:1").sort("compound_id")
latest = latest.filter(pl.col("endpoint") == endpoint).sort(keys)
current = current.filter((pl.col("endpoint") == endpoint) & (pl.col("component_count") == 2)).sort(
    keys
)
for column in keys + ["y_true", "conf_low", "conf_high"]:
    if not latest[column].equals(current[column]):
        raise ValueError(f"Baseline mismatch in {column}")

rows = []
fit_rows = []
memberships = []
threadpool_limits(limits=cfg["threads"])
torch.set_num_threads(cfg["threads"])
for fold in range(5):
    data = (
        base.join(static, on="compound_id", validate="1:1")
        .join(
            nested.filter(pl.col("outer_fold") == fold).select(
                "compound_id", *[f"aux_{index}" for index in range(4)]
            ),
            on="compound_id",
            validate="1:1",
        )
        .sort("compound_id")
    )
    validation = (data["fold"].to_numpy() == fold) & ~data["is_anchor"].to_numpy()
    train = ~validation
    y = data["value"].to_numpy().astype(float)
    lower = data["conf_low"].to_numpy().astype(float)
    upper = data["conf_high"].to_numpy().astype(float)
    score_x = data.select(score_columns + [f"aux_{index}" for index in range(4)]).to_numpy()
    fp_x = data.select(fingerprint_columns).to_numpy().astype(np.float32)
    if not np.isfinite(score_x).all() or not np.isfinite(fp_x).all():
        raise ValueError("Non-finite features")
    scaler = StandardScaler().fit(score_x[train])
    scaled = scaler.transform(score_x).astype(np.float32)
    fold_dir = out / f"fold{fold}"
    fold_dir.mkdir()
    joblib.dump(scaler, fold_dir / "scores-scaler.joblib")
    high, low, pair_weight, available = confident_pairs(
        lower[train], upper[train], cfg["maximum_pairs"], 42 + fold
    )
    memberships.append(
        data.select("compound_id", "endpoint", "fold", "group", "is_anchor").with_columns(
            pl.lit(fold).alias("outer_fold"), pl.Series("validation", validation)
        )
    )

    predictions: dict[tuple[str, str], np.ndarray] = {}
    calibrations: dict[tuple[str, str], tuple[float | None, float | None]] = {}

    difference = scaled[train][high] - scaled[train][low]
    pair_x = np.vstack([difference, -difference])
    pair_y = np.concatenate([np.ones(len(high)), np.zeros(len(high))])
    weights = np.concatenate([pair_weight, pair_weight])
    pair_model = LogisticRegression(C=1.0, max_iter=3000, random_state=42 + fold)
    pair_model.fit(pair_x, pair_y, sample_weight=weights)
    pair_utility = scaled @ pair_model.coef_.reshape(-1)
    prediction, slope, intercept = positive_affine(pair_utility[train], y[train], pair_utility)
    key = ("scores51", "pairwise_logistic_weighted")
    predictions[key] = prediction[validation]
    calibrations[key] = (slope, intercept)
    joblib.dump(pair_model, fold_dir / "pairwise-logistic-weighted.joblib")

    y_mean = y[train].mean()
    y_scale = y[train].std()
    standardized_y = (y[train] - y_mean) / y_scale
    standardized_lower = (lower[train] - y_mean) / y_scale
    standardized_upper = (upper[train] - y_mean) / y_scale
    for strength in cfg["interval_pair_weights"]:
        method = f"interval_pair{int(round(strength * 100)):03d}"
        utility = fit_interval_pair(
            scaled[train],
            scaled,
            standardized_y,
            standardized_lower,
            standardized_upper,
            high,
            low,
            pair_weight,
            strength,
            42 + fold,
            fold_dir / f"{method}.pt",
        )
        predictions[("scores51", method)] = (y_mean + y_scale * utility)[validation]
        calibrations[("scores51", method)] = (y_scale, y_mean)

    train_fp = fp_x[train]
    train_kernel = tanimoto(train_fp, train_fp)
    validation_kernel = tanimoto(fp_x[validation], train_fp)
    for alpha in cfg["tanimoto_alphas"]:
        suffix = f"{alpha:g}".replace(".", "p")
        method = f"tanimoto_krr_a{suffix}"
        model = KernelRidge(alpha=alpha, kernel="precomputed")
        model.fit(train_kernel, y[train] - y_mean)
        predictions[("morgan2048", method)] = y_mean + model.predict(validation_kernel)
        calibrations[("morgan2048", method)] = (None, y_mean)
        joblib.dump(model, fold_dir / f"{method}.joblib")

    for bins in cfg["ordinal_bins"]:
        thresholds = np.quantile(y[train], np.arange(1, bins) / bins)
        models = []
        train_utility = np.zeros(int(train.sum()))
        validation_utility = np.zeros(int(validation.sum()))
        for threshold in thresholds:
            model = LogisticRegression(C=1.0, max_iter=3000, random_state=42 + fold)
            model.fit(scaled[train], y[train] > threshold)
            train_utility += model.predict_proba(scaled[train])[:, 1]
            validation_utility += model.predict_proba(scaled[validation])[:, 1]
            models.append(model)
        prediction, slope, intercept = positive_affine(train_utility, y[train], validation_utility)
        method = f"ordinal_logit{bins}"
        predictions[("scores51", method)] = prediction
        calibrations[("scores51", method)] = (slope, intercept)
        joblib.dump((thresholds, models), fold_dir / f"{method}.joblib")

    gam = make_pipeline(
        StandardScaler(),
        SplineTransformer(n_knots=4, degree=2, knots="quantile", include_bias=False),
        RidgeCV(alphas=np.asarray([0.1, 1.0, 10.0, 100.0])),
    )
    gam.fit(score_x[train], y[train])
    predictions[("scores51", "spline_gam")] = gam.predict(score_x[validation])
    calibrations[("scores51", "spline_gam")] = (None, None)
    joblib.dump(gam, fold_dir / "spline-gam.joblib")

    huber = HuberRegressor(epsilon=1.35, alpha=1.0, max_iter=3000)
    huber.fit(scaled[train], y[train])
    predictions[("scores51", "huber")] = huber.predict(scaled[validation])
    calibrations[("scores51", "huber")] = (None, None)
    joblib.dump(huber, fold_dir / "huber.joblib")

    held = data.filter(pl.Series(validation))
    for (view, method), prediction in predictions.items():
        row = held.select(
            "compound_id",
            "endpoint",
            "fold",
            "group",
            "is_neighbor",
            "conf_low",
            "conf_high",
            pl.col("value").alias("y_true"),
        ).with_columns(
            pl.lit(view).alias("view"),
            pl.lit(method).alias("method"),
            pl.lit(cfg["split_run"]).alias("split_id"),
            pl.Series("prediction", prediction, dtype=pl.Float64),
        )
        rows.append(row)
        slope, intercept = calibrations[(view, method)]
        record = {
            "fold": fold,
            "view": view,
            "method": method,
            "train": int(train.sum()),
            "validation": int(validation.sum()),
            "available_confident_pairs": available,
            "sampled_confident_pairs": len(high),
            "calibration_slope": slope,
            "calibration_intercept": intercept,
            "fold_ST_RAE": st_rae(row, prediction),
        }
        fit_rows.append(record)
        print(json.dumps(record), flush=True)
    pl.concat(rows).write_parquet(out / "oof-progress.parquet")

oof = pl.concat(rows).sort("view", "method", "fold", "compound_id")
expected_methods = {tuple(value) for value in cfg["methods"]}
observed_methods = set(zip(oof["view"], oof["method"], strict=True))
if observed_methods != expected_methods:
    raise ValueError(f"Method inventory mismatch: {observed_methods ^ expected_methods}")
if oof.height != latest.height * len(expected_methods):
    raise ValueError("OOF row count mismatch")
oof.write_parquet(out / "oof.parquet")
pl.concat(memberships).write_parquet(out / "memberships.parquet")

baseline = latest.select(*keys, "prediction").rename({"prediction": "baseline"})
metric_rows = []
diversity_rows = []
for (view, method), frame in oof.group_by("view", "method"):
    frame = frame.sort(keys).join(baseline, on=keys, validate="1:1")
    prediction = frame["prediction"].to_numpy()
    base_prediction = frame["baseline"].to_numpy()
    y = frame["y_true"].to_numpy()
    diversity_rows.append(
        {
            "view": view,
            "method": method,
            "prediction_correlation": float(np.corrcoef(prediction, base_prediction)[0, 1]),
            "residual_correlation": float(np.corrcoef(y - prediction, y - base_prediction)[0, 1]),
        }
    )
    for weight in cfg["blend_weights"]:
        blended = (1 - weight) * base_prediction + weight * prediction
        metric_rows.append(
            {
                "view": view,
                "method": method,
                "new_weight": weight,
                "ST_RAE": st_rae(frame, blended),
                "MAE": float(np.abs(y - blended).mean()),
                "Spearman": float(spearmanr(y, blended).statistic),
            }
        )

metrics = pl.DataFrame(metric_rows).sort("ST_RAE", "view", "method")
diversity = pl.DataFrame(diversity_rows).sort("residual_correlation", "view", "method")
metrics.write_csv(out / "metrics.csv")
diversity.write_csv(out / "diversity.csv")
pl.DataFrame(fit_rows).sort("fold", "view", "method").write_csv(out / "fits.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "baseline_latest_ST_RAE": st_rae(latest, latest["prediction"].to_numpy()),
            "baseline_current_ST_RAE": st_rae(current, current["prediction"].to_numpy()),
            "pair_definition": (
                "Outer-training-only confidence-disjoint pairs are sampled proportional to "
                "interval gap and weighted by clipped normalized gap."
            ),
            "selection_note": "Broad fixed-recipe discovery screen; no automatic adoption.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics.head(30), flush=True)
print(diversity, flush=True)
