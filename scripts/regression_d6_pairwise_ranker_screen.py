"""CYP2D6 pairwise-ranking screen using interval-confident molecular pairs."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from catboost import CatBoostRanker
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from threadpoolctl import threadpool_limits
from xgboost import XGBRanker

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
    y = frame["y_true"].to_numpy()
    lower = frame["conf_low"].to_numpy()
    upper = frame["conf_high"].to_numpy()
    error = np.maximum(lower - prediction, 0) + np.maximum(prediction - upper, 0)
    center = np.full_like(y, y.mean())
    denominator = (np.maximum(lower - center, 0) + np.maximum(center - upper, 0)).sum()
    if denominator <= 0:
        raise ValueError("Non-positive ST-RAE denominator")
    return float(error.sum() / denominator)


def confident_pairs(
    lower: np.ndarray, upper: np.ndarray, maximum: int, seed: int
) -> tuple[np.ndarray, np.ndarray, int]:
    first, second = np.triu_indices(len(lower), k=1)
    first_higher = lower[first] > upper[second]
    second_higher = lower[second] > upper[first]
    high = np.concatenate([first[first_higher], second[second_higher]])
    low = np.concatenate([second[first_higher], first[second_higher]])
    total = len(high)
    if total > maximum:
        chosen = np.random.default_rng(seed).choice(total, size=maximum, replace=False)
        high, low = high[chosen], low[chosen]
    return high, low, total


def calibrate(
    train_score: np.ndarray, train_y: np.ndarray, score: np.ndarray
) -> tuple[np.ndarray, dict]:
    variance = float(np.var(train_score))
    if variance <= 0:
        raise ValueError("Constant rank utility")
    slope = float(np.cov(train_score, train_y, ddof=0)[0, 1] / variance)
    if slope <= 0:
        raise ValueError(f"Non-positive rank calibration slope: {slope}")
    intercept = float(train_y.mean() - slope * train_score.mean())
    return intercept + slope * score, {"slope": slope, "intercept": intercept}


class RankNet(torch.nn.Module):
    def __init__(self, dimensions: int):
        super().__init__()
        self.network = torch.nn.Sequential(
            torch.nn.Linear(dimensions, 64),
            torch.nn.SiLU(),
            torch.nn.Dropout(0.1),
            torch.nn.Linear(64, 32),
            torch.nn.SiLU(),
            torch.nn.Linear(32, 1),
        )

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return self.network(values).squeeze(-1)


def fit_ranknet(
    features: np.ndarray,
    train_indices: np.ndarray,
    high: np.ndarray,
    low: np.ndarray,
    seed: int,
    path: Path,
) -> np.ndarray:
    settings = cfg["ranknet"]
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = RankNet(features.shape[1]).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=settings["learning_rate"], weight_decay=settings["weight_decay"]
    )
    tensor = torch.from_numpy(features).to(device)
    high_global = train_indices[high]
    low_global = train_indices[low]
    generator = np.random.default_rng(seed)
    model.train()
    for _epoch in range(settings["epochs"]):
        order = generator.permutation(len(high_global))
        for start in range(0, len(order), settings["batch_size"]):
            batch = order[start : start + settings["batch_size"]]
            high_index = torch.from_numpy(high_global[batch]).to(device)
            low_index = torch.from_numpy(low_global[batch]).to(device)
            optimizer.zero_grad(set_to_none=True)
            difference = model(tensor[high_index]) - model(tensor[low_index])
            loss = torch.nn.functional.softplus(-difference).mean()
            loss.backward()
            optimizer.step()
    model.eval()
    with torch.no_grad():
        utility = model(tensor).cpu().numpy().astype(float)
    torch.save(model.cpu().state_dict(), path)
    return utility


assignments, labels, medium, scores, classic, nested, latest, current = [
    read(index) for index in range(8)
]
labels = labels.filter(pl.col("endpoint") == endpoint)
medium_columns = sorted(column for column in medium.columns if column.startswith("emb_"))
score_columns = sorted(column for column in scores.columns if column.startswith("score_"))
fingerprint_columns = sorted(column for column in classic.columns if column.startswith("morgan_"))
if not fingerprint_columns:
    fingerprint_columns = sorted(column for column in classic.columns if column.startswith("fp_"))
if (len(medium_columns), len(score_columns), len(fingerprint_columns)) != (600, 47, 2048):
    raise ValueError("Unexpected feature dimensions")
static = (
    medium.select("compound_id", *medium_columns)
    .join(scores.select("compound_id", *score_columns), on="compound_id", validate="1:1")
    .join(classic.select("compound_id", *fingerprint_columns), on="compound_id", validate="1:1")
)
base = labels.join(assignments, on="compound_id", validate="1:1").sort("compound_id")
latest = latest.filter(pl.col("endpoint") == endpoint).sort(keys)
current = current.filter((pl.col("endpoint") == endpoint) & (pl.col("component_count") == 2)).sort(
    keys
)
for column in keys + ["y_true", "conf_low", "conf_high"]:
    if not latest[column].equals(current[column]):
        raise ValueError(f"Baseline mismatch in {column}")

view_columns = {
    "scores51": score_columns,
    "medium604": medium_columns,
    "hybrid2699": medium_columns + score_columns + fingerprint_columns,
}
rows = []
fit_records = []
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
    train_indices = np.flatnonzero(train)
    y = data["value"].to_numpy().astype(float)
    lower = data["conf_low"].to_numpy().astype(float)
    upper = data["conf_high"].to_numpy().astype(float)
    high, low, available_pairs = confident_pairs(
        lower[train], upper[train], cfg["maximum_pairs"], 42 + fold
    )
    if len(high) < 1000:
        raise ValueError(f"Too few confident pairs in fold {fold}: {len(high)}")
    memberships.append(
        data.select("compound_id", "endpoint", "fold", "group", "is_anchor").with_columns(
            pl.lit(fold).alias("outer_fold"), pl.Series("validation", validation)
        )
    )
    for view, columns in view_columns.items():
        columns = columns + [f"aux_{index}" for index in range(4)]
        x = data.select(columns).to_numpy().astype(np.float32)
        if not np.isfinite(x).all():
            raise ValueError(f"Non-finite features in {view}, fold {fold}")
        directory = out / f"fold{fold}-{view}"
        directory.mkdir()
        scaler = StandardScaler().fit(x[train])
        scaled = scaler.transform(x)
        dimensions = min(cfg["pca_dimensions"], int(train.sum()) - 1, x.shape[1])
        pca = PCA(n_components=dimensions, svd_solver="randomized", random_state=42 + fold).fit(
            scaled[train]
        )
        projected = pca.transform(scaled)
        latent_scaler = StandardScaler().fit(projected[train])
        compact = latent_scaler.transform(projected).astype(np.float32)
        joblib.dump((scaler, pca, latent_scaler), directory / "pair-preprocessing.joblib")

        difference = compact[train_indices[high]] - compact[train_indices[low]]
        pair_x = np.vstack([difference, -difference])
        pair_y = np.concatenate([np.ones(len(difference)), np.zeros(len(difference))])
        rank_models = {
            "pairwise_logistic": LogisticRegression(C=1.0, max_iter=3000, random_state=42 + fold),
            "rank_svm": LinearSVC(C=0.1, max_iter=10000, random_state=42 + fold),
        }
        utilities = {}
        for method, model in rank_models.items():
            model.fit(pair_x, pair_y)
            utilities[method] = compact @ model.coef_.reshape(-1)
            joblib.dump(model, directory / f"{method}.joblib")
        utilities["ranknet"] = fit_ranknet(
            compact,
            train_indices,
            high,
            low,
            42 + fold,
            directory / "ranknet.pt",
        )

        xgb = XGBRanker(
            objective="rank:pairwise",
            n_estimators=500,
            max_depth=3,
            learning_rate=0.03,
            min_child_weight=5,
            subsample=0.8,
            colsample_bytree=0.5,
            reg_lambda=5.0,
            tree_method="hist",
            n_jobs=cfg["threads"],
            random_state=42 + fold,
        )
        xgb.fit(x[train], y[train], qid=np.zeros(int(train.sum()), dtype=np.int32))
        utilities["xgb_pairwise"] = xgb.predict(x).astype(float)
        xgb.save_model(directory / "xgb-pairwise.json")

        cat = CatBoostRanker(
            iterations=500,
            depth=6,
            learning_rate=0.03,
            l2_leaf_reg=5.0,
            loss_function="PairLogit",
            random_seed=42 + fold,
            thread_count=cfg["threads"],
            verbose=False,
            allow_writing_files=False,
        )
        cat.fit(
            x[train],
            y[train],
            group_id=np.zeros(int(train.sum()), dtype=np.int32),
            pairs=np.column_stack([high, low]),
        )
        utilities["cat_pairlogit"] = cat.predict(x).astype(float)
        cat.save_model(directory / "cat-pairlogit.cbm")

        for method, utility in utilities.items():
            prediction, calibration = calibrate(utility[train], y[train], utility[validation])
            held = data.filter(pl.Series(validation))
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
            record = {
                "fold": fold,
                "view": view,
                "method": method,
                "train": int(train.sum()),
                "validation": int(validation.sum()),
                "dimensions": x.shape[1],
                "compact_dimensions": compact.shape[1]
                if method in cfg["compact_methods"]
                else None,
                "available_confident_pairs": available_pairs,
                "sampled_confident_pairs": len(high),
                "calibration_slope": calibration["slope"],
                "calibration_intercept": calibration["intercept"],
                "fold_ST_RAE": st_rae(row, prediction),
            }
            fit_records.append(record)
            print(json.dumps(record), flush=True)
            pl.concat(rows).write_parquet(out / "oof-progress.parquet")

oof = pl.concat(rows).sort("view", "method", "fold", "compound_id")
expected = latest.height * len(view_columns) * len(cfg["methods"])
if oof.height != expected:
    raise ValueError(f"Expected {expected} OOF rows, found {oof.height}")
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
pl.DataFrame(fit_records).sort("fold", "view", "method").write_csv(out / "fits.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "baseline_latest_ST_RAE": st_rae(latest, latest["prediction"].to_numpy()),
            "baseline_current_ST_RAE": st_rae(current, current["prediction"].to_numpy()),
            "pair_definition": (
                "Logistic, RankSVM, RankNet, and CatBoost training-only pairs require "
                "conf_low(high) > conf_high(low); XGBoost rank:pairwise orders the full "
                "outer-training group."
            ),
            "calibration": "Positive affine utility-to-pIC50 fit on outer-training labels only.",
            "selection_note": "Broad fixed-recipe discovery screen; no automatic adoption.",
        },
        indent=2,
    )
    + "\n"
)
print(metrics.head(30), flush=True)
print(diversity, flush=True)
