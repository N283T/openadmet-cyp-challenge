"""Broad fixed-recipe CYP2D6 readout screen for complementary ranking signal."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from catboost import CatBoostRegressor
from scipy.stats import spearmanr
from sklearn.base import clone
from sklearn.compose import TransformedTargetRegressor
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.kernel_approximation import Nystroem
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from threadpoolctl import threadpool_limits
from xgboost import XGBRegressor

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


class TriangularKANLayer(torch.nn.Module):
    """KAN layer with learned piecewise-linear univariate edge functions."""

    def __init__(self, in_features: int, out_features: int, grid_size: int, bound: float):
        super().__init__()
        centers = torch.linspace(-bound, bound, grid_size)
        self.register_buffer("centers", centers)
        self.spacing = float(centers[1] - centers[0])
        self.base = torch.nn.Linear(in_features, out_features)
        self.spline = torch.nn.Parameter(torch.empty(out_features, in_features, grid_size))
        torch.nn.init.normal_(self.spline, mean=0.0, std=0.02 / np.sqrt(in_features))

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        clipped = values.clamp(float(self.centers[0]), float(self.centers[-1]))
        basis = (1 - (clipped.unsqueeze(-1) - self.centers).abs() / self.spacing).clamp_min(0)
        return self.base(torch.nn.functional.silu(values)) + torch.einsum(
            "bik,oik->bo", basis, self.spline
        )


class CompactKAN(torch.nn.Module):
    def __init__(self, dimensions: int, hidden: int, grid_size: int, bound: float):
        super().__init__()
        self.first = TriangularKANLayer(dimensions, hidden, grid_size, bound)
        self.second = TriangularKANLayer(hidden, 1, grid_size, bound)

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return self.second(self.first(values)).squeeze(-1)


def fit_kan(
    x_train: np.ndarray, y_train: np.ndarray, x_valid: np.ndarray, seed: int, directory: Path
) -> np.ndarray:
    settings = cfg["kan"]
    scaler = StandardScaler().fit(x_train)
    x_train_scaled = scaler.transform(x_train)
    x_valid_scaled = scaler.transform(x_valid)
    dimensions = min(settings["pca_dimensions"], x_train.shape[0] - 1, x_train.shape[1])
    pca = PCA(n_components=dimensions, svd_solver="randomized", random_state=seed).fit(
        x_train_scaled
    )
    x_train_pca = pca.transform(x_train_scaled).astype(np.float32)
    x_valid_pca = pca.transform(x_valid_scaled).astype(np.float32)
    latent_scaler = StandardScaler().fit(x_train_pca)
    x_train_pca = np.clip(
        latent_scaler.transform(x_train_pca), -settings["bound"], settings["bound"]
    ).astype(np.float32)
    x_valid_pca = np.clip(
        latent_scaler.transform(x_valid_pca), -settings["bound"], settings["bound"]
    ).astype(np.float32)
    y_mean, y_scale = float(y_train.mean()), float(y_train.std())
    if y_scale <= 0:
        raise ValueError("Constant KAN training target")
    target = ((y_train - y_mean) / y_scale).astype(np.float32)

    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = CompactKAN(dimensions, settings["hidden"], settings["grid_size"], settings["bound"]).to(
        device
    )
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=settings["learning_rate"], weight_decay=settings["weight_decay"]
    )
    generator = torch.Generator().manual_seed(seed)
    train_tensor = torch.from_numpy(x_train_pca)
    target_tensor = torch.from_numpy(target)
    dataset = torch.utils.data.TensorDataset(train_tensor, target_tensor)
    loader = torch.utils.data.DataLoader(
        dataset,
        batch_size=settings["batch_size"],
        shuffle=True,
        generator=generator,
    )
    model.train()
    for _epoch in range(settings["epochs"]):
        for batch_x, batch_y in loader:
            batch_x, batch_y = batch_x.to(device), batch_y.to(device)
            optimizer.zero_grad(set_to_none=True)
            prediction = model(batch_x)
            loss = torch.nn.functional.smooth_l1_loss(prediction, batch_y, beta=0.5)
            loss.backward()
            optimizer.step()
    model.eval()
    with torch.no_grad():
        prediction = (
            model(torch.from_numpy(x_valid_pca).to(device)).cpu().numpy() * y_scale + y_mean
        )
    torch.save(
        {
            "state_dict": model.cpu().state_dict(),
            "dimensions": dimensions,
            "hidden": settings["hidden"],
            "grid_size": settings["grid_size"],
            "bound": settings["bound"],
            "y_mean": y_mean,
            "y_scale": y_scale,
        },
        directory / "model.pt",
    )
    joblib.dump((scaler, pca, latent_scaler), directory / "preprocessing.joblib")
    return prediction.astype(float)


def make_models(seed: int, dimensions: int) -> dict[str, object]:
    pca_dimensions = min(cfg["pca_dimensions"], dimensions)
    scaled_pca = [
        StandardScaler(),
        PCA(n_components=pca_dimensions, svd_solver="randomized", random_state=seed),
        StandardScaler(),
    ]
    return {
        "ridge": make_pipeline(StandardScaler(), Ridge(alpha=10.0)),
        "elastic_net": make_pipeline(
            StandardScaler(), ElasticNet(alpha=0.01, l1_ratio=0.1, max_iter=20000)
        ),
        "pls16": make_pipeline(
            StandardScaler(), PLSRegression(n_components=min(16, dimensions), scale=False)
        ),
        "extra_trees": ExtraTreesRegressor(
            n_estimators=cfg["tree_estimators"],
            min_samples_leaf=3,
            max_features=0.5,
            n_jobs=cfg["threads"],
            random_state=seed,
        ),
        "random_forest": RandomForestRegressor(
            n_estimators=cfg["tree_estimators"],
            min_samples_leaf=3,
            max_features=0.5,
            n_jobs=cfg["threads"],
            random_state=seed,
        ),
        "rbf_svr": make_pipeline(StandardScaler(), SVR(C=10.0, epsilon=0.1, gamma="scale")),
        "nystroem_ridge": make_pipeline(
            StandardScaler(),
            Nystroem(
                kernel="rbf",
                gamma=1.0 / dimensions,
                n_components=min(256, dimensions * 4),
                random_state=seed,
            ),
            Ridge(alpha=1.0),
        ),
        "knn": make_pipeline(
            *scaled_pca,
            KNeighborsRegressor(n_neighbors=15, weights="distance", p=2),
        ),
        "mlp": TransformedTargetRegressor(
            regressor=make_pipeline(
                *scaled_pca,
                MLPRegressor(
                    hidden_layer_sizes=(64, 32),
                    activation="relu",
                    alpha=0.01,
                    learning_rate_init=0.001,
                    max_iter=500,
                    early_stopping=True,
                    validation_fraction=0.15,
                    n_iter_no_change=30,
                    random_state=seed,
                ),
            ),
            transformer=StandardScaler(),
        ),
        "xgboost": XGBRegressor(
            n_estimators=500,
            max_depth=3,
            learning_rate=0.03,
            min_child_weight=5,
            subsample=0.8,
            colsample_bytree=0.5,
            reg_lambda=5.0,
            objective="reg:squarederror",
            tree_method="hist",
            n_jobs=cfg["threads"],
            random_state=seed,
        ),
        "catboost": CatBoostRegressor(
            iterations=500,
            depth=6,
            learning_rate=0.03,
            l2_leaf_reg=5.0,
            loss_function="MAE",
            random_seed=seed,
            thread_count=cfg["threads"],
            verbose=False,
            allow_writing_files=False,
        ),
    }


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
if latest.height != current.height or latest.height != base.filter(~pl.col("is_anchor")).height:
    raise ValueError("Baseline coverage mismatch")
for column in keys + ["y_true", "conf_low", "conf_high"]:
    if not latest[column].equals(current[column]):
        raise ValueError(f"Current/latest mismatch in {column}")

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
    y = data["value"].to_numpy().astype(float)
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
        models = make_models(42 + fold, x.shape[1])
        for method in cfg["methods"]:
            directory = out / f"fold{fold}-{view}-{method}"
            directory.mkdir()
            if method == "compact_kan":
                prediction = fit_kan(x[train], y[train], x[validation], 42 + fold, directory)
            else:
                model = clone(models[method])
                model.fit(x[train], y[train])
                prediction = np.asarray(model.predict(x[validation])).reshape(-1).astype(float)
                joblib.dump({"model": model, "columns": columns}, directory / "model.joblib")
            if not np.isfinite(prediction).all():
                raise ValueError(f"Non-finite predictions for {view}/{method}/fold{fold}")
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
                pl.Series("prediction", prediction),
            )
            rows.append(row)
            record = {
                "fold": fold,
                "view": view,
                "method": method,
                "train": int(train.sum()),
                "validation": int(validation.sum()),
                "dimensions": x.shape[1],
                "fold_ST_RAE": st_rae(row, prediction),
            }
            fit_records.append(record)
            (directory / "fit.json").write_text(json.dumps(record, indent=2) + "\n")
            print(json.dumps(record), flush=True)
            pl.concat(rows).write_parquet(out / "oof-progress.parquet")

oof = pl.concat(rows).sort("view", "method", "fold", "compound_id")
expected_rows = latest.height * len(view_columns) * len(cfg["methods"])
if oof.height != expected_rows:
    raise ValueError(f"Expected {expected_rows} OOF rows, found {oof.height}")
oof.write_parquet(out / "oof.parquet")
pl.concat(memberships).write_parquet(out / "memberships.parquet")

baseline = latest.select(*keys, "prediction").rename({"prediction": "latest_prediction"})
current_predictions = current.select(*keys, "prediction").rename(
    {"prediction": "current_prediction"}
)
metric_rows = []
diversity_rows = []
for (view, method), frame in oof.group_by("view", "method"):
    frame = (
        frame.sort(keys)
        .join(baseline, on=keys, validate="1:1")
        .join(current_predictions, on=keys, validate="1:1")
    )
    prediction = frame["prediction"].to_numpy()
    latest_prediction = frame["latest_prediction"].to_numpy()
    current_prediction = frame["current_prediction"].to_numpy()
    y = frame["y_true"].to_numpy()
    matched = (prediction - prediction.mean()) / prediction.std() * latest_prediction.std()
    matched += latest_prediction.mean()
    diversity_rows.append(
        {
            "view": view,
            "method": method,
            "prediction_correlation_latest": float(
                np.corrcoef(prediction, latest_prediction)[0, 1]
            ),
            "residual_correlation_latest": float(
                np.corrcoef(y - prediction, y - latest_prediction)[0, 1]
            ),
            "prediction_correlation_current": float(
                np.corrcoef(prediction, current_prediction)[0, 1]
            ),
            "residual_correlation_current": float(
                np.corrcoef(y - prediction, y - current_prediction)[0, 1]
            ),
        }
    )
    for placement, candidate in [("raw", prediction), ("match_latest_moments", matched)]:
        for weight in cfg["blend_weights"]:
            blended = (1 - weight) * latest_prediction + weight * candidate
            metric_rows.append(
                {
                    "view": view,
                    "method": method,
                    "placement": placement,
                    "new_weight": weight,
                    "ST_RAE": st_rae(frame, blended),
                    "MAE": float(np.abs(y - blended).mean()),
                    "Spearman": float(spearmanr(y, blended).statistic),
                }
            )

metrics = pl.DataFrame(metric_rows).sort("ST_RAE", "view", "method")
diversity = pl.DataFrame(diversity_rows).sort("residual_correlation_latest", "view", "method")
metrics.write_csv(out / "metrics.csv")
diversity.write_csv(out / "diversity.csv")
pl.DataFrame(fit_records).sort("fold", "view", "method").write_csv(out / "fits.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "feature_dimensions": {
                view: len(columns) + 4 for view, columns in view_columns.items()
            },
            "baseline_latest_ST_RAE": st_rae(latest, latest["prediction"].to_numpy()),
            "baseline_current_ST_RAE": st_rae(current, current["prediction"].to_numpy()),
            "kan_definition": (
                "Two learned triangular-spline KAN layers on outer-training-fitted PCA64; "
                "compact exploratory implementation, not the pykan package."
            ),
            "selection_note": (
                "Fixed recipes and blend weights form a discovery screen. Best observed rows "
                "require a separate OOF-safe confirmation before adoption."
            ),
        },
        indent=2,
    )
    + "\n"
)
print(metrics.head(30), flush=True)
print(diversity.head(20), flush=True)
