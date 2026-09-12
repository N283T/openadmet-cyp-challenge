"""Focused 1A2/2D6 readouts with grouped inner-only model selection."""

import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import joblib
import lightgbm as lgb
import numpy as np
import polars as pl
import rustuna
from scipy.stats import spearmanr
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from threadpoolctl import threadpool_limits

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(key):
    return pl.read_parquet(inp / cfg[key])


assignments, labels, aux = read("assignments"), read("measurements"), read("aux")
base = read("baseline").select(
    "compound_id", "endpoint", "fold", pl.col("prediction").alias("baseline")
)
medium, scores, classic = read("medium"), read("scores"), read("classic")
mcols = sorted(c for c in medium.columns if c.startswith("emb_"))
scols = sorted(c for c in scores.columns if c.startswith("score_"))
# Preserve only finite binary fingerprint columns; no global descriptor imputation.
fcols = [c for c in classic.columns if c.startswith("morgan_")]
assert len(mcols) == 600 and len(scols) == 47
if not fcols:
    fcols = [c for c in classic.columns if c.startswith("fp_")]
assert len(fcols) == 2048, classic.columns[:10]
features = (
    medium.select("compound_id", *mcols)
    .join(scores.select("compound_id", *scols), on="compound_id", validate="1:1")
    .join(classic.select("compound_id", *fcols), on="compound_id", validate="1:1")
)
threadpool_limits(limits=2)


def score(h, p):
    y, low, high = (h[c].to_numpy() for c in ["value", "conf_low", "conf_high"])
    denominator = np.mean(np.maximum(low - y.mean(), 0) + np.maximum(y.mean() - high, 0))
    assert denominator > 0
    return float(np.mean(np.maximum(low - p, 0) + np.maximum(p - high, 0)) / denominator)


def fit_one(endpoint, fold, view):
    table = (
        labels.filter(pl.col("endpoint") == endpoint)
        .join(assignments, on="compound_id", validate="1:1")
        .join(features, on="compound_id", validate="1:1")
        .join(
            aux.filter(pl.col("outer_fold") == fold).select(
                "compound_id", *[f"aux_{j}" for j in range(4)]
            ),
            on="compound_id",
            validate="1:1",
        )
        .sort("compound_id")
    )
    cols = (
        scols if view == "scores_svr" else mcols if view == "medium_svr" else mcols + scols + fcols
    ) + [f"aux_{j}" for j in range(4)]
    x = table.select(cols).to_numpy().astype(np.float32)
    assert np.isfinite(x).all()
    y = table["value"].to_numpy()
    outer = (table["fold"].to_numpy() == fold) & ~table["is_anchor"].to_numpy()
    eligible_groups = sorted(
        set(table.filter(pl.Series(~outer))["group"])
        - set(assignments.filter(pl.col("is_anchor"))["group"])
    )
    rng = np.random.default_rng(42 + fold)
    inner_groups = set(
        rng.permutation(eligible_groups)[: max(1, round(len(eligible_groups) * 0.2))]
    )
    valid = ~outer & np.array([g in inner_groups for g in table["group"]])
    train = ~outer & ~valid
    assert train.sum() > 100 and valid.sum() > 20
    assert set(table.filter(pl.Series(train))["group"]).isdisjoint(
        table.filter(pl.Series(valid))["group"]
    )
    directory = out / f"{endpoint}-fold{fold}-{view}"
    directory.mkdir()
    table.select("compound_id", "endpoint", "fold", "group", "is_anchor").with_columns(
        pl.Series(
            "role", np.where(outer, "outer", np.where(valid, "inner_validation", "inner_train"))
        )
    ).write_parquet(directory / "membership.parquet")
    records = []

    def make_model(params):
        if view.endswith("svr"):
            return make_pipeline(StandardScaler(), SVR(**params, cache_size=256))
        return lgb.LGBMRegressor(
            **params,
            random_state=42 + fold,
            n_jobs=2,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
        )

    def objective(trial):
        if view.endswith("svr"):
            params = dict(
                kernel="rbf",
                C=trial.suggest_float("C", 0.1, 100, log=True),
                epsilon=trial.suggest_float("epsilon", 0.02, 0.4, log=True),
                gamma=trial.suggest_float("gamma_factor", 0.01, 10, log=True) / len(cols),
            )
        elif trial.number == 0:
            params = dict(
                objective="regression_l1",
                n_estimators=300,
                num_leaves=15,
                min_child_samples=20,
                learning_rate=0.05,
                reg_lambda=1,
            )
        else:
            params = dict(
                objective=trial.suggest_categorical("objective", ["regression_l1", "quantile"]),
                n_estimators=trial.suggest_int("n_estimators", 200, 800),
                num_leaves=trial.suggest_categorical("num_leaves", [7, 15, 31]),
                min_child_samples=trial.suggest_int("min_child_samples", 10, 60),
                learning_rate=trial.suggest_float("learning_rate", 0.02, 0.1, log=True),
                reg_lambda=trial.suggest_float("reg_lambda", 0.1, 20, log=True),
                colsample_bytree=trial.suggest_float("colsample_bytree", 0.3, 1),
            )
            if params["objective"] == "quantile":
                params["alpha"] = trial.suggest_float("alpha", 0.3, 0.7)
        model = make_model(params)
        model.fit(x[train], y[train])
        loss = score(table.filter(pl.Series(valid)), model.predict(x[valid]))
        record = dict(trial=trial.number, params=params, inner_score=loss)
        records.append(record)
        with (directory / "trials.jsonl").open("a") as f:
            f.write(json.dumps(record) + "\n")
        return loss

    study = rustuna.create_study(
        direction="minimize",
        sampler=rustuna.samplers.TPESampler(seed=42 + fold, n_startup_trials=6),
    )
    study.optimize(objective, n_trials=cfg["trials"])
    winner = min(records, key=lambda r: r["inner_score"])
    model = make_model(winner["params"])
    model.fit(x[~outer], y[~outer])
    p = model.predict(x[outer])
    joblib.dump(dict(model=model, columns=cols), directory / "model.joblib")
    h = table.filter(pl.Series(outer))
    result = h.select(
        "compound_id",
        "endpoint",
        "fold",
        "is_neighbor",
        "conf_low",
        "conf_high",
        pl.col("value").alias("y_true"),
    ).with_columns(
        pl.Series("prediction", p),
        pl.lit(view).alias("method"),
        pl.lit(cfg["split_run"]).alias("split_id"),
    )
    summary = dict(
        endpoint=endpoint,
        fold=fold,
        method=view,
        train=int((~outer).sum()),
        valid=int(outer.sum()),
        **winner,
        outer_score=score(h, p),
    )
    (directory / "selection.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)
    return result, summary


rows, selections = [], []
with ThreadPoolExecutor(max_workers=3) as pool:
    futures = [
        pool.submit(fit_one, ep, fold, view)
        for ep in cfg["endpoints"]
        for fold in range(5)
        for view in cfg["views"]
    ]
    for future in as_completed(futures):
        result, selection = future.result()
        rows.append(result)
        selections.append(selection)
        pl.concat(rows).write_parquet(out / "oof-progress.parquet")
oof = pl.concat(rows)
assert oof.height == 2829 * len(cfg["views"])
oof.write_parquet(out / "focus-oof.parquet")
metrics = []
for (ep, method), h in oof.group_by("endpoint", "method"):
    h = h.join(base, on=["compound_id", "endpoint", "fold"], validate="1:1")
    for subset in ["all", "neighbors"]:
        g = h if subset == "all" else h.filter(pl.col("is_neighbor"))
        for blend in [0, 0.25, 0.5, 1]:
            p = (1 - blend) * g["baseline"].to_numpy() + blend * g["prediction"].to_numpy()
            metrics.append(
                dict(
                    endpoint=ep,
                    method=method,
                    subset=subset,
                    new_weight=blend,
                    n=g.height,
                    ST_RAE=score(g.rename({"y_true": "value"}), p),
                    Spearman=float(spearmanr(g["y_true"].to_numpy(), p).statistic),
                )
            )
pl.DataFrame(metrics).write_csv(out / "metrics.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            selections=selections,
            inner_aux_reuse=(
                "Auxiliaries crossfit outer folds; no additional independent inner encoder fit"
            ),
            blends="Prespecified diagnostic weights; no automatic adoption",
        ),
        indent=2,
    )
)
