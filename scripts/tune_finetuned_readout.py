"""Outer-safe Rustuna tuning of FT embeddings; inner encoder-label reuse is explicit."""

import importlib.util
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import rustuna

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name):
    return pl.read_parquet(next(inp.glob("*-" + name + ".parquet")))


a = read("parent-assignments")
labels = read("measurements").filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))
aux = read("nested-log2fc")
embeddings = [pl.read_parquet(inp / name) for name in cfg["embeddings"]]
ecols = sorted(c for c in embeddings[0].columns if c.startswith("emb_"))
assert len(ecols) == 600
encoder_checks = []
for fold, emb in enumerate(embeddings):
    assert emb["fold"].unique().to_list() == [fold]
    assert emb["split_id"].unique().to_list() == [cfg["split_run"]]
    assert set(a["compound_id"]).issubset(set(emb["compound_id"]))
    member = pl.read_parquet(inp / cfg["encoder_memberships"][fold])
    expected = set(a.filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))["compound_id"])
    assert set(member.filter(pl.col("role") == "outer_holdout")["compound_id"]) == expected
    assert expected.isdisjoint(
        set(member.filter(pl.col("role").is_in(["training", "inner_validation"]))["compound_id"])
    )
    embeddings[fold] = emb.select("compound_id", *ecols)
    encoder_checks.append(
        dict(fold=fold, outer_compounds=len(expected), outer_labels_excluded=True)
    )
(out / "encoder-exclusion-audit.json").write_text(json.dumps(encoder_checks, indent=2))
spec = importlib.util.spec_from_file_location(
    "official", next(inp.glob("*-custom_scoring_functions.py"))
)
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)


def fit_endpoint(outer, endpoint):
    emb = embeddings[outer]
    features = a.join(emb, on="compound_id", validate="1:1").join(
        aux.filter(pl.col("outer_fold") == outer).select(
            "compound_id", *[f"aux_{i}" for i in range(4)]
        ),
        on="compound_id",
        validate="1:1",
    )
    table = (
        labels.filter(pl.col("endpoint") == endpoint)
        .join(features, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    x = table.select(ecols + [f"aux_{i}" for i in range(4)]).to_numpy().astype(np.float32)
    assert np.isfinite(x).all()
    outer_valid = (table["fold"].to_numpy() == outer) & ~table["is_anchor"].to_numpy()
    outer_train = ~outer_valid
    # Hold out whole parent groups; groups containing anchors remain inner-training.
    anchor_groups = set(a.filter(pl.col("is_anchor"))["group"])
    eligible = sorted(set(table.filter(pl.Series(outer_train))["group"]) - anchor_groups)
    rng = np.random.default_rng(cfg["seed"] + outer)
    selected_groups = set(
        rng.permutation(eligible)[: max(1, round(len(eligible) * cfg["inner_fraction"]))]
    )
    inner_valid = outer_train & np.array([g in selected_groups for g in table["group"]])
    inner_train = outer_train & ~inner_valid
    assert inner_valid.sum() >= 20 and inner_train.sum() >= 100
    assert not (inner_valid & table["is_anchor"].to_numpy()).any()
    assert set(table.filter(pl.Series(inner_train))["group"]).isdisjoint(
        set(table.filter(pl.Series(inner_valid))["group"])
    )
    # Slice inner labels/bounds once; the objective never reads outer-validation labels.
    ytrain = table.filter(pl.Series(inner_train))["value"].to_numpy()
    ival = table.filter(pl.Series(inner_valid))
    yvalid, low, high = (ival[c].to_numpy() for c in ["value", "conf_low", "conf_high"])
    run_dir = out / f"fold-{outer}" / endpoint
    run_dir.mkdir(parents=True)
    membership = table.select("compound_id", "endpoint", "fold", "group", "is_anchor").with_columns(
        pl.lit(outer).alias("outer_fold"),
        pl.lit(cfg["split_run"]).alias("split_id"),
        pl.Series(
            "role",
            np.where(
                outer_valid,
                "outer_validation",
                np.where(inner_valid, "inner_validation", "inner_train"),
            ),
        ),
    )
    membership.write_parquet(run_dir / "membership.parquet")
    records = []
    study = rustuna.create_study(
        direction="minimize",
        sampler=rustuna.samplers.TPESampler(seed=cfg["seed"] + outer, n_startup_trials=6),
    )

    def objective(
        trial,
        outer=outer,
        x=x,
        inner_train=inner_train,
        inner_valid=inner_valid,
        ytrain=ytrain,
        yvalid=yvalid,
        low=low,
        high=high,
        records=records,
        run_dir=run_dir,
    ):
        if trial.number == 0:
            params = dict(
                objective="regression_l1",
                num_leaves=15,
                min_child_samples=20,
                learning_rate=0.05,
                n_estimators=300,
                reg_lambda=1.0,
                colsample_bytree=1.0,
            )
        else:
            obj = trial.suggest_categorical("objective", ["regression_l1", "quantile"])
            params = dict(
                objective=obj,
                num_leaves=trial.suggest_categorical("num_leaves", [7, 15, 31, 63]),
                min_child_samples=trial.suggest_int("min_child_samples", 10, 80),
                learning_rate=trial.suggest_float("learning_rate", 0.02, 0.1, log=True),
                n_estimators=trial.suggest_int("n_estimators", 200, 1000),
                reg_lambda=trial.suggest_float("reg_lambda", 0.1, 20.0, log=True),
                colsample_bytree=trial.suggest_float("colsample_bytree", 0.5, 1.0),
            )
            if obj == "quantile":
                params["alpha"] = trial.suggest_float("alpha", 0.35, 0.75)
        model = lgb.LGBMRegressor(
            **params,
            n_jobs=cfg["threads"],
            random_state=cfg["seed"] + outer,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
        )
        model.fit(x[inner_train], ytrain)
        prediction = model.predict(x[inner_valid])
        score = float(
            scorer.rae_soft_threshold_absolute_error(
                yvalid, prediction, y_true_lower=low, y_true_upper=high
            )
        )
        record = dict(trial=trial.number, score=score, params=params)
        records.append(record)
        with (run_dir / "trials.jsonl").open("a") as f:
            f.write(json.dumps(record) + "\n")
        return score

    study.optimize(objective, n_trials=cfg["n_trials"])
    winner = min(records, key=lambda r: r["score"])
    final = lgb.LGBMRegressor(
        **winner["params"],
        n_jobs=cfg["threads"],
        random_state=cfg["seed"] + outer,
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
    )
    final.fit(x[outer_train], table.filter(pl.Series(outer_train))["value"].to_numpy())
    pred = final.predict(x[outer_valid])
    assert np.isfinite(pred).all()
    final.booster_.save_model(str(run_dir / "model.txt"))
    (run_dir / "selected.json").write_text(json.dumps(winner, indent=2))
    prediction_table = (
        table.filter(pl.Series(outer_valid))
        .select(
            "compound_id",
            "endpoint",
            "fold",
            "is_neighbor",
            "conf_low",
            "conf_high",
            pl.col("value").alias("y_true"),
        )
        .with_columns(
            pl.Series("prediction", pred),
            pl.lit(cfg["method"]).alias("method"),
            pl.lit(cfg["split_run"]).alias("split_id"),
        )
    )
    selection = dict(
        endpoint=endpoint,
        fold=outer,
        **winner,
        inner_train=int(inner_train.sum()),
        inner_valid=int(inner_valid.sum()),
        outer_train=int(outer_train.sum()),
        outer_valid=int(outer_valid.sum()),
    )
    print(json.dumps(selection), flush=True)
    return prediction_table, selection, membership


rows, selections, memberships = [], [], []
with ThreadPoolExecutor(max_workers=cfg["workers"]) as pool:
    futures = [
        pool.submit(fit_endpoint, outer, ep)
        for outer in range(5)
        for ep in sorted(labels["endpoint"].unique())
    ]
    for future in as_completed(futures):
        prediction_table, selection, membership = future.result()
        rows.append(prediction_table)
        selections.append(selection)
        memberships.append(membership)
        pl.concat(rows).write_parquet(out / "external-oof-progress.parquet")


result = pl.concat(rows)
assert result.height == result.select("compound_id", "endpoint", "fold").n_unique() == 6384
result.write_parquet(out / "external-oof.parquet")
pl.concat(memberships).write_parquet(out / "membership.parquet")
(out / "summary.json").write_text(
    json.dumps(
        dict(config=cfg, features=ecols + [f"aux_{i}" for i in range(4)], selections=selections),
        indent=2,
    )
)
