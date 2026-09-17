"""Nested grouped LightGBM selection for atomic raw embedding views."""

import importlib.util
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    exact = inp / name
    paths = (
        [exact]
        if exact.exists()
        else [path for path in inp.iterdir() if path.name.partition("-")[2] == name]
    )
    if len(paths) != 1:
        raise ValueError(f"Expected one input for {name}: {paths}")
    return pl.read_parquet(paths[0])


spec = importlib.util.spec_from_file_location("official", inp / "custom_scoring_functions.py")
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)

assignments = read("parent-assignments.parquet")
measurements = read("measurements.parquet").filter(
    pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
)
anchor_groups = set(assignments.filter(pl.col("is_anchor"))["group"])


def score(rows: pl.DataFrame, prediction: np.ndarray) -> float:
    return float(
        scorer.rae_soft_threshold_absolute_error(
            rows["value"].to_numpy(),
            prediction,
            y_true_lower=rows["conf_low"].to_numpy(),
            y_true_upper=rows["conf_high"].to_numpy(),
        )
    )


def model(params: dict, random_state: int) -> lgb.LGBMRegressor:
    resolved = dict(
        objective="regression_l1",
        n_estimators=300,
        num_leaves=15,
        min_child_samples=20,
        learning_rate=0.05,
        reg_lambda=1.0,
        n_jobs=cfg["num_threads"],
        random_state=random_state,
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
    )
    resolved.update(params)
    return lgb.LGBMRegressor(**resolved)


oof_rows = []
selection_rows = []
inner_rows = []
for candidate in cfg["candidates"]:
    embeddings = read(candidate["file"])
    columns = sorted(column for column in embeddings.columns if column.startswith("emb_"))
    if len(columns) != 600 or embeddings["compound_id"].n_unique() != embeddings.height:
        raise ValueError(f"Invalid raw600 embedding: {candidate['name']}")
    for endpoint in candidate["endpoints"]:
        table = (
            measurements.filter(pl.col("endpoint") == endpoint)
            .join(assignments, on="compound_id", validate="1:1")
            .join(embeddings, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        x = table.select(columns).to_numpy().astype(np.float32)
        y = table["value"].to_numpy()
        groups = table["group"].to_numpy()
        folds = table["fold"].to_numpy()
        anchors = table["is_anchor"].to_numpy()
        if not np.isfinite(x).all():
            raise ValueError(f"Non-finite embedding: {candidate['name']}")
        for outer_fold in range(5):
            outer_valid = (folds == outer_fold) & ~anchors
            outer_train = ~outer_valid
            eligible_groups = sorted(set(groups[outer_train]) - anchor_groups)
            rng = np.random.default_rng(cfg["seed"] + outer_fold)
            shuffled_groups = rng.permutation(eligible_groups)
            inner_group_fold = {
                group: index % cfg["inner_folds"]
                for index, group in enumerate(shuffled_groups)
            }
            trial_records = []
            for grid_index, grid in enumerate(cfg["lgbm_grid"]):
                scores = []
                for inner_fold in range(cfg["inner_folds"]):
                    inner_valid = outer_train & np.array(
                        [inner_group_fold.get(group) == inner_fold for group in groups]
                    )
                    inner_train = outer_train & ~inner_valid
                    if anchors[inner_valid].any():
                        raise AssertionError("Anchor entered inner validation")
                    fitted = model(grid["params"], cfg["seed"] + outer_fold)
                    fitted.fit(x[inner_train], y[inner_train])
                    inner_table = table.filter(pl.Series(inner_valid))
                    inner_score = score(inner_table, fitted.predict(x[inner_valid]))
                    scores.append(inner_score)
                    inner_rows.append(
                        {
                            "method": candidate["name"],
                            "endpoint": endpoint,
                            "outer_fold": outer_fold,
                            "inner_fold": inner_fold,
                            "grid": grid["name"],
                            "ST_RAE": inner_score,
                            "train_n": int(inner_train.sum()),
                            "validation_n": int(inner_valid.sum()),
                        }
                    )
                trial_records.append((float(np.mean(scores)), grid_index, grid))
            mean_score, _, winner = min(trial_records)
            final_model = model(winner["params"], cfg["seed"] + outer_fold)
            final_model.fit(x[outer_train], y[outer_train])
            prediction = final_model.predict(x[outer_valid]).astype(float)
            output_dir = out / candidate["name"] / endpoint / f"fold-{outer_fold}"
            output_dir.mkdir(parents=True)
            final_model.booster_.save_model(str(output_dir / "model.txt"))
            outer_table = table.filter(pl.Series(outer_valid))
            oof_rows.append(
                outer_table.select(
                    "compound_id",
                    "endpoint",
                    "fold",
                    "is_neighbor",
                    "conf_low",
                    "conf_high",
                    pl.col("value").alias("y_true"),
                ).with_columns(
                    pl.lit(candidate["name"]).alias("method"),
                    pl.lit(cfg["split_run"]).alias("split_id"),
                    pl.Series("prediction", prediction),
                )
            )
            selection = {
                "method": candidate["name"],
                "endpoint": endpoint,
                "outer_fold": outer_fold,
                "grid": winner["name"],
                "inner_mean_ST_RAE": mean_score,
                "outer_ST_RAE": score(outer_table, prediction),
            }
            selection_rows.append(selection)
            pl.concat(oof_rows).write_parquet(out / "external-oof-progress.parquet")
            pl.DataFrame(selection_rows).write_csv(out / "selections-progress.csv")
            pl.DataFrame(inner_rows).write_csv(out / "inner-scores-progress.csv")
            print(json.dumps(selection), flush=True)

oof = pl.concat(oof_rows)
oof.write_parquet(out / "external-oof.parquet")
pl.DataFrame(selection_rows).write_csv(out / "selections.csv")
pl.DataFrame(inner_rows).write_csv(out / "inner-scores.csv")

metrics = []
for candidate in cfg["candidates"]:
    for endpoint in candidate["endpoints"]:
        rows = oof.filter(
            (pl.col("method") == candidate["name"]) & (pl.col("endpoint") == endpoint)
        )
        metrics.append(
            {
                "method": candidate["name"],
                "endpoint": endpoint,
                "ST_RAE": score(
                    rows.rename({"y_true": "value"}), rows["prediction"].to_numpy()
                ),
                "n": rows.height,
            }
        )
pl.DataFrame(metrics).write_csv(out / "metrics.csv")
(out / "summary.json").write_text(
    json.dumps({"config": cfg, "metrics": metrics, "selections": selection_rows}, indent=2)
    + "\n"
)
