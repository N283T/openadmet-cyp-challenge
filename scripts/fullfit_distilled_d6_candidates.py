"""Full-fit atomic D6 readouts and package several placement-matched candidates."""

import importlib.util
import itertools
import json
import os
from pathlib import Path

import joblib
import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from catboost import CatBoostRegressor
from scipy.stats import spearmanr
from tabicl import TabICLRegressor

from cyp_challenge.tracks import validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoint = "CYP2D6_pIC50_direct_inhibition"
keys = ["fold", "compound_id"]

measurements = pl.read_parquet(inp / "measurements.parquet")
train_labels = measurements.filter(pl.col("endpoint") == endpoint).select(
    "compound_id", "value"
)
aliases = pl.read_parquet(inp / "aliases.parquet")
test_ids = (
    aliases.filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
reference = pl.read_csv(inp / "reference.csv")
test = reference.select("SMILES", "Molecule_Name").join(
    test_ids, on="Molecule_Name", validate="1:1"
)
assert train_labels.height == cfg["expected_train_rows"]
assert test.height == test["compound_id"].n_unique() == 750
assert set(train_labels["compound_id"]).isdisjoint(test["compound_id"])


def load_feature(filename: str, prefixes: list[str]) -> tuple[pl.DataFrame, list[str]]:
    frame = pl.read_parquet(inp / filename)
    columns = sorted(
        column
        for column in frame.columns
        if any(column.startswith(prefix) for prefix in prefixes)
    )
    assert columns and frame["compound_id"].n_unique() == frame.height
    assert set(train_labels["compound_id"]).union(test["compound_id"]).issubset(
        set(frame["compound_id"])
    )
    assert np.isfinite(frame.select(columns).to_numpy()).all()
    return frame.select("compound_id", *columns), columns


def matrices(filename: str, prefixes: list[str]):
    feature, columns = load_feature(filename, prefixes)
    train = train_labels.join(feature, on="compound_id", validate="1:1").sort("compound_id")
    blind = test.join(feature, on="compound_id", validate="1:1").sort("compound_id")
    return (
        train.select(columns).to_numpy().astype(np.float32),
        train["value"].to_numpy(),
        blind.select(columns).to_numpy().astype(np.float32),
        blind,
        columns,
    )


def fit_lgbm(name: str, filename: str, prefixes: list[str]) -> pl.DataFrame:
    x, y, xt, blind, columns = matrices(filename, prefixes)
    model = lgb.LGBMRegressor(
        objective="regression_l1",
        n_estimators=300,
        num_leaves=15,
        min_child_samples=20,
        learning_rate=0.05,
        reg_lambda=1.0,
        n_jobs=cfg["num_threads"],
        random_state=42,
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
    )
    model.fit(x, y)
    prediction = model.predict(xt).astype(float)
    model.booster_.save_model(str(out / f"{name}.lgbm.txt"))
    (out / f"{name}-columns.json").write_text(json.dumps(columns, indent=2))
    return blind.select("compound_id", "Molecule_Name").with_columns(
        pl.Series(name, prediction)
    )


def fit_catboost(name: str, filename: str, prefixes: list[str]) -> pl.DataFrame:
    x, y, xt, blind, columns = matrices(filename, prefixes)
    model = CatBoostRegressor(
        iterations=500,
        depth=6,
        learning_rate=0.03,
        l2_leaf_reg=5,
        loss_function="MAE",
        random_seed=42,
        thread_count=cfg["num_threads"],
        verbose=False,
        allow_writing_files=False,
        task_type="CPU",
    )
    model.fit(x, y)
    prediction = np.asarray(model.predict(xt), dtype=float)
    model.save_model(str(out / f"{name}.catboost.cbm"))
    (out / f"{name}-columns.json").write_text(json.dumps(columns, indent=2))
    return blind.select("compound_id", "Molecule_Name").with_columns(
        pl.Series(name, prediction)
    )


def existing_prediction(name: str, filename: str) -> pl.DataFrame:
    frame = (
        pl.read_parquet(inp / filename)
        .filter(pl.col("endpoint") == endpoint)
        .select("compound_id", pl.col("prediction").alias(name))
    )
    assert frame.height == frame["compound_id"].n_unique() == 750
    return frame


components = test.select("compound_id", "Molecule_Name").sort("compound_id")
for name, filename in [
    ("jeremy_frozen_tabpfn", "jeremy-frozen-test.parquet"),
    ("jeremy_scores_log2fc_tabpfn", "jeremy-scores-test.parquet"),
    ("jeremy_finetuned_log2fc_lgbm", "jeremy-finetuned-test.parquet"),
]:
    components = components.join(
        existing_prediction(name, filename), on="compound_id", validate="1:1"
    )

chemeleon, chemeleon_columns = load_feature("chemeleon-raw.parquet", ["emb_"])
auxiliary = pl.read_parquet(inp / "full-aux.parquet").select(
    "compound_id", *[f"aux_{index}" for index in range(4)]
)
chemeleon = chemeleon.join(auxiliary, on="compound_id", validate="1:1")
chemeleon_columns += [f"aux_{index}" for index in range(4)]
cheme_train = train_labels.join(chemeleon, on="compound_id", validate="1:1").sort("compound_id")
cheme_test = test.join(chemeleon, on="compound_id", validate="1:1").sort("compound_id")
tabicl = TabICLRegressor(
    n_estimators=8,
    batch_size=1,
    device="cuda",
    n_jobs=cfg["num_threads"],
    model_path=str(inp / "tabicl-regressor.ckpt"),
    allow_auto_download=False,
    use_amp="auto",
    use_fa3="auto",
    offload_mode="auto",
    disk_offload_dir=os.environ["TMPDIR"],
    random_state=42,
)
tabicl.fit(
    cheme_train.select(chemeleon_columns).to_numpy().astype(np.float32),
    cheme_train["value"].to_numpy(),
)
chemeleon_prediction = np.asarray(
    tabicl.predict(cheme_test.select(chemeleon_columns).to_numpy().astype(np.float32)),
    dtype=float,
)
joblib.dump(tabicl, out / "chemeleon_raw2048_log2fc_tabicl.joblib")
components = components.join(
    cheme_test.select("compound_id").with_columns(
        pl.Series("chemeleon_raw2048_log2fc_tabicl", chemeleon_prediction)
    ),
    on="compound_id",
    validate="1:1",
)
del tabicl
torch.cuda.empty_cache()

new_specs = [
    ("owned_raw_scores_lgbm", "owned-raw-scores.parquet", ["raw_", "score_"], "lgbm"),
    ("anchor59_raw_scores_lgbm", "anchor59-raw-scores.parquet", ["raw_", "score_"], "lgbm"),
    ("anchor59_scores_lgbm", "anchor59-scores.parquet", ["score_"], "lgbm"),
    ("pubonly_raw_lgbm", "pubonly-raw.parquet", ["raw_"], "lgbm"),
    ("scratch_raw_catboost", "scratch-raw.parquet", ["raw_"], "catboost"),
    ("anchor59_raw_scores_catboost", "anchor59-raw-scores.parquet", ["raw_", "score_"], "catboost"),
]
for name, filename, prefixes, backend in new_specs:
    fitted = (
        fit_lgbm(name, filename, prefixes)
        if backend == "lgbm"
        else fit_catboost(name, filename, prefixes)
    )
    components = components.join(
        fitted.select("compound_id", name), on="compound_id", validate="1:1"
    )
components = components.sort("compound_id")

component_names = cfg["component_order"]
assert components.columns == ["compound_id", "Molecule_Name", *component_names]
assert np.isfinite(components.select(component_names).to_numpy()).all()
components.write_parquet(out / "atomic-test-predictions.parquet")

spec = importlib.util.spec_from_file_location("official", inp / "custom_scoring_functions.py")
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
oof_frames = []
for item in cfg["oof_components"]:
    frame = (
        pl.read_parquet(inp / item["file"])
        .filter(pl.col("endpoint") == endpoint)
        .sort(keys)
    )
    assert frame.select(keys).n_unique() == frame.height
    oof_frames.append(frame)
base_oof = oof_frames[0]
for frame in oof_frames[1:]:
    for column in [*keys, "y_true", "conf_low", "conf_high"]:
        assert frame[column].equals(base_oof[column])
oof_matrix = np.column_stack([frame["prediction"].to_numpy() for frame in oof_frames])
y_true = base_oof["y_true"].to_numpy()
low = base_oof["conf_low"].to_numpy()
high = base_oof["conf_high"].to_numpy()


def st_rae(prediction: np.ndarray) -> float:
    return float(
        scorer.rae_soft_threshold_absolute_error(
            y_true, prediction, y_true_lower=low, y_true_upper=high
        )
    )


subset_scores = []
for size in range(1, len(component_names) + 1):
    for subset in itertools.combinations(range(len(component_names)), size):
        subset_scores.append((st_rae(oof_matrix[:, subset].mean(axis=1)), subset))
best_score, best_subset = min(subset_scores, key=lambda item: (item[0], len(item[1]), item[1]))
best_names = [component_names[index] for index in best_subset]
assert best_names == cfg["recipes"]["stable5"]

base_candidate = pl.read_csv(inp / "base-candidate.csv").sort("Molecule_Name")
assert base_candidate.height == 750
base_columns = base_candidate.columns
old_by_component = components.select("Molecule_Name").join(
    base_candidate.select("Molecule_Name", endpoint),
    on="Molecule_Name",
    validate="1:1",
    maintain_order="left",
)
old_d6 = old_by_component[endpoint].to_numpy()
recipe_rows = []
for recipe_name, members in cfg["recipes"].items():
    raw = components.select(members).mean_horizontal().to_numpy()
    placed = (raw - raw.mean()) / raw.std() * old_d6.std() + old_d6.mean()
    replacement = components.select("Molecule_Name").with_columns(
        pl.Series(endpoint, placed)
    )
    candidate = (
        base_candidate.drop(endpoint)
        .join(
            replacement,
            on="Molecule_Name",
            validate="1:1",
            maintain_order="left",
        )
        .select(base_columns)
    )
    for other_endpoint in [
        "CYP1A2_pIC50_direct_inhibition",
        "CYP2C9_pIC50_direct_inhibition",
        "CYP3A4_pIC50_direct_inhibition",
    ]:
        assert candidate[other_endpoint].equals(base_candidate[other_endpoint])
    aligned_placed = components.select("Molecule_Name").join(
        candidate.select("Molecule_Name", endpoint),
        on="Molecule_Name",
        validate="1:1",
        maintain_order="left",
    )[endpoint].to_numpy()
    assert np.array_equal(aligned_placed, placed)
    parquet = out / f"candidate-{recipe_name}.parquet"
    csv = out / f"candidate-{recipe_name}.csv"
    candidate.write_parquet(parquet)
    candidate.write_csv(csv)
    validation = validate(parquet, "regression", inp / "reference.csv")
    if not validation["valid"] or not validate(csv, "regression", inp / "reference.csv")["valid"]:
        raise ValueError(validation)
    indices = [component_names.index(member) for member in members]
    oof_score = st_rae(oof_matrix[:, indices].mean(axis=1))
    recipe_rows.append(
        {
            "recipe": recipe_name,
            "members": "+".join(members),
            "member_count": len(members),
            "oof_ST_RAE": oof_score,
            "test_mean": float(placed.mean()),
            "test_std": float(placed.std()),
            "base_spearman": float(spearmanr(old_d6, placed).statistic),
            "mean_absolute_change": float(np.abs(placed - old_d6).mean()),
            "maximum_absolute_change": float(np.abs(placed - old_d6).max()),
            "changed_ranks": int(
                np.sum(
                    pl.Series(old_d6).rank("average").to_numpy()
                    != pl.Series(placed).rank("average").to_numpy()
                )
            ),
        }
    )

pl.DataFrame(recipe_rows).write_csv(out / "candidate-comparison.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "full_oof_best_subset": best_names,
            "full_oof_best_ST_RAE": best_score,
            "placement": {
                "source": cfg["base_candidate"],
                "mean": float(old_d6.mean()),
                "std": float(old_d6.std()),
            },
            "recipes": recipe_rows,
            "all_train_labels_used": True,
            "test_labels_used": False,
            "other_endpoints_bitwise_unchanged": True,
        },
        indent=2,
    )
    + "\n"
)
print(json.dumps(recipe_rows), flush=True)
