"""Select external pretraining objectives by target-specific raw LGBM transfer."""

import gc
import importlib.util
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import torch
import train_external_medium as medium

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
root_config_path = Path(os.environ["CYP_CONFIG_PATH"])
cfg = json.loads(root_config_path.read_text())
medium.trainer.AssayModel = medium.MediumAssayModel


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / name)


spec = importlib.util.spec_from_file_location(
    "official", inp / "custom_scoring_functions.py"
)
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
assignments = read("parent-assignments.parquet")
labels = read("measurements.parquet").filter(
    pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
)
endpoints = sorted(labels["endpoint"].unique())


def score_embedding(path: Path, fold: int, selected_endpoints: list[str]) -> list[dict]:
    emb = pl.read_parquet(path)
    columns = sorted(column for column in emb.columns if column.startswith("emb_"))
    results = []
    for endpoint in selected_endpoints:
        table = (
            labels.filter(pl.col("endpoint") == endpoint)
            .join(assignments, on="compound_id", validate="1:1")
            .join(emb, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        x = table.select(columns).to_numpy()
        y = table["value"].to_numpy()
        valid = (table["fold"].to_numpy() == fold) & ~table["is_anchor"].to_numpy()
        model = lgb.LGBMRegressor(
            objective="regression_l1",
            n_estimators=300,
            num_leaves=15,
            min_child_samples=20,
            learning_rate=0.05,
            reg_lambda=1.0,
            n_jobs=4,
            random_state=42 + fold,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
        )
        model.fit(x[~valid], y[~valid])
        prediction = model.predict(x[valid])
        selected = table.filter(pl.Series(valid))
        value = float(
            scorer.rae_soft_threshold_absolute_error(
                selected["value"].to_numpy(),
                prediction,
                y_true_lower=selected["conf_low"].to_numpy(),
                y_true_upper=selected["conf_high"].to_numpy(),
            )
        )
        results.append(
            {"endpoint": endpoint, "fold": fold, "ST_RAE": value, "n": selected.height}
        )
    return results


def reference_scores(name: str, fold: int) -> list[dict]:
    frame = read(name)
    rows = []
    for endpoint in endpoints:
        table = frame.filter(
            (pl.col("endpoint") == endpoint) & (pl.col("fold") == fold)
        )
        value = float(
            scorer.rae_soft_threshold_absolute_error(
                table["y_true"].to_numpy(),
                table["prediction"].to_numpy(),
                y_true_lower=table["conf_low"].to_numpy(),
                y_true_upper=table["conf_high"].to_numpy(),
            )
        )
        rows.append({"endpoint": endpoint, "fold": fold, "ST_RAE": value})
    return rows


records = []
for index, variant in enumerate(cfg["variants"]):
    trial_dir = out / f"trial-{index:02d}-{variant['name']}"
    trial_dir.mkdir()
    trial_cfg = {**cfg["base"], **{key: value for key, value in variant.items() if key != "name"}}
    trial_cfg["variant"] = variant["name"]
    trial_config_path = trial_dir / "config.json"
    trial_config_path.write_text(json.dumps(trial_cfg, indent=2) + "\n")
    os.environ["CYP_OUTPUT_DIR"] = str(trial_dir)
    os.environ["CYP_CONFIG_PATH"] = str(trial_config_path)
    medium.trainer.main()
    screen_endpoints = variant.get("screen_endpoints", endpoints)
    if not screen_endpoints or not set(screen_endpoints).issubset(set(endpoints)):
        raise ValueError(f"Invalid screen_endpoints for {variant['name']}: {screen_endpoints}")
    selection = score_embedding(
        trial_dir / "external-embeddings.parquet", cfg["selection_fold"], screen_endpoints
    )
    summary = json.loads((trial_dir / "summary.json").read_text())
    record = {
        "variant": variant["name"],
        "trial": index,
        "encoder_validation_loss": summary["best_validation_loss"],
        "best_epoch": summary["best_epoch"],
        "selection": selection,
        "selection_macro": float(np.mean([row["ST_RAE"] for row in selection])),
        "parameters": variant,
    }
    records.append(record)
    (out / "trials.json").write_text(json.dumps(records, indent=2) + "\n")
    print(json.dumps(record), flush=True)
    gc.collect()
    torch.cuda.empty_cache()

os.environ["CYP_OUTPUT_DIR"] = str(out)
os.environ["CYP_CONFIG_PATH"] = str(root_config_path)
selection_rows = [
    {"variant": record["variant"], **metric}
    for record in records
    for metric in record["selection"]
]
selection_frame = pl.DataFrame(selection_rows)
selection_frame.write_csv(out / "selection-fold3.csv")
winners = {}
confirmation_rows = []
for endpoint in endpoints:
    endpoint_rows = selection_frame.filter(pl.col("endpoint") == endpoint).sort("ST_RAE")
    winner = endpoint_rows.row(0, named=True)
    winners[endpoint] = winner["variant"]
    trial = next(record for record in records if record["variant"] == winner["variant"])
    path = out / f"trial-{trial['trial']:02d}-{trial['variant']}" / "external-embeddings.parquet"
    metric = score_embedding(path, cfg["confirmation_fold"], [endpoint])[0]
    confirmation_rows.append({"variant": winner["variant"], **metric})

references = []
for alias, name in [("previous_owned", "previous-owned.parquet"), ("jeremy", "jeremy-oof.parquet")]:
    for fold in [cfg["selection_fold"], cfg["confirmation_fold"]]:
        references.extend(
            {"reference": alias, **row} for row in reference_scores(name, fold)
        )
pl.DataFrame(confirmation_rows).write_csv(out / "confirmation-fold4.csv")
pl.DataFrame(references).write_csv(out / "references.csv")
(out / "selection.json").write_text(
    json.dumps(
        {
            "protocol": cfg["protocol"],
            "target_winners": winners,
            "confirmation": confirmation_rows,
            "references": references,
        },
        indent=2,
    )
    + "\n"
)
print(json.dumps({"target_winners": winners, "confirmation": confirmation_rows}, indent=2))
