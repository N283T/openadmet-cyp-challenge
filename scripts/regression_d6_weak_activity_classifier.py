"""Estimate weak-activity classification signal with outer-fold exclusion."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from catboost import CatBoostClassifier
from sklearn.metrics import average_precision_score, roc_auc_score


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    labels = pl.read_parquet(inp / "measurements.parquet").filter(
        pl.col("endpoint") == cfg["endpoint"]
    ).select("compound_id", "endpoint", pl.col("value").alias("y_true"), "conf_high")
    scores = pl.read_parquet(inp / "scores.parquet")
    medium = pl.read_parquet(inp / "medium.parquet")
    score_cols = sorted(c for c in scores.columns if c.startswith("score_"))
    medium_cols = sorted(c for c in medium.columns if c.startswith("emb_"))
    if len(score_cols) != 47 or len(medium_cols) != 600 or labels.height != 1493:
        raise ValueError("Unexpected feature or label coverage")
    features = scores.select("compound_id", *score_cols).join(
        medium.select("compound_id", *medium_cols), on="compound_id", validate="1:1"
    )
    views = {
        "scores47": score_cols,
        "scores47_medium600": score_cols + medium_cols,
    }
    outputs = []
    for split in cfg["splits"]:
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet")
        data = labels.join(
            assignments.select("compound_id", "fold", "group", "is_anchor"),
            on="compound_id", validate="1:1",
        ).join(features, on="compound_id", validate="1:1").sort("compound_id")
        if data.height != 1493 or data.filter(~pl.col("is_anchor")).height != 1457:
            raise ValueError(f"Incomplete {split} alignment")
        fold_ids = data["fold"].to_numpy()
        anchors = data["is_anchor"].to_numpy()
        y = data["y_true"].to_numpy()
        hi = data["conf_high"].to_numpy()
        for view in cfg["feature_views"]:
            matrix = data.select(views[view]).to_numpy().astype(np.float32)
            if not np.isfinite(matrix).all():
                raise ValueError(f"Nonfinite {view} features")
            for rule in cfg["label_rules"]:
                weak = (y if rule == "point" else hi) < cfg["weak_cutoff"]
                for fold in range(5):
                    train = (fold_ids != fold) | anchors
                    valid = (fold_ids == fold) & ~anchors
                    if weak[train].sum() < 20 or weak[valid].sum() < 1:
                        raise ValueError(f"Insufficient {split}/{rule}/fold{fold} positives")
                    settings = cfg["model"]
                    model = CatBoostClassifier(
                        iterations=settings["iterations"], depth=settings["depth"],
                        learning_rate=settings["learning_rate"],
                        l2_leaf_reg=settings["l2_leaf_reg"],
                        loss_function="Logloss", random_seed=settings["seed"] + fold,
                        thread_count=settings["thread_count"], verbose=False,
                        allow_writing_files=False,
                    )
                    model.fit(matrix[train], weak[train].astype(int))
                    probability = model.predict_proba(matrix[valid])[:, 1]
                    outputs.append(data.filter(pl.Series(valid)).select(
                        "compound_id", "endpoint", "fold", "group", "y_true", "conf_high"
                    ).with_columns(
                        pl.lit(split).alias("split_id"),
                        pl.lit(view).alias("view"),
                        pl.lit(rule).alias("label_rule"),
                        pl.Series("probability", probability),
                    ))
                    print(split, view, rule, fold, flush=True)
    oof = pl.concat(outputs).sort("split_id", "view", "label_rule", "compound_id")
    expected = len(cfg["splits"]) * len(cfg["feature_views"]) * len(cfg["label_rules"]) * 1457
    if oof.height != expected or oof.select(
        "split_id", "view", "label_rule", "compound_id"
    ).is_duplicated().any():
        raise ValueError("Incomplete classifier OOF")
    rows = []
    for (split, view, rule), frame in oof.group_by("split_id", "view", "label_rule"):
        truth = (
            frame["y_true"].to_numpy() if rule == "point" else frame["conf_high"].to_numpy()
        ) < cfg["weak_cutoff"]
        probability = frame["probability"].to_numpy()
        selected = np.argsort(probability)[-int(np.ceil(len(probability) * 0.1)):]
        rows.append({
            "split_id": split, "view": view, "label_rule": rule,
            "n": len(truth), "positives": int(truth.sum()),
            "prevalence": float(truth.mean()),
            "AUROC": float(roc_auc_score(truth, probability)),
            "AP": float(average_precision_score(truth, probability)),
            "top_decile_n": len(selected),
            "top_decile_positive_rate": float(truth[selected].mean()),
            "top_decile_lift": float(truth[selected].mean() / truth.mean()),
            "top_decile_mean_pic50": float(frame["y_true"].to_numpy()[selected].mean()),
        })
    pl.DataFrame(rows).sort("split_id", "label_rule", "view").write_csv(out / "metrics.csv")
    oof.write_parquet(out / "oof.parquet")
    (out / "summary.json").write_text(json.dumps({
        "config": cfg, "rows": oof.height,
        "qualification": "Exploratory classification signal, not a submission recipe",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
