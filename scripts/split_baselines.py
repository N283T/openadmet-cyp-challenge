"""Fixed LightGBM split probe, with durable keyed OOF and structural diagnostics."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

from cyp_challenge.storage import write_json


def metric_table(rows):
    """Regression rows precede TDI rows; infer the complete union of metric columns."""
    return pl.DataFrame(rows, infer_schema_length=None)


def st_rae(y, prediction, lower, upper):
    def error(p):
        return np.maximum(p - upper, 0) + np.maximum(lower - p, 0)

    denominator = error(np.mean(y)).sum()
    return float(error(prediction).sum() / denominator) if denominator > 0 else None


def metrics(y, p, lower, upper, classification):
    if classification:
        return {
            "MCC": float(matthews_corrcoef(y, p >= 0.5)),
            "AUROC": float(roc_auc_score(y, p)),
            "AP": float(average_precision_score(y, p)),
        }
    return {
        "ST-RAE": st_rae(y, p, lower, upper),
        "MAE": float(np.abs(y - p).mean()),
        "Spearman": float(spearmanr(y, p).statistic),
    }


def nearest(query, reference):
    reference = reference.astype(np.float32)
    sizes = reference.sum(axis=1)
    result = []
    for start in range(0, len(query), 128):
        q = query[start : start + 128].astype(np.float32)
        intersection = q @ reference.T
        union = q.sum(axis=1)[:, None] + sizes - intersection
        similarity = np.divide(intersection, union, out=np.ones_like(union), where=union != 0)
        result.extend(similarity.max(axis=1).tolist())
    return np.array(result)


def main():
    inputs = Path(os.environ["CYP_INPUT_DIR"])
    output = Path(os.environ["CYP_OUTPUT_DIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def read(name):
        return pl.read_parquet(next(inputs.glob("*-" + name + ".parquet")))

    labels = read("measurements")
    aliases = read("aliases")
    assignments = read("assignments")
    desc = pl.concat([pl.read_parquet(p) for p in sorted(inputs.glob("*-mordred-*.parquet"))])
    assert desc["compound_id"].n_unique() == desc.height
    features = [c for c in desc.columns if c != "compound_id"]
    fp = np.load(next(inputs.glob("*-fingerprints.npz")), allow_pickle=False)
    fpindex = {cid: i for i, cid in enumerate(fp["ids"])}
    fingerprint = fp["fingerprints"]
    test_ids = sorted(set(aliases.filter(pl.col("role") == "test")["compound_id"]))
    testfp = fingerprint[[fpindex[c] for c in test_ids]]
    names = (
        aliases.select("compound_id", "molecule_name")
        .drop_nulls()
        .unique()
        .group_by("compound_id")
        .agg(pl.col("molecule_name").sort().str.join("|"))
    )
    predictions, scores, diagnostics = [], [], []
    (output / "models").mkdir()
    write_json(output / "effective_config.json", config)
    for endpoint in config["endpoints"]:
        frame = labels.filter(pl.col("endpoint") == endpoint).sort("compound_id")
        assert frame["compound_id"].n_unique() == frame.height
        frame = frame.join(desc, on="compound_id", how="left", validate="1:1").join(
            names, on="compound_id", how="left", validate="1:1"
        )
        ids = frame["compound_id"].to_list()
        y = frame["value"].to_numpy()
        lower, upper = frame["conf_low"].to_numpy(), frame["conf_high"].to_numpy()
        classification = endpoint.endswith("is_TDI")
        assert np.isfinite(y).all()
        if classification:
            assert set(y) == {0, 1}
        else:
            assert np.isfinite(lower).all() and np.isfinite(upper).all()
            assert (lower <= upper).all()
        x = frame.select(features).to_numpy().astype(np.float32)
        x[~np.isfinite(x)] = np.nan
        f = fingerprint[[fpindex[c] for c in ids]]
        test_nn = nearest(testfp, f)
        diagnostics.append(
            {
                "endpoint": endpoint,
                "strategy": "test_vs_all_labeled",
                "fold": -1,
                "n": len(test_nn),
                "nn_mean": float(test_nn.mean()),
                "nn_median": float(np.median(test_nn)),
                "nn_gt_05": float((test_nn > 0.5).mean()),
            }
        )
        for strategy in config["strategies"]:
            split = assignments.filter(pl.col("strategy") == strategy)
            assert split["compound_id"].n_unique() == split.height
            if strategy != "random":
                assert split.group_by("group").agg(pl.col("fold").n_unique())["fold"].max() == 1
            folds = frame.select("compound_id").join(
                split, on="compound_id", how="left", validate="1:1"
            )["fold"]
            assert folds.null_count() == 0
            folds = folds.to_numpy()
            for fold in range(config["n_folds"]):
                train, val = folds != fold, folds == fold
                nn = nearest(f[val], f[train])
                diagnostics.append(
                    {
                        "endpoint": endpoint,
                        "strategy": strategy,
                        "fold": fold,
                        "n": int(val.sum()),
                        "train_n": int(train.sum()),
                        "y_mean": float(y[val].mean()),
                        "y_std": float(y[val].std()),
                        "nn_mean": float(nn.mean()),
                        "nn_median": float(np.median(nn)),
                        "nn_gt_05": float((nn > 0.5).mean()),
                    }
                )
            for variant in config["features"]:
                matrix = x if variant == "mordred" else np.concatenate([x, f], axis=1)
                feature_names = (
                    features
                    if variant == "mordred"
                    else features + [f"Morgan_{i}" for i in range(f.shape[1])]
                )
                oof = np.full(len(y), np.nan)
                for fold in range(config["n_folds"]):
                    train, val = folds != fold, folds == fold
                    if classification:
                        assert len(set(y[train])) == len(set(y[val])) == 2
                    model_type = lgb.LGBMClassifier if classification else lgb.LGBMRegressor
                    model = model_type(
                        **config["lightgbm"],
                        random_state=config["seed"],
                        objective="binary" if classification else config["regression_objective"],
                    )
                    model.fit(matrix[train], y[train], feature_name=feature_names)
                    p = (
                        model.predict_proba(matrix[val])[:, 1]
                        if classification
                        else model.predict(matrix[val])
                    )
                    oof[val] = p
                    prefix = f"{endpoint}-{strategy}-{variant}-{fold}"
                    model.booster_.save_model(str(output / "models" / (prefix + ".txt")))
                    rows = (
                        frame.filter(pl.Series(val))
                        .select("compound_id", "molecule_name")
                        .with_columns(
                            pl.lit(endpoint).alias("endpoint"),
                            pl.lit(config["split_run"]).alias("split_id"),
                            pl.lit(strategy).alias("strategy"),
                            pl.lit(variant).alias("features"),
                            pl.lit(fold).alias("fold"),
                            pl.Series("y_true", y[val]),
                            pl.Series("prediction", p),
                        )
                    )
                    rows.write_parquet(output / (prefix + ".parquet"))
                    predictions.append(rows)
                    scores.append(
                        {
                            "endpoint": endpoint,
                            "strategy": strategy,
                            "features": variant,
                            "fold": fold,
                            "n": int(val.sum()),
                            **metrics(y[val], p, lower[val], upper[val], classification),
                        }
                    )
                assert np.isfinite(oof).all()
                score = {
                    "endpoint": endpoint,
                    "strategy": strategy,
                    "features": variant,
                    "fold": -1,
                    "n": len(y),
                    **metrics(y, oof, lower, upper, classification),
                }
                scores.append(score)
                print(json.dumps(score), flush=True)
                metric_table(scores).write_parquet(output / "metrics-progress.parquet")
    pl.concat(predictions).write_parquet(output / "oof.parquet")
    metric_table(scores).write_parquet(output / "metrics.parquet")
    pl.DataFrame(diagnostics).write_parquet(output / "diagnostics.parquet")
    write_json(
        output / "summary.json",
        {
            "config": config,
            "pooled_metrics": [s for s in scores if s["fold"] == -1],
            "model_count": len(predictions),
            "feature_names": features,
        },
    )


if __name__ == "__main__":
    main()
