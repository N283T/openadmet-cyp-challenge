"""LightGBM screen with source-to-hidden class-prior importance weights."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
target_priors = json.loads((inp / cfg["lb_summary"]).read_text())["leaderboard_priors"]


def read(name):
    return pl.read_parquet(inp / cfg[name])


def matrix(ids, table, prefix=None):
    columns = [
        column
        for column in table.columns
        if column != "compound_id" and (prefix is None or column.startswith(prefix))
    ]
    result = (
        ids.join(
            table.select("compound_id", *columns),
            on="compound_id",
            how="left",
            validate="1:1",
            maintain_order="left",
        )
        .select(columns)
        .to_numpy()
        .astype(np.float32)
    )
    result[~np.isfinite(result)] = np.nan
    return result


def expected_mcc(y, binary, prevalence):
    positive = y == 1
    tpr = float(binary[positive].mean())
    fpr = float(binary[~positive].mean())
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator > 0 else 0.0


def choose_threshold(y, probability, prevalence):
    order = np.argsort(-probability, kind="stable")
    ordered_probability = probability[order]
    positions = np.r_[
        0,
        np.flatnonzero(ordered_probability[:-1] != ordered_probability[1:]) + 1,
        len(probability),
    ]
    choices = []
    for position in positions:
        binary = np.zeros(len(y), dtype=bool)
        binary[order[:position]] = True
        if position == 0:
            threshold = np.nextafter(ordered_probability[0], np.inf)
        else:
            threshold = ordered_probability[position - 1]
        choices.append(
            (
                expected_mcc(y, binary, prevalence),
                -abs(float(threshold) - 0.5),
                float(threshold),
            )
        )
    return max(choices)[2]


compounds = read("compounds").sort("compound_id")
ids = compounds.select("compound_id")
labels = read("labels").filter(pl.col("endpoint").is_in(cfg["endpoints"]))
folds = (
    read("assignments")
    .filter(pl.col("strategy") == "umap")
    .select("compound_id", "fold", "group")
)
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
index = {compound_id: row for row, compound_id in enumerate(ids["compound_id"])}
test_index = np.asarray([index[compound_id] for compound_id in test["compound_id"]])

mordred = matrix(
    ids,
    pl.concat([pl.read_parquet(inp / filename) for filename in cfg["mordred"]]),
)
class_embedding = matrix(ids, read("class_emb"), "emb_")
mixed_embedding = matrix(ids, read("mixed_emb"), "emb_")
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
).astype(np.float32)
views = {
    "morgan": morgan,
    "mordred": mordred,
    "fp_mordred": np.column_stack([morgan, mordred]),
    "external_class": class_embedding,
    "external_mixed": mixed_embedding,
    "external_class_mordred": np.column_stack([class_embedding, mordred]),
}

oof_rows = []
inner_rows = []
test_rows = []
model_rows = []
for (endpoint,), endpoint_labels in labels.group_by("endpoint"):
    data = endpoint_labels.join(folds, on="compound_id", validate="1:1").sort("compound_id")
    y = data["value"].to_numpy().astype(int)
    fold_ids = data["fold"].to_numpy()
    data_index = np.asarray([index[compound_id] for compound_id in data["compound_id"]])
    source_prior = float(y.mean())
    target_prior = float(target_priors[endpoint])
    positive_weight = (target_prior / (1 - target_prior)) / (
        source_prior / (1 - source_prior)
    )
    for method, features in views.items():
        for fold in range(5):
            train = fold_ids != fold
            validation = ~train
            inner = fold_ids == (fold + 1) % 5
            parameters = {
                "objective": "binary",
                "n_estimators": 400,
                "num_leaves": 15,
                "min_child_samples": 20,
                "learning_rate": 0.03,
                "reg_lambda": 2,
                "n_jobs": 3,
                "random_state": 42 + fold,
                "verbosity": -1,
                "deterministic": True,
                "force_col_wise": True,
                "scale_pos_weight": positive_weight,
            }
            provisional = lgb.LGBMClassifier(**parameters).fit(
                features[data_index[train & ~inner]], y[train & ~inner]
            )
            inner_probability = provisional.predict_proba(features[data_index[inner]])[:, 1]
            provisional.booster_.save_model(str(out / f"{method}-{endpoint}-{fold}-inner.txt"))
            inner_rows.append(
                data.filter(pl.Series(inner))
                .select("compound_id", "endpoint", pl.col("value").alias("y_true"))
                .with_columns(
                    pl.lit(fold).alias("fold"),
                    pl.lit(method).alias("method"),
                    pl.Series("prediction", inner_probability),
                )
            )
            model = lgb.LGBMClassifier(**parameters).fit(features[data_index[train]], y[train])
            probability = model.predict_proba(features[data_index[validation]])[:, 1]
            test_probability = model.predict_proba(features[test_index])[:, 1]
            model.booster_.save_model(str(out / f"{method}-{endpoint}-{fold}.txt"))
            oof_rows.append(
                data.filter(pl.Series(validation))
                .select(
                    "compound_id",
                    "endpoint",
                    "fold",
                    "group",
                    pl.col("value").alias("y_true"),
                )
                .with_columns(
                    pl.lit(method).alias("method"),
                    pl.Series("prediction", probability),
                )
            )
            test_rows.append(
                test.with_columns(
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit(fold).alias("fold"),
                    pl.lit(method).alias("method"),
                    pl.Series("prediction", test_probability),
                )
            )
            model_rows.append(
                {
                    "endpoint": endpoint,
                    "method": method,
                    "fold": fold,
                    "source_prior": source_prior,
                    "target_prior": target_prior,
                    "scale_pos_weight": positive_weight,
                    "dimensions": features.shape[1],
                    "params": parameters,
                }
            )
            print(endpoint, method, fold, flush=True)

oof = pl.concat(oof_rows)
inner = pl.concat(inner_rows)
fold_test = pl.concat(test_rows)
for table_name, table in (("oof", oof), ("inner", inner), ("test", fold_test)):
    extra_keys = (
        ("group", "y_true")
        if table_name == "oof"
        else (("y_true",) if table_name == "inner" else ())
    )
    ensemble = (
        table.group_by("compound_id", "endpoint", "fold", *extra_keys)
        .agg(pl.col("prediction").mean())
        .with_columns(pl.lit("six_view_mean").alias("method"))
    )
    if table_name == "oof":
        oof = pl.concat([oof, ensemble.select(oof.columns)])
    elif table_name == "inner":
        inner = pl.concat([inner, ensemble.select(inner.columns)])
    else:
        fold_test = pl.concat([fold_test, ensemble.select(fold_test.columns)])

metric_rows = []
cut_rows = []
deployment_rows = []
evaluated_rows = []
for (endpoint, method), method_oof in oof.group_by("endpoint", "method"):
    prevalence = float(target_priors[endpoint])
    parts = []
    cuts = []
    for fold in range(5):
        fit = inner.filter(
            (pl.col("endpoint") == endpoint)
            & (pl.col("method") == method)
            & (pl.col("fold") == fold)
        )
        validation = method_oof.filter(pl.col("fold") == fold)
        cut = choose_threshold(
            fit["y_true"].to_numpy().astype(int), fit["prediction"].to_numpy(), prevalence
        )
        cuts.append(cut)
        cut_rows.append({"endpoint": endpoint, "method": method, "fold": fold, "threshold": cut})
        parts.append(validation.with_columns(pl.lit(cut).alias("threshold")))
    evaluated = pl.concat(parts)
    evaluated_rows.append(evaluated)
    y = evaluated["y_true"].to_numpy().astype(int)
    probability = evaluated["prediction"].to_numpy()
    binary = probability >= evaluated["threshold"].to_numpy()
    metric_rows.append(
        {
            "endpoint": endpoint,
            "method": method,
            "source_MCC": float(matthews_corrcoef(y, binary)),
            "target_expected_MCC": expected_mcc(y, binary, prevalence),
            "AUROC": float(roc_auc_score(y, probability)),
            "AP": float(average_precision_score(y, probability)),
        }
    )
    deployment_cut = float(np.median(cuts))
    blind = (
        fold_test.filter((pl.col("endpoint") == endpoint) & (pl.col("method") == method))
        .group_by("compound_id", "endpoint", "method")
        .agg(pl.col("prediction").mean())
        .with_columns(
            pl.lit(deployment_cut).alias("threshold"),
            (pl.col("prediction") >= deployment_cut).alias("binary"),
        )
    )
    deployment_rows.append(blind)

pl.concat(evaluated_rows).write_parquet(out / "oof.parquet")
inner.write_parquet(out / "inner.parquet")
fold_test.write_parquet(out / "fold-test-probabilities.parquet")
pl.concat(deployment_rows).write_parquet(out / "deployment.parquet")
metrics = pl.DataFrame(metric_rows).sort(
    "endpoint", "target_expected_MCC", descending=[False, True]
)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(cut_rows).write_csv(out / "thresholds.csv")
(out / "models.json").write_text(json.dumps(model_rows, indent=2) + "\n")
(out / "summary.json").write_text(
    json.dumps(
        {
            "target_priors": target_priors,
            "methods": list(views) + ["six_view_mean"],
            "weighting": "target odds divided by source odds",
            "threshold": "Expected target-prior MCC on isolated inner fold",
            "split_run": cfg["split_run"],
        },
        indent=2,
    )
    + "\n"
)
print(metrics.write_csv(), flush=True)
