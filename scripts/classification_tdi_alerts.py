"""Matched LightGBM ablation for explicit TDI structural-alert features."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

from cyp_challenge.tdi_alerts import TDI_ALERTS, alert_matrix

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name):
    return pl.read_parquet(inp / cfg[name])


def table_matrix(ids, table):
    columns = [column for column in table.columns if column != "compound_id"]
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


def select_threshold(y, probability):
    order = np.argsort(-probability, kind="stable")
    ordered_y = y[order]
    ordered_probability = probability[order]
    cut_positions = np.r_[
        np.flatnonzero(ordered_probability[:-1] != ordered_probability[1:]) + 1,
        len(probability),
    ]
    true_positive = np.cumsum(ordered_y)[cut_positions - 1]
    false_positive = cut_positions - true_positive
    false_negative = y.sum() - true_positive
    true_negative = len(y) - y.sum() - false_positive
    denominator = np.sqrt(
        (true_positive + false_positive)
        * (true_positive + false_negative)
        * (true_negative + false_positive)
        * (true_negative + false_negative)
    )
    scores = np.divide(
        true_positive * true_negative - false_positive * false_negative,
        denominator,
        out=np.zeros_like(denominator, dtype=float),
        where=denominator > 0,
    )
    candidates = [
        (float(score), -abs(float(cut) - 0.5), float(cut))
        for score, cut in zip(scores, ordered_probability[cut_positions - 1], strict=True)
    ]
    candidates.append((float(matthews_corrcoef(y, probability >= 0.5)), 0.0, 0.5))
    return max(candidates)[2]


def metrics(y, probability, binary):
    return {
        "MCC": float(matthews_corrcoef(y, binary)),
        "AUROC": float(roc_auc_score(y, probability)),
        "AP": float(average_precision_score(y, probability)),
        "positives": int(binary.sum()),
        "n": len(y),
    }


reference = (inp / cfg["alert_reference"]).read_text()
assert "Bioactivation alerts for TDI" in reference

compounds = read("compounds").sort("compound_id")
ids = compounds.select("compound_id")
labels = read("labels").filter(pl.col("endpoint").is_in(cfg["endpoints"]))
folds = (
    read("assignments")
    .filter(pl.col("strategy") == "umap")
    .select("compound_id", "fold", "group")
)
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
assert len(labels) == 5081
assert len(test) == 750

compound_index = {compound_id: index for index, compound_id in enumerate(ids["compound_id"])}
test_index = np.asarray([compound_index[compound_id] for compound_id in test["compound_id"]])

mordred_table = pl.concat([pl.read_parquet(inp / filename) for filename in cfg["mordred"]])
mordred = table_matrix(ids, mordred_table)

generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
molecules = [Chem.MolFromSmiles(smiles) for smiles in compounds["standardized_smiles"]]
assert all(molecule is not None for molecule in molecules)
morgan = np.stack(
    [generator.GetFingerprintAsNumPy(molecule) for molecule in molecules]
).astype(np.float32)
alerts, alert_names = alert_matrix(compounds["standardized_smiles"].to_list())

views = {
    "morgan_balanced": morgan,
    "morgan_alerts_balanced": np.column_stack([morgan, alerts]),
    "fp_mordred_balanced": np.column_stack([morgan, mordred]),
    "fp_mordred_alerts_balanced": np.column_stack([morgan, mordred, alerts]),
}

alert_table = ids.with_columns(
    *[pl.Series(f"alert_{name}", alerts[:, index]) for index, name in enumerate(alert_names)]
)
alert_table.write_parquet(out / "alert-features.parquet")

prevalence_rows = []
for (endpoint,), endpoint_labels in labels.group_by("endpoint"):
    joined = endpoint_labels.join(alert_table, on="compound_id", validate="m:1")
    y = joined["value"].to_numpy().astype(int)
    base_rate = float(y.mean())
    for name in alert_names[:-2]:
        fired = joined[f"alert_{name}"].to_numpy() > 0
        positive = int(y[fired].sum())
        precision = float(y[fired].mean()) if fired.any() else None
        prevalence_rows.append(
            {
                "endpoint": endpoint,
                "alert": name,
                "n_fired": int(fired.sum()),
                "positive": positive,
                "base_rate": base_rate,
                "precision": precision,
                "enrichment": precision / base_rate if precision is not None else None,
            }
        )
pl.DataFrame(prevalence_rows).write_csv(out / "alert-prevalence.csv")

oof_rows = []
test_rows = []
metric_rows = []
model_rows = []
membership_rows = []
for (endpoint,), endpoint_labels in labels.group_by("endpoint"):
    endpoint_data = endpoint_labels.join(folds, on="compound_id", validate="1:1").sort(
        "compound_id"
    )
    y = endpoint_data["value"].to_numpy().astype(int)
    fold_ids = endpoint_data["fold"].to_numpy()
    endpoint_index = np.asarray(
        [compound_index[compound_id] for compound_id in endpoint_data["compound_id"]]
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
                "class_weight": "balanced",
            }
            provisional = lgb.LGBMClassifier(**parameters).fit(
                features[endpoint_index[train & ~inner]], y[train & ~inner]
            )
            inner_probability = provisional.predict_proba(features[endpoint_index[inner]])[:, 1]
            threshold = select_threshold(y[inner], inner_probability)
            provisional.booster_.save_model(str(out / f"{method}-{endpoint}-{fold}-inner.txt"))

            model = lgb.LGBMClassifier(**parameters).fit(
                features[endpoint_index[train]], y[train]
            )
            probability = model.predict_proba(features[endpoint_index[validation]])[:, 1]
            test_probability = model.predict_proba(features[test_index])[:, 1]
            assert np.isfinite(probability).all() and np.isfinite(test_probability).all()
            model.booster_.save_model(str(out / f"{method}-{endpoint}-{fold}.txt"))

            oof_rows.append(
                endpoint_data.filter(pl.Series(validation))
                .select(
                    "compound_id",
                    "endpoint",
                    "fold",
                    "group",
                    pl.col("value").alias("y_true"),
                )
                .with_columns(
                    pl.Series("prediction", probability),
                    pl.lit(method).alias("method"),
                    pl.lit(threshold).alias("threshold"),
                    pl.lit(cfg["split_run"]).alias("split_id"),
                )
            )
            test_rows.append(
                test.with_columns(
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit(fold).alias("fold"),
                    pl.lit(method).alias("method"),
                    pl.Series("prediction", test_probability),
                    pl.lit(threshold).alias("threshold"),
                    pl.lit(cfg["split_run"]).alias("split_id"),
                )
            )
            for policy, cut in (("fixed_05", 0.5), ("inner_mcc", threshold)):
                metric_rows.append(
                    {
                        "endpoint": endpoint,
                        "method": method,
                        "fold": fold,
                        "policy": policy,
                        **metrics(y[validation], probability, probability >= cut),
                    }
                )
            model_rows.append(
                {
                    "endpoint": endpoint,
                    "method": method,
                    "fold": fold,
                    "threshold": threshold,
                    "n_train": int(train.sum()),
                    "n_inner_fit": int((train & ~inner).sum()),
                    "n_inner_validation": int(inner.sum()),
                    "dimensions": features.shape[1],
                    "params": parameters,
                }
            )
            if method == "morgan_balanced":
                membership_rows.append(
                    endpoint_data.select("compound_id", "endpoint", "fold").with_columns(
                        pl.lit(fold).alias("outer_fold"),
                        pl.Series("outer_validation", validation),
                        pl.Series("inner_validation", inner),
                    )
                )
            print(endpoint, method, fold, flush=True)

oof = pl.concat(oof_rows)
oof.write_parquet(out / "oof.parquet")
pl.concat(test_rows).write_parquet(out / "fold-test-probabilities.parquet")
for (endpoint, method), method_data in oof.group_by("endpoint", "method"):
    y = method_data["y_true"].to_numpy().astype(int)
    probability = method_data["prediction"].to_numpy()
    for policy, cut in (
        ("fixed_05", 0.5),
        ("inner_mcc", method_data["threshold"].to_numpy()),
    ):
        metric_rows.append(
            {
                "endpoint": endpoint,
                "method": method,
                "fold": -1,
                "policy": policy,
                **metrics(y, probability, probability >= cut),
            }
        )

metric_table = pl.DataFrame(metric_rows)
metric_table.write_csv(out / "metrics.csv")
pl.concat(membership_rows).write_parquet(out / "membership.parquet")
(out / "models.json").write_text(json.dumps(model_rows, indent=2))
(out / "summary.json").write_text(
    json.dumps(
        {
            "methods": list(views),
            "folds": 5,
            "models": len(model_rows),
            "alert_features": alert_names,
            "alerts": [alert.__dict__ for alert in TDI_ALERTS],
            "alert_implementation": (
                "Independent RDKit approximation of the categories in the archived "
                "OpenADMET article; not the article's unavailable exact SMARTS"
            ),
            "threshold": "Fixed inner UMAP fold; outer fold excluded from fit and threshold choice",
            "split_run": cfg["split_run"],
        },
        indent=2,
    )
    + "\n"
)
print(
    metric_table.filter(pl.col("fold") == -1)
    .sort("endpoint", "policy", "MCC", descending=[False, False, True])
    .write_csv(),
    flush=True,
)
