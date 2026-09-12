"""Twenty fixed external assay tree models; no challenge labels or GPU."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, log_loss, roc_auc_score


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def path(name):
        files = [p for p in inp.iterdir() if p.name.partition("-")[2] == name]
        assert len(files) == 1
        return files[0]

    features = np.load(path("molecular-features.npy"), mmap_mode="r")
    ids = pl.read_parquet(path("feature-ids.parquet"))
    compounds = pl.read_parquet(path("external-compounds.parquet")).sort("compound_id")
    observations = pl.read_parquet(path("observations.parquet")).filter(
        pl.col("source").is_in(["pubchem1851", "tox21_p450_glo"])
    )
    assert not observations["challenge_overlap"].any()
    challenge = ids.filter(pl.col("role") == "challenge").sort("compound_id")
    assert not set(compounds["compound_id"]) & set(challenge["compound_id"])
    assert features.shape[0] == ids.height and ids["compound_id"].n_unique() == ids.height
    lookup = {key: i for i, key in enumerate(ids["compound_id"])}
    challenge_indices = np.array([lookup[key] for key in challenge["compound_id"]])
    rng = np.random.default_rng(cfg.get("seed", 42))
    groups = rng.permutation(sorted(compounds["connectivity"].unique()))
    val_groups = set(groups[: max(1, int(len(groups) * 0.1))])
    membership = compounds.select("compound_id", "connectivity").with_columns(
        pl.Series(
            "role",
            [
                "validation" if key in val_groups else "training"
                for key in compounds["connectivity"]
            ],
        )
    )
    membership.write_parquet(out / "membership.parquet")
    assays = observations.select("assay_id", "source", "isoform").unique().sort("assay_id")
    assert assays.height == 10 and assays["assay_id"].n_unique() == 10
    scores, metrics, prediction_rows, metadata = {}, [], [], []
    params = dict(
        n_estimators=400,
        num_leaves=15,
        min_child_samples=20,
        learning_rate=0.05,
        n_jobs=4,
        random_state=cfg.get("seed", 42),
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
    )
    for assay in assays.iter_rows(named=True):
        subset = observations.filter(pl.col("assay_id") == assay["assay_id"])
        # Exact potency replicates are median-aggregated; censored values are
        # not invented. Any contradictory or fractional binary labels are masked.
        subset = (
            subset.group_by("compound_id")
            .agg(
                pl.col("value")
                .filter((pl.col("relation") == "=") & pl.col("value").is_finite())
                .median(),
                pl.when(
                    (pl.col("binary_label").drop_nulls().n_unique() == 1)
                    & pl.col("binary_label").drop_nulls().is_in([0, 1]).all()
                )
                .then(pl.col("binary_label").drop_nulls().first())
                .otherwise(None)
                .alias("binary_label"),
            )
            .join(membership, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        for objective, label in [("binary", "binary_label"), ("continuous", "value")]:
            data = subset.filter(pl.col(label).is_not_null() & pl.col(label).is_finite())
            indices = np.array([lookup[key] for key in data["compound_id"]])
            validation = data["role"].to_numpy() == "validation"
            y = data[label].to_numpy()
            assert validation.any() and (~validation).any()
            name = "score_" + assay["source"] + "_" + assay["isoform"] + "_" + objective
            if objective == "binary":
                assert len(np.unique(y[~validation])) == 2
                model = lgb.LGBMClassifier(objective="binary", **params)
            else:
                model = lgb.LGBMRegressor(objective="regression_l1", **params)
            model.fit(features[indices[~validation]], y[~validation])
            if objective == "binary":
                pv = model.predict_proba(features[indices[validation]])[:, 1]
                pc = model.predict_proba(features[challenge_indices])[:, 1]
                stats = dict(
                    log_loss=float(log_loss(y[validation], pv, labels=[0, 1])),
                    roc_auc=float(roc_auc_score(y[validation], pv)),
                    average_precision=float(average_precision_score(y[validation], pv)),
                )
            else:
                pv = model.predict(features[indices[validation]])
                pc = model.predict(features[challenge_indices])
                stats = dict(
                    mae=float(np.abs(pv - y[validation]).mean()),
                    rmse=float(np.sqrt(((pv - y[validation]) ** 2).mean())),
                    spearman=float(spearmanr(y[validation], pv).statistic),
                )
            assert np.isfinite(pv).all() and np.isfinite(pc).all()
            scores[name] = pc
            model.booster_.save_model(str(out / (name + ".lgbm.txt")))
            info = dict(
                name=name,
                **assay,
                objective=objective,
                training=int((~validation).sum()),
                validation=int(validation.sum()),
            )
            metrics.extend(
                dict(**info, metric=metric, value=value) for metric, value in stats.items()
            )
            metadata.append(info)
            prediction_rows.append(
                data.filter(pl.Series(validation))
                .select("compound_id", "connectivity")
                .with_columns(
                    pl.lit(assay["assay_id"]).alias("assay_id"),
                    pl.lit(objective).alias("objective"),
                    pl.lit("external_connectivity_10pct").alias("split_id"),
                    pl.lit("validation").alias("role"),
                    pl.Series("y_true", y[validation]),
                    pl.Series("prediction", pv),
                )
            )
            print(json.dumps(dict(**info, **stats)), flush=True)
    assert len(scores) == 20
    pl.DataFrame(scores).insert_column(0, challenge["compound_id"]).write_parquet(
        out / "external-scores.parquet"
    )
    pl.DataFrame(metrics).write_csv(out / "external-metrics.csv")
    pl.concat(prediction_rows).write_parquet(out / "external-validation-predictions.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                models=metadata,
                parameters=params,
                challenge_rows=challenge.height,
                score_columns=list(scores),
                challenge_labels_used=False,
                external_refit=False,
                training_policy="Fixed 400 trees on external training 90%; holdout not refitted",
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
