"""Evaluate distance-based D6 applicability domain against fold-held soft error."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


def fingerprints(compounds: pl.DataFrame) -> dict[str, object]:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    result = {}
    for compound_id, smiles in compounds.select("compound_id", "standardized_smiles").iter_rows():
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"Invalid SMILES for {compound_id}")
        result[compound_id] = generator.GetFingerprint(molecule)
    return result


def tanimoto_distance(
    query: list[str], train: list[str], fps: dict[str, object], neighbors: int
) -> tuple[np.ndarray, np.ndarray]:
    training_fps = [fps[key] for key in train]
    first = np.empty(len(query), dtype=float)
    topk = np.empty(len(query), dtype=float)
    for index, key in enumerate(query):
        values = np.asarray(DataStructs.BulkTanimotoSimilarity(fps[key], training_fps))
        best = np.partition(values, -neighbors)[-neighbors:]
        first[index] = 1 - float(best.max())
        topk[index] = 1 - float(best.mean())
    return first, topk


def latent_distance(
    train: pl.DataFrame, query: pl.DataFrame, columns: list[str], components: int, neighbors: int
) -> np.ndarray:
    scaler = StandardScaler().fit(train.select(columns).to_numpy())
    train_x = scaler.transform(train.select(columns).to_numpy())
    query_x = scaler.transform(query.select(columns).to_numpy())
    pca = PCA(n_components=components, svd_solver="randomized", random_state=42)
    train_z = pca.fit_transform(train_x)
    query_z = pca.transform(query_x)
    model = NearestNeighbors(n_neighbors=neighbors, metric="euclidean", n_jobs=3).fit(train_z)
    return model.kneighbors(query_z, return_distance=True)[0].mean(axis=1)


def ranked_error(frame: pl.DataFrame, indicator: str) -> dict:
    metric = frame[indicator].to_numpy()
    error = frame["soft_error"].to_numpy()
    absolute = frame["absolute_error"].to_numpy()
    order = np.argsort(metric)
    width = max(1, len(order) // 5)
    return {
        "indicator": indicator,
        "n": len(order),
        "rho_soft_error": float(spearmanr(metric, error).statistic),
        "rho_absolute_error": float(spearmanr(metric, absolute).statistic),
        "low_fifth_mean_soft_error": float(error[order[:width]].mean()),
        "high_fifth_mean_soft_error": float(error[order[-width:]].mean()),
    }


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    endpoint = cfg["endpoint"]
    compounds = pl.read_parquet(inp / "compounds.parquet")
    test_ids = pl.read_parquet(inp / "test-mapping.parquet").select("compound_id").unique()
    direct = (
        pl.read_parquet(inp / "measurements.parquet")
        .filter(pl.col("endpoint") == endpoint)
        .select("compound_id", pl.col("value").alias("y_true"))
    )
    embedding = pl.read_parquet(inp / "medium.parquet")
    emb_cols = sorted(c for c in embedding.columns if c.startswith("emb_"))
    if direct.height != 1493 or test_ids.height != 750 or len(emb_cols) != 600:
        raise ValueError("Unexpected D6 direct, test, or embedding coverage")
    relevant_ids = pl.concat([direct.select("compound_id"), test_ids]).unique()
    relevant = relevant_ids.join(compounds, on="compound_id", validate="1:1").join(
        embedding.select("compound_id", *emb_cols), on="compound_id", validate="1:1"
    )
    if relevant.height != 2243:
        raise ValueError("Missing chemistry or embeddings")
    fps = fingerprints(relevant)
    latent = relevant.select("compound_id", *emb_cols)
    parent_members = (
        pl.read_parquet(inp / "parent-members.parquet")
        .filter(
            (pl.col("endpoint") == endpoint)
            & pl.col("scenario").is_in(["current_tabpfn", "current_ft_lgbm"])
        )
        .pivot(on="scenario", index=["compound_id", "fold"], values="prediction")
    )
    if parent_members.height != 1457:
        raise ValueError("Incomplete submitted parent model members")
    parent_members = parent_members.with_columns(
        (pl.col("current_tabpfn") - pl.col("current_ft_lgbm")).abs().alias("member_disagreement")
    ).select("compound_id", "fold", "member_disagreement")
    oof_rows = []
    test_rows = []
    fold_rows = []
    k = cfg["neighbors"]
    for split in cfg["splits"]:
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet").select(
            "compound_id", "fold", "group", "is_anchor"
        )
        labelled = direct.join(assignments, on="compound_id", validate="1:1")
        oof = (
            pl.read_parquet(inp / f"{split}-oof.parquet")
            .filter(pl.col("endpoint") == endpoint)
            .select("compound_id", "fold", "prediction", "y_true", "conf_low", "conf_high")
        )
        if oof.height != 1457 or oof.select("compound_id", "fold").is_duplicated().any():
            raise ValueError(f"Unexpected {split} OOF coverage")
        if split == "parent":
            oof = oof.join(parent_members, on=["compound_id", "fold"], validate="1:1")
        for fold in range(5):
            training = labelled.filter((pl.col("fold") != fold) | pl.col("is_anchor")).sort(
                "compound_id"
            )
            held = (
                oof.filter(pl.col("fold") == fold)
                .join(
                    assignments.select("compound_id", "group", "is_anchor"),
                    on="compound_id",
                    validate="1:1",
                )
                .sort("compound_id")
            )
            if held.filter(pl.col("is_anchor")).height or training.height + held.height != 1493:
                raise ValueError(f"Invalid {split}/fold{fold} membership")
            train_latent = (
                training.select("compound_id")
                .join(latent, on="compound_id", validate="1:1")
                .sort("compound_id")
            )
            held_latent = (
                held.select("compound_id")
                .join(latent, on="compound_id", validate="1:1")
                .sort("compound_id")
            )
            test_latent = test_ids.join(latent, on="compound_id", validate="1:1").sort(
                "compound_id"
            )
            train_keys = training["compound_id"].to_list()
            held_keys = held["compound_id"].to_list()
            test_keys = test_latent["compound_id"].to_list()
            held_top1, held_top5 = tanimoto_distance(held_keys, train_keys, fps, k)
            test_top1, test_top5 = tanimoto_distance(test_keys, train_keys, fps, k)
            held_latent_d = latent_distance(
                train_latent, held_latent, emb_cols, cfg["latent_components"], k
            )
            test_latent_d = latent_distance(
                train_latent, test_latent, emb_cols, cfg["latent_components"], k
            )
            y_fit = training["y_true"].to_numpy()
            pred = held["prediction"].to_numpy()
            low = held["conf_low"].to_numpy()
            high = held["conf_high"].to_numpy()
            truth = held["y_true"].to_numpy()
            oof_rows.append(
                held.select(
                    "compound_id",
                    "fold",
                    "group",
                    "y_true",
                    "conf_low",
                    "conf_high",
                    "prediction",
                    *(["member_disagreement"] if split == "parent" else []),
                ).with_columns(
                    pl.lit(split).alias("split_id"),
                    pl.Series("tanimoto_top1_dissimilarity", held_top1),
                    pl.Series("tanimoto_top5_dissimilarity", held_top5),
                    pl.Series("latent_top5_distance", held_latent_d),
                    pl.Series("prediction_abs_z", np.abs((pred - y_fit.mean()) / y_fit.std())),
                    pl.Series("soft_error", np.maximum(low - pred, 0) + np.maximum(pred - high, 0)),
                    pl.Series("absolute_error", np.abs(truth - pred)),
                )
            )
            test_rows.append(
                test_latent.select("compound_id").with_columns(
                    pl.lit(split).alias("split_id"),
                    pl.lit(fold).alias("fold"),
                    pl.Series("tanimoto_top1_dissimilarity", test_top1),
                    pl.Series("tanimoto_top5_dissimilarity", test_top5),
                    pl.Series("latent_top5_distance", test_latent_d),
                )
            )
            for metric, held_values, test_values in [
                ("tanimoto_top1_dissimilarity", held_top1, test_top1),
                ("tanimoto_top5_dissimilarity", held_top5, test_top5),
                ("latent_top5_distance", held_latent_d, test_latent_d),
            ]:
                fold_rows.append(
                    {
                        "split_id": split,
                        "fold": fold,
                        "indicator": metric,
                        "oof_n": len(held_values),
                        "test_n": len(test_values),
                        "oof_median": float(np.median(held_values)),
                        "test_median": float(np.median(test_values)),
                        "oof_p90": float(np.quantile(held_values, 0.9)),
                        "test_above_oof_p90": float(
                            (test_values > np.quantile(held_values, 0.9)).mean()
                        ),
                    }
                )
            print(split, fold, flush=True)
    oof_ad = pl.concat(oof_rows, how="diagonal").sort("split_id", "compound_id")
    test_ad = pl.concat(test_rows).sort("split_id", "fold", "compound_id")
    if oof_ad.height != 2914 or test_ad.height != 7500:
        raise ValueError("Incomplete AD outputs")
    association_rows = []
    risk_rows = []
    for split in cfg["splits"]:
        frame = oof_ad.filter(pl.col("split_id") == split)
        indicators = [
            "tanimoto_top1_dissimilarity",
            "tanimoto_top5_dissimilarity",
            "latent_top5_distance",
            "prediction_abs_z",
        ] + (["member_disagreement"] if split == "parent" else [])
        for indicator in indicators:
            association_rows.append(
                {"split_id": split, "fold": -1, **ranked_error(frame, indicator)}
            )
            for fold in range(5):
                part = frame.filter(pl.col("fold") == fold)
                association_rows.append(
                    {"split_id": split, "fold": fold, **ranked_error(part, indicator)}
                )
        xcols = cfg["error_features"]
        for fold in range(5):
            fit = frame.filter(pl.col("fold") != fold)
            valid = frame.filter(pl.col("fold") == fold)
            settings = cfg["error_model"]
            model = RandomForestRegressor(
                n_estimators=settings["n_estimators"],
                min_samples_leaf=settings["min_samples_leaf"],
                max_features=settings["max_features"],
                n_jobs=settings["n_jobs"],
                random_state=settings["seed"] + fold,
            )
            model.fit(fit.select(xcols).to_numpy(), fit["soft_error"].to_numpy())
            risk_rows.append(
                valid.select("compound_id", "fold", "soft_error").with_columns(
                    pl.lit(split).alias("split_id"),
                    pl.Series("predicted_error", model.predict(valid.select(xcols).to_numpy())),
                )
            )
    risk = pl.concat(risk_rows).sort("split_id", "compound_id")
    for split in cfg["splits"]:
        part = (
            risk.filter(pl.col("split_id") == split)
            .rename({"predicted_error": "learned_error_score"})
            .join(
                oof_ad.filter(pl.col("split_id") == split).select("compound_id", "absolute_error"),
                on="compound_id",
                validate="1:1",
            )
        )
        association_rows.append(
            {
                "split_id": split,
                "fold": -1,
                **ranked_error(part, "learned_error_score"),
            }
        )
    oof_ad.write_parquet(out / "oof-ad.parquet")
    test_ad.write_parquet(out / "test-ad.parquet")
    risk.write_parquet(out / "crossfit-error-score.parquet")
    pl.DataFrame(association_rows).sort("split_id", "indicator", "fold").write_csv(
        out / "error-associations.csv"
    )
    pl.DataFrame(fold_rows).sort("split_id", "indicator", "fold").write_csv(
        out / "blind-vs-oof-coverage.csv"
    )
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "oof_rows": oof_ad.height,
                "test_fold_rows": test_ad.height,
                "scope": "UNIQUE-style indicators without calibrated uncertainty",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
