"""Direct local label/interval predictions blended with the fixed TabPFN mean."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr

from cyp_challenge.external_training import snapshot_input
from cyp_challenge.storage import write_json


def weighted_median(values, weights):
    order = np.argsort(values)
    index = np.searchsorted(np.cumsum(weights[order]), weights.sum() / 2, side="left")
    return values[order[min(index, len(order) - 1)]]


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def read(name):
        return pl.read_parquet(snapshot_input(inp, name + ".parquet"))

    oof = read("regression-oof-diagnostics")
    test = read("ensemble_probabilities").filter(pl.col("endpoint").is_in(oof["endpoint"].unique()))
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps, keys = {}, {}
    for cid, s in read("compounds").select("compound_id", "standardized_smiles").iter_rows():
        mol = Chem.MolFromSmiles(s)
        fps[cid] = gen.GetFingerprint(mol)
        keys[cid] = Chem.MolToInchiKey(mol).split("-")[0]
    spec = importlib.util.spec_from_file_location(
        "official", snapshot_input(inp, "custom_scoring_functions.py")
    )
    official = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(official)
    predictions, testpred, stats, rows = [], [], [], []
    for ep in sorted(oof["endpoint"].unique()):
        f = oof.filter(pl.col("endpoint") == ep).sort("compound_id")
        t = test.filter(pl.col("endpoint") == ep).sort("compound_id")
        ids, tids = f["compound_id"].to_list(), t["compound_id"].to_list()
        y, lo, hi, p = [
            f[c].to_numpy().astype(float) for c in ["y_true", "conf_low", "conf_high", "prediction"]
        ]
        tp = t["prediction"].to_numpy().astype(float)
        folds = f["fold"].to_numpy()
        methods = [
            (k, power, target)
            for k in cfg["neighbors"]
            for power in cfg["powers"]
            for target in ["label", "interval"]
        ]
        qmaps = {str(m): np.zeros(len(y)) for m in methods}
        tmaps = {str(m): np.zeros(len(tids)) for m in methods}
        for fold in range(5):
            va = folds == fold
            forbidden = {keys[c] for c in np.array(ids)[va]}
            tr = np.flatnonzero((folds != fold) & np.array([keys[c] not in forbidden for c in ids]))
            refs = [fps[ids[i]] for i in tr]
            for split, query, dest in [
                ("umap", np.array(ids)[va].tolist(), np.flatnonzero(va)),
                ("test", tids, np.arange(len(tids))),
            ]:
                sim = np.array([DataStructs.BulkTanimotoSimilarity(fps[c], refs) for c in query])
                index = np.argsort(sim, axis=1)[:, ::-1][:, : max(cfg["neighbors"])]
                stats.append(
                    dict(
                        endpoint=ep,
                        fold=fold,
                        split=split,
                        mean_nearest=float(sim.max(axis=1).mean()),
                        fraction_above05=float((sim.max(axis=1) > 0.5).mean()),
                    )
                )
                for k, power, target in methods:
                    pred = []
                    for i, nn in enumerate(index[:, :k]):
                        source = tr[nn]
                        w = np.maximum(sim[i, nn], 1e-8) ** power
                        if target == "label":
                            value = weighted_median(y[source], w)
                        else:
                            value = weighted_median(np.r_[lo[source], hi[source]], np.r_[w, w])
                        pred.append(value)
                    key = str((k, power, target))
                    if split == "umap":
                        qmaps[key][dest] = pred
                    else:
                        tmaps[key] += np.asarray(pred) / 5
        for key, q in qmaps.items():
            for strength in cfg["strengths"]:
                name = f"knn{key}:weight{strength}"
                pred = (1 - strength) * p + strength * q
                target = (1 - strength) * tp + strength * tmaps[key]
                assert np.isfinite(pred).all() and np.isfinite(target).all()
                for fold in [-1, 0, 1, 2, 3, 4]:
                    ix = np.ones(len(y), dtype=bool) if fold == -1 else folds == fold
                    st = official.rae_soft_threshold_absolute_error(
                        y[ix], pred[ix], y_true_lower=lo[ix], y_true_upper=hi[ix]
                    )
                    rows.append(
                        dict(
                            endpoint=ep,
                            method=name,
                            fold=fold,
                            ST_RAE=float(st),
                            MAE=float(np.abs(pred[ix] - y[ix]).mean()),
                            Spearman=float(spearmanr(y[ix], pred[ix]).statistic),
                        )
                    )
                predictions.append(
                    f.select("compound_id", "endpoint", "fold", "split_id", "y_true").with_columns(
                        pl.lit(name).alias("method"), pl.Series("prediction", pred)
                    )
                )
                testpred.append(
                    t.select("compound_id", "endpoint").with_columns(
                        pl.lit(name).alias("method"),
                        pl.Series("prediction", target),
                        pl.lit(-1).alias("fold"),
                        pl.lit("test:five-fold-point-mean").alias("split_id"),
                    )
                )
    score = pl.DataFrame(rows)
    summary = (
        score.filter(pl.col("fold") == -1)
        .group_by("method")
        .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
        .sort("ST_RAE")
    )
    score.write_csv(out / "neighborhood-metrics.csv")
    summary.write_csv(out / "neighborhood-summary.csv")
    pl.DataFrame(stats).write_csv(out / "neighborhood-similarity.csv")
    pl.concat(predictions).write_parquet(out / "neighborhood-oof.parquet")
    pl.concat(testpred).write_parquet(out / "neighborhood-test.parquet")
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            audit_passed=True,
            results=summary.to_dicts(),
            scope="Fixed neighborhood rules use outer training labels only; "
            "test averages five training-pool predictions. Exploratory grid, no LB feedback.",
        ),
    )
    print(summary.head(10))


if __name__ == "__main__":
    main()
