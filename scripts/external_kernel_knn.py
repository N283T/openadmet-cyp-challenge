"""Nested inner-selected local SAR baselines on fixed outer folds, CPU only."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(name):
    return next(p for p in inp.iterdir() if p.name.partition("-")[2] == name)


def read(name):
    return pl.read_parquet(path(name + ".parquet"))


spec = importlib.util.spec_from_file_location("official", path("custom_scoring_functions.py"))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
threadpool_limits(limits=cfg["threads"])
a, labels, c, emb = (
    read("parent-assignments"),
    read("measurements"),
    read("compounds"),
    read("external-embeddings"),
)
embcols = sorted(x for x in emb.columns if x.startswith("emb_"))
assert len(embcols) == 600
fpgen = rdFingerprintGenerator.GetMorganGenerator(
    radius=cfg["morgan_radius"], fpSize=cfg["morgan_bits"]
)
outputs = {name: [] for name in ["morgan_tanimoto_knn", "external3_classification_cosine_knn"]}
selections, inner_scores, references, neighbors, metrics = [], [], [], [], []


def score(h, p):
    return float(
        scorer.rae_soft_threshold_absolute_error(
            h["value"].to_numpy(),
            p,
            y_true_lower=h["conf_low"].to_numpy(),
            y_true_upper=h["conf_high"].to_numpy(),
        )
    )


def predict(distance, query, reference, y, k, weight):
    dist = distance[np.ix_(query, reference)]
    k = min(k, len(reference))
    nearest = np.argsort(dist, axis=1, kind="stable")[:, :k]
    d = np.take_along_axis(dist, nearest, axis=1)
    yy = y[reference[nearest]]
    w = np.ones_like(d) if weight == "uniform" else 1 / np.maximum(d, cfg["distance_epsilon"])
    pred = (w * yy).sum(axis=1) / w.sum(axis=1)
    return pred, reference[nearest], d


for (endpoint,), target in labels.filter(
    pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
).group_by("endpoint", maintain_order=True):
    h = (
        target.join(a, on="compound_id", validate="1:1")
        .join(c, on="compound_id", validate="1:1")
        .join(emb, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    y = h["value"].to_numpy()
    ids = h["compound_id"].to_numpy()
    groups = h["group"].to_numpy()
    anchors = h["is_anchor"].to_numpy()
    fingerprints = [fpgen.GetFingerprint(Chem.MolFromSmiles(s)) for s in h["standardized_smiles"]]
    similarity = np.asarray(
        [DataStructs.BulkTanimotoSimilarity(f, fingerprints) for f in fingerprints],
        dtype=np.float32,
    )
    x = h.select(embcols).to_numpy().astype(np.float32)
    x = x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-12)
    distances = {
        "morgan_tanimoto_knn": np.clip(1 - similarity, 0, 1),
        "external3_classification_cosine_knn": np.clip(1 - x @ x.T, 0, 2),
    }
    for fold in range(5):
        outer_validation = (h["fold"].to_numpy() == fold) & ~anchors
        tr = np.flatnonzero(~outer_validation)
        va = np.flatnonzero(outer_validation)
        anchor_groups = set(groups[tr[anchors[tr]]])
        eligible = sorted(set(groups[tr]) - anchor_groups)
        rng = np.random.default_rng(cfg["seed"] + fold)
        inner_groups = set(
            rng.permutation(eligible)[: max(1, int(cfg["inner_fraction"] * len(eligible)))]
        )
        inner_val = np.array([i for i in tr if groups[i] in inner_groups], dtype=int)
        inner_train = np.array([i for i in tr if groups[i] not in inner_groups], dtype=int)
        assert set(groups[inner_train]).isdisjoint(set(groups[inner_val]))
        assert set(va).isdisjoint(inner_train) and set(va).isdisjoint(inner_val)
        assert not anchors[inner_val].any()
        for method, distance in distances.items():
            trials = []
            for k in cfg["k_grid"]:
                for weight in cfg["weight_grid"]:
                    p, _, _ = predict(distance, inner_val, inner_train, y, k, weight)
                    s = score(h[inner_val], p)
                    trials.append((s, k, weight))
                    inner_scores.append(
                        dict(
                            endpoint=endpoint,
                            fold=fold,
                            method=method,
                            k=k,
                            weight=weight,
                            ST_RAE=s,
                            inner_train_n=len(inner_train),
                            inner_validation_n=len(inner_val),
                        )
                    )
            best, k, weight = min(trials)
            prediction, near, dist = predict(distance, va, tr, y, k, weight)
            assert np.isfinite(prediction).all()
            part = (
                h[va]
                .select(
                    "compound_id",
                    "endpoint",
                    "fold",
                    pl.lit(cfg["cv_run"]).alias("split_id"),
                    "is_neighbor",
                    pl.col("value").alias("y_true"),
                    "conf_low",
                    "conf_high",
                )
                .with_columns(pl.Series("prediction", prediction), pl.lit(method).alias("method"))
            )
            outputs[method].append(part)
            selections.append(
                dict(
                    endpoint=endpoint,
                    fold=fold,
                    method=method,
                    k=k,
                    weight=weight,
                    inner_ST_RAE=best,
                    outer_reference_n=len(tr),
                    outer_validation_n=len(va),
                    inner_train_n=len(inner_train),
                    inner_validation_n=len(inner_val),
                    reference_model=(
                        "Keyed reference-membership.parquet plus declared immutable "
                        "input features and labels"
                    ),
                )
            )
            for row in tr:
                references.append(
                    dict(
                        endpoint=endpoint,
                        outer_fold=fold,
                        method=method,
                        compound_id=ids[row],
                        role="inner_validation" if row in inner_val else "inner_training",
                        is_anchor=bool(anchors[row]),
                        group=str(groups[row]),
                    )
                )
            for j, row in enumerate(va):
                for rank, neighbor in enumerate(near[j, :5]):
                    neighbors.append(
                        dict(
                            endpoint=endpoint,
                            fold=fold,
                            method=method,
                            compound_id=ids[row],
                            neighbor_compound_id=ids[neighbor],
                            rank=rank + 1,
                            distance=float(dist[j, rank]),
                        )
                    )
            print(
                json.dumps(
                    dict(
                        endpoint=endpoint,
                        fold=fold,
                        method=method,
                        k=k,
                        weight=weight,
                        inner_score=best,
                        outer_score=score(h[va], prediction),
                    )
                ),
                flush=True,
            )
for method, parts in outputs.items():
    oof = pl.concat(parts).sort("compound_id", "endpoint", "fold")
    assert oof.height == 6384 and oof.select("compound_id", "endpoint", "fold").n_unique() == 6384
    name = "morgan-oof.parquet" if method == "morgan_tanimoto_knn" else "embedding-oof.parquet"
    oof.write_parquet(out / name)
    for (endpoint,), h in oof.group_by("endpoint"):
        for subset in ["all", "neighbors"]:
            g = h if subset == "all" else h.filter("is_neighbor")
            metrics.append(
                dict(
                    method=method,
                    endpoint=endpoint,
                    subset=subset,
                    n=g.height,
                    ST_RAE=score(g.rename({"y_true": "value"}), g["prediction"].to_numpy()),
                    MAE=float((g["prediction"] - g["y_true"]).abs().mean()),
                    Spearman=float(spearmanr(g["y_true"], g["prediction"]).statistic),
                )
            )
pl.DataFrame(selections).write_json(out / "selected-parameters.json")
pl.DataFrame(inner_scores).write_csv(out / "inner-selection-scores.csv")
pl.DataFrame(references).write_parquet(out / "reference-membership.parquet")
pl.DataFrame(neighbors).write_parquet(out / "nearest-neighbor-audit.parquet")
m = pl.DataFrame(metrics)
m.write_csv(out / "metrics.csv")
macro = (
    m.group_by("method", "subset")
    .agg(pl.col("ST_RAE").mean(), pl.col("MAE").mean(), pl.col("Spearman").mean())
    .sort("method", "subset")
)
macro.write_csv(out / "macro.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            macro=macro.to_dicts(),
            selection="Inner group holdout only; no outer label selection",
        ),
        indent=2,
    )
)
print(macro)
