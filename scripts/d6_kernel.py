"""2D6 count-fingerprint kernel SVR, selected on grouped inner ST-RAE."""

import importlib.util
import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator as fg
from scipy.stats import spearmanr
from sklearn.svm import SVR
from threadpoolctl import threadpool_limits

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


spec = importlib.util.spec_from_file_location("official", path(5))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
threadpool_limits(3)
comp, labels, assignments, emb, old = [pl.read_parquet(path(i)) for i in range(5)]
ep = "CYP2D6_pIC50_direct_inhibition"
h = (
    labels.filter(pl.col("endpoint") == ep)
    .join(assignments, on="compound_id", validate="1:1")
    .join(comp, on="compound_id", validate="1:1")
    .join(emb, on="compound_id", validate="1:1")
    .sort("compound_id")
)
y = h["value"].to_numpy()
anc, groups = h["is_anchor"].to_numpy(), h["group"].to_numpy()
mols = [Chem.MolFromSmiles(s) for s in h["standardized_smiles"]]
assert all(m is not None for m in mols)
generators = {
    "morgan_count": fg.GetMorganGenerator(radius=2, fpSize=4096),
    "feature_morgan": fg.GetMorganGenerator(
        radius=2, fpSize=4096, atomInvariantsGenerator=fg.GetMorganFeatureAtomInvGen()
    ),
    "atom_pair": fg.GetAtomPairGenerator(fpSize=4096),
    "torsion": fg.GetTopologicalTorsionGenerator(fpSize=4096),
}
kernels = {}
for name, gen in generators.items():
    fps = [gen.GetSparseCountFingerprint(m) for m in mols]
    kernels[name] = np.asarray([DataStructs.BulkTanimotoSimilarity(f, fps) for f in fps])
kernels["fingerprint_mix"] = sum(kernels.values()) / len(kernels)
x = h.select(sorted(c for c in h.columns if c.startswith("emb_"))).to_numpy()
x = x / np.maximum(np.linalg.norm(x, axis=1, keepdims=True), 1e-12)
kernels["fp_emb_mix"] = 0.5 * kernels["fingerprint_mix"] + 0.5 * (x @ x.T)
np.savez_compressed(out / "kernels.npz", **kernels)
h.select("compound_id").write_parquet(out / "kernel-ids.parquet")


def score(frame, pred):
    return float(
        scorer.rae_soft_threshold_absolute_error(
            frame["value"].to_numpy(),
            pred,
            y_true_lower=frame["conf_low"].to_numpy(),
            y_true_upper=frame["conf_high"].to_numpy(),
        )
    )


rows, history, memberships = [], [], []
for fold in range(5):
    val = (h["fold"].to_numpy() == fold) & ~anc
    tr, va = np.flatnonzero(~val), np.flatnonzero(val)
    eligible = sorted(set(groups[tr]) - set(groups[tr[anc[tr]]]))
    inner_groups = set(
        np.random.default_rng(42 + fold).permutation(eligible)[: max(1, int(0.2 * len(eligible)))]
    )
    it = np.array([i for i in tr if groups[i] not in inner_groups])
    iv = np.array([i for i in tr if groups[i] in inner_groups])
    assert set(groups[it]).isdisjoint(groups[iv])
    memberships.append(
        h.select("compound_id").with_columns(
            pl.lit(fold).alias("outer_fold"),
            pl.Series("outer_validation", val),
            pl.Series("inner_validation", np.isin(np.arange(len(h)), iv)),
        )
    )
    for name, kernel in kernels.items():
        trials = []
        for c in cfg["C"]:
            for epsilon in cfg["epsilon"]:
                model = SVR(kernel="precomputed", C=c, epsilon=epsilon)
                model.fit(kernel[np.ix_(it, it)], y[it])
                loss = score(h[iv], model.predict(kernel[np.ix_(iv, it)]))
                trials.append((loss, c, epsilon))
                history.append(dict(method=name, fold=fold, C=c, epsilon=epsilon, ST_RAE=loss))
        loss, c, epsilon = min(trials)
        model = SVR(kernel="precomputed", C=c, epsilon=epsilon)
        model.fit(kernel[np.ix_(tr, tr)], y[tr])
        pred = model.predict(kernel[np.ix_(va, tr)])
        joblib.dump(
            dict(model=model, reference_ids=h[tr]["compound_id"].to_list(), inner_ST_RAE=loss),
            out / f"{name}-{fold}.joblib",
        )
        rows.append(
            h[va]
            .select(
                "compound_id",
                "endpoint",
                "fold",
                "is_neighbor",
                "conf_low",
                "conf_high",
                pl.col("value").alias("y_true"),
            )
            .with_columns(
                pl.lit(name).alias("method"),
                pl.lit(cfg["split_run"]).alias("split_id"),
                pl.Series("prediction", pred),
            )
        )
        print(name, fold, c, epsilon, score(h[va], pred), flush=True)
oof = pl.concat(rows)
oof.write_parquet(out / "d6-kernel-oof.parquet")
pl.DataFrame(history).write_csv(out / "inner-trials.csv")
pl.concat(memberships).write_parquet(out / "memberships.parquet")
base = old.filter(
    (pl.col("endpoint") == ep) & (pl.col("method") == "caruana_weight_crossfit")
).sort("compound_id")
metrics = []
for method in oof["method"].unique().sort():
    p = oof.filter(pl.col("method") == method).sort("compound_id")
    for col in ["compound_id", "fold", "y_true", "conf_low", "conf_high"]:
        assert p[col].equals(base[col]), col
    for weight in [0, 0.25, 0.5, 1]:
        pred = (1 - weight) * base["prediction"].to_numpy() + weight * p["prediction"].to_numpy()
        for subset in ["all", "neighbors"]:
            mask = np.ones(len(p), dtype=bool) if subset == "all" else p["is_neighbor"].to_numpy()
            f = p.filter(pl.Series(mask)).rename({"y_true": "value"})
            metrics.append(
                dict(
                    method=method,
                    weight=weight,
                    subset=subset,
                    ST_RAE=score(f, pred[mask]),
                    Spearman=float(spearmanr(f["value"], pred[mask]).statistic),
                )
            )
pl.DataFrame(metrics).write_csv(out / "metrics.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            n=len(h),
            oof_per_method=len(base),
            selection="grouped inner 20%; outer validation never selected C/epsilon",
        ),
        indent=2,
    )
)
