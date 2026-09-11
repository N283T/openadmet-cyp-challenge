"""Local-label ablations and observed-primary diagnostic under corrected parent CV."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(n):
    return pl.read_parquet(next(inp.glob("*-" + n + ".parquet")))


a, labels, primary = read("parent-assignments"), read("measurements"), read("targets")
base = read("oof").filter((pl.col("mode") == "parent_retained") & (pl.col("features") == "full"))
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = {
    c: gen.GetFingerprint(Chem.MolFromSmiles(s))
    for c, s in read("compounds").select("compound_id", "standardized_smiles").iter_rows()
}
spec = importlib.util.spec_from_file_location(
    "official", next(inp.glob("*-custom_scoring_functions.py"))
)
official = importlib.util.module_from_spec(spec)
spec.loader.exec_module(official)
rows = []
diagnostics = []
preds = []
primaryrows = []
for ep in sorted(base["endpoint"].unique()):
    enzyme = ep.split("_")[0]
    pr = primary.filter(pl.col("head").str.starts_with(enzyme + "@")).select(
        "compound_id", pl.col("value").alias("primary")
    )
    f = (
        labels.filter(pl.col("endpoint") == ep)
        .join(a, on="compound_id", validate="1:1")
        .join(pr, on="compound_id", how="left", validate="1:1")
        .sort("compound_id")
    )
    for fold in range(5):
        train = f.filter((pl.col("fold") != fold) | pl.col("is_anchor"))
        val = base.filter((pl.col("endpoint") == ep) & (pl.col("fold") == fold)).sort("compound_id")
        assert set(train["compound_id"]).isdisjoint(val["compound_id"])
        qmaps = {}
        for pool in ["all_train", "anchors"]:
            t = train if pool == "all_train" else train.filter(pl.col("is_anchor"))
            values = t["value"].to_numpy()
            tfps = [fps[c] for c in t["compound_id"]]
            sims = np.array(
                [DataStructs.BulkTanimotoSimilarity(fps[c], tfps) for c in val["compound_id"]]
            )
            order = np.argsort(-sims, axis=1, kind="stable")
            nearest = values[order[:, 0]]
            maxsim = np.take_along_axis(sims, order[:, :1], axis=1)[:, 0]
            if pool == "all_train":
                for i, r in enumerate(val.iter_rows(named=True)):
                    js = order[i, : min(5, len(t))]
                    diagnostics.append(
                        dict(
                            compound_id=r["compound_id"],
                            endpoint=ep,
                            fold=fold,
                            is_neighbor=r["is_neighbor"],
                            similarity=float(maxsim[i]),
                            nearest_y=float(nearest[i]),
                            y_true=r["y_true"],
                            prediction=r["prediction"],
                            abs_error=abs(r["prediction"] - r["y_true"]),
                            nearest_gap=abs(nearest[i] - r["y_true"]),
                            local_label_sd=float(values[js].std()),
                        )
                    )
            for k in cfg["k"]:
                js = order[:, : min(k, len(t))]
                w = np.take_along_axis(sims, js, axis=1) ** 4
                p = (values[js] * w).sum(axis=1) / w.sum(axis=1)
                for strength in cfg["strength"]:
                    name = f"{pool}:k{k}:w{strength}"
                    qmaps[name] = (1 - strength) * val["prediction"].to_numpy() + strength * p
                    if strength == 0.5:
                        gated = np.clip((maxsim - 0.3) / 0.3, 0, 1) * strength
                        qmaps[name + ":similarity_gated"] = (1 - gated) * val[
                            "prediction"
                        ].to_numpy() + gated * p
        qmaps["baseline"] = val["prediction"].to_numpy()
        for name, p in qmaps.items():
            preds.append(val.with_columns(pl.lit(name).alias("method"), pl.Series("prediction", p)))
        t = train.filter(pl.col("primary").is_not_null())
        v = val.join(
            pr, on="compound_id", how="left", validate="1:1", maintain_order="left"
        ).filter(pl.col("primary").is_not_null())
        iso = IsotonicRegression(increasing=False, out_of_bounds="clip").fit(
            t["primary"].to_numpy(), t["value"].to_numpy()
        )
        prediction = iso.predict(v["primary"].to_numpy())
        for name, p in [
            ("observed_primary_diagnostic", prediction),
            ("baseline_matched_primary", v["prediction"].to_numpy()),
        ]:
            primaryrows.append(
                v.drop("primary").with_columns(
                    pl.lit(name).alias("method"), pl.Series("prediction", p)
                )
            )
allpred = pl.concat(preds + primaryrows)
for (ep, method), f in allpred.group_by("endpoint", "method"):
    for subset in ["all", "neighbors"]:
        g = f if subset == "all" else f.filter(pl.col("is_neighbor"))
        y, p, lo, hi = [g[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]]
        rows.append(
            dict(
                endpoint=ep,
                method=method,
                subset=subset,
                n=g.height,
                ST_RAE=float(
                    official.rae_soft_threshold_absolute_error(
                        y, p, y_true_lower=lo, y_true_upper=hi
                    )
                ),
                MAE=float(np.abs(y - p).mean()),
                Spearman=float(spearmanr(y, p).statistic),
            )
        )
m = pl.DataFrame(rows)
m.write_csv(out / "metrics.csv")
m.group_by("method", "subset").agg(pl.col("ST_RAE", "MAE", "Spearman").mean()).sort(
    "subset", "ST_RAE"
).write_csv(out / "macro.csv")
d = pl.DataFrame(diagnostics)
d.write_parquet(out / "diagnostics.parquet")
d.with_columns(
    pl.when(pl.col("similarity") < 0.4)
    .then(pl.lit("<0.4"))
    .when(pl.col("similarity") < 0.6)
    .then(pl.lit("0.4-0.6"))
    .otherwise(pl.lit(">=0.6"))
    .alias("similarity_bin")
).group_by("endpoint", "is_neighbor", "similarity_bin").agg(
    pl.len().alias("n"), pl.col("abs_error", "nearest_gap", "local_label_sd").mean()
).write_csv(out / "error-strata.csv")
allpred.write_parquet(out / "predictions.parquet")
write_json(
    out / "summary.json",
    dict(
        config=cfg,
        scope=(
            "Local blends use outer-training labels. Observed-primary diagnostic is "
            "not deployable; comparisons use matched primary-available rows."
        ),
    ),
)
