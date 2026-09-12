"""Structural coverage, prevalence-controlled label smoothness and observed assay regimes."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import matthews_corrcoef, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


comp, measure, assign, test, obs, oof = [pl.read_parquet(path(i)) for i in range(6)]
assign = assign.filter(pl.col("strategy") == "umap")
comp = comp.sort("compound_id").join(
    assign.select("compound_id", "fold"),
    on="compound_id",
    how="left",
    validate="1:1",
    maintain_order="left",
)
ids = comp["compound_id"].to_numpy()
folds = comp["fold"].fill_null(-1).to_numpy()
lookup = {c: i for i, c in enumerate(ids)}
assert len(comp) == 6896 and len(test) == 750
ends = ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]
y = np.full((len(comp), 2), np.nan)
for j, ep in enumerate(ends):
    rows = measure.filter(pl.col("endpoint") == ep)
    for r in rows.iter_rows(named=True):
        y[lookup[r["compound_id"]], j] = r["value"]
assert np.isfinite(y).sum() == 5081
fg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = [fg.GetFingerprint(Chem.MolFromSmiles(s)) for s in comp["standardized_smiles"]]
shared = np.isfinite(y).all(1)
records = []
permutation = []
coverage = []
rng = np.random.default_rng(20260912)
for cohort in ["endpoint_all", "shared_only"]:
    for j, ep in enumerate(ends):
        available = np.isfinite(y[:, j]) if cohort == "endpoint_all" else shared
        nearest = np.full(len(comp), -1)
        sim = np.zeros(len(comp))
        for fold in [-1, 0, 1, 2, 3, 4]:
            query = np.flatnonzero(folds == fold)
            pool = np.flatnonzero(available & (folds >= 0) & (folds != fold))
            poolfp = [fps[i] for i in pool]
            for ix in query:
                similarities = np.array(DataStructs.BulkTanimotoSimilarity(fps[ix], poolfp))
                nearest[ix] = pool[similarities.argmax()]
                sim[ix] = similarities.max()
        assert np.all(nearest >= 0)
        bins = np.digitize(sim, [0.3, 0.5, 0.7, 0.9])
        querylabel = available & (folds >= 0)
        disagreement = y[:, j] != y[nearest, j]
        t = pl.DataFrame(
            dict(
                compound_id=ids,
                endpoint=ep,
                cohort=cohort,
                fold=folds,
                y_true=y[:, j],
                label_available=np.isfinite(y[:, j]),
                shared_labelled=shared,
                nearest_compound_id=ids[nearest],
                nearest_label=y[nearest, j],
                similarity=sim,
                similarity_bin=bins,
                disagreement=disagreement,
            )
        )
        records.append(t)
        for role, mask in [
            ("labelled", querylabel),
            ("missing_label", (folds >= 0) & ~np.isfinite(y[:, j])),
            ("test", folds < 0),
            ("labelled_positive", querylabel & (y[:, j] == 1)),
            ("labelled_negative", querylabel & (y[:, j] == 0)),
        ]:
            coverage.append(
                dict(
                    endpoint=ep,
                    cohort=cohort,
                    role=role,
                    n=int(mask.sum()),
                    median_similarity=float(np.median(sim[mask])),
                    mean_similarity=float(sim[mask].mean()),
                    fraction_ge05=float((sim[mask] >= 0.5).mean()),
                    fraction_ge07=float((sim[mask] >= 0.7).mean()),
                )
            )
        masks = {"all": querylabel}
        masks.update({f"sim_bin{k}": querylabel & (bins == k) for k in range(5)})
        null = {k: [] for k, m in masks.items() if m.any()}
        for _rep in range(500):
            yp = y[:, j].copy()
            for fold in range(5):
                ix = np.flatnonzero(available & (folds == fold))
                yp[ix] = rng.permutation(yp[ix])
            d = yp != yp[nearest]
            for k in null:
                null[k].append(float(d[masks[k]].mean()))
        for k, dist in null.items():
            m = masks[k]
            observed = float(disagreement[m].mean())
            permutation.append(
                dict(
                    endpoint=ep,
                    cohort=cohort,
                    stratum=k,
                    n=int(m.sum()),
                    positive_rate=float(y[m, j].mean()),
                    observed_disagreement=observed,
                    null_mean=float(np.mean(dist)),
                    null_low=float(np.quantile(dist, 0.025)),
                    null_high=float(np.quantile(dist, 0.975)),
                    agreement_gain=float(np.mean(dist) - observed),
                )
            )
        print(ep, cohort, "done", flush=True)
near = pl.concat(records)
near.write_parquet(out / "neighbors.parquet")
pl.DataFrame(coverage).write_csv(out / "structural-coverage.csv")
pl.DataFrame(permutation).write_csv(out / "prevalence-controlled-disagreement.csv")
# Official calls are never replaced. Bounds define descriptive observed-potency regimes only.
routes = []
for j, ep in enumerate(ends):
    enzyme = ep.split("_")[0]
    t = comp.select("compound_id", "fold").with_columns(
        pl.Series("y_true", y[:, j]), pl.lit(ep).alias("endpoint")
    )
    for kind in ["direct", "tdi"]:
        q = obs.filter((pl.col("enzyme") == enzyme) & (pl.col("kind") == kind)).select(
            "compound_id",
            *[pl.col(c).alias(kind + "_" + c) for c in ["value", "conf_low", "conf_high"]],
        )
        t = t.join(q, on="compound_id", how="left", validate="1:1")
    t = t.with_columns(
        (pl.col("tdi_value") - pl.col("direct_value")).alias("shift"),
        (pl.col("tdi_conf_low") - pl.col("direct_conf_high")).alias("shift_lower"),
        (pl.col("tdi_conf_high") - pl.col("direct_conf_low")).alias("shift_upper"),
    )
    t = t.with_columns(
        pl.when(pl.col("direct_value").is_not_null() & pl.col("tdi_value").is_not_null())
        .then(pl.lit("both"))
        .when(pl.col("direct_value").is_not_null())
        .then(pl.lit("direct_only"))
        .when(pl.col("tdi_value").is_not_null())
        .then(pl.lit("tdi_only"))
        .otherwise(pl.lit("neither"))
        .alias("availability"),
        pl.when(pl.col("shift").is_null())
        .then(pl.lit("unavailable"))
        .when(pl.col("shift") < 0)
        .then(pl.lit("negative"))
        .when(pl.col("shift") < 0.3)
        .then(pl.lit("0_to_0.3"))
        .when(pl.col("shift") < 1)
        .then(pl.lit("0.3_to_1"))
        .otherwise(pl.lit("ge1"))
        .alias("shift_regime"),
        pl.when(pl.col("shift_lower").is_null() | pl.col("shift_upper").is_null())
        .then(pl.lit("unavailable"))
        .when(pl.col("shift_lower") > 0)
        .then(pl.lit("positive_nonoverlap"))
        .when(pl.col("shift_upper") < 0)
        .then(pl.lit("negative_nonoverlap"))
        .otherwise(pl.lit("overlap_zero"))
        .alias("bounds_regime"),
    )
    routes.append(t)
route = pl.concat(routes)
route.write_parquet(out / "assay-route.parquet")
rows = []
for ep in ends:
    selected = "all_plain_pfn_graph_mean" if "2D6" in ep else "lgbm_balanced_mean"
    pred = oof.filter((pl.col("endpoint") == ep) & (pl.col("method") == selected))
    h = (
        route.filter(pl.col("endpoint") == ep)
        .join(
            pred.select("compound_id", "fold", "prediction", "threshold"),
            on=["compound_id", "fold"],
            how="inner",
            validate="1:1",
        )
        .join(
            near.filter((pl.col("endpoint") == ep) & (pl.col("cohort") == "endpoint_all")).select(
                "compound_id", "similarity_bin"
            ),
            on="compound_id",
            validate="1:1",
        )
    )
    assert len(h) == (1497 if "2D6" in ep else 3584)
    for dimension in [
        "availability",
        "shift_regime",
        "bounds_regime",
        "similarity_bin",
        "shared_labelled",
    ]:
        if dimension == "shared_labelled":
            h = h.join(
                comp.select("compound_id").with_columns(pl.Series("shared_labelled", shared)),
                on="compound_id",
                validate="1:1",
            )
        for (stratum,), g in h.group_by(dimension):
            yy = g["y_true"].to_numpy()
            pp = g["prediction"].to_numpy()
            call = pp >= g["threshold"].to_numpy()
            rows.append(
                dict(
                    endpoint=ep,
                    method=selected,
                    dimension=dimension,
                    stratum=str(stratum),
                    n=len(g),
                    positive_rate=float(yy.mean()),
                    error_rate=float((call != yy).mean()),
                    brier=float(((pp - yy) ** 2).mean()),
                    MCC=float(matthews_corrcoef(yy, call)),
                    AUC=float(roc_auc_score(yy, pp)) if len(np.unique(yy)) == 2 else None,
                )
            )
pl.DataFrame(rows).write_csv(out / "selected-ensemble-regime-errors.csv")
counts = []
for dimension in ["availability", "shift_regime", "bounds_regime"]:
    g = (
        route.filter(pl.col("fold").is_not_null())
        .with_columns(pl.col("y_true").is_finite().alias("label_available"))
        .group_by("endpoint", dimension, "label_available")
        .agg(
            pl.len().alias("n"),
            pl.col("y_true").filter(pl.col("y_true").is_finite()).mean().alias("positive_rate"),
        )
        .rename({dimension: "stratum"})
        .with_columns(pl.lit(dimension).alias("dimension"))
    )
    counts.append(g)
pl.concat(counts).write_csv(out / "label-route-counts.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            class_labels=5081,
            shared_labelled_compounds=int(shared.sum()),
            compounds=6896,
            test750=True,
            permutations=500,
            permutation="Shuffle official labels within UMAP fold and cohort "
            "preserving per-fold prevalence; fixed graph",
            fingerprint="Morgan radius2 2048 bits",
            heldout_exclusion="All query-fold compounds excluded from labelled neighbor pool",
            diagnostic_only="Observed assay values and bounds used only for retrospective "
            "stratification; no labels replaced or deploy features produced",
        ),
        indent=2,
    )
)
