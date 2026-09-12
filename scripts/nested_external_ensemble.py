"""Nested Caruana selection on inner groups, with fresh outer-training refits."""

import importlib.util
import json
import os
from collections import Counter
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(key):
    return pl.read_parquet(inp / cfg[key])


a, labels, compounds = read("assignments"), read("measurements"), read("compounds")
a = a.sort("compound_id")
ids = a.select("compound_id")
class_emb, mixed_emb, tree, nested = (
    read("class_emb"),
    read("mixed_emb"),
    read("tree"),
    read("nested"),
)


def matrix(table):
    cols = sorted(c for c in table.columns if c != "compound_id")
    x = (
        ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
        .select(cols)
        .to_numpy()
        .astype(np.float32)
    )
    assert np.isfinite(x).all()
    return x, cols


xc, classcols = matrix(class_emb)
xm, mixedcols = matrix(mixed_emb)
xt, treecols = matrix(tree)
assert xc.shape[1] == xm.shape[1] == 600 and xt.shape[1] == 20
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
smiles = ids.join(compounds, on="compound_id", validate="1:1", maintain_order="left")[
    "standardized_smiles"
]
xf = np.stack([fpgen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in smiles])
index = {c: i for i, c in enumerate(a["compound_id"])}
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
rows, summary, memberships = [], [], []
for outer in range(5):
    xn = (
        ids.join(
            nested.filter(pl.col("outer_fold") == outer),
            on="compound_id",
            validate="1:1",
            maintain_order="left",
        )
        .select([f"aux_{i}" for i in range(4)])
        .to_numpy()
        .astype(np.float32)
    )
    views = {
        "morgan": xf,
        "logfc": xn,
        "class": xc,
        "mixed": xm,
        "tree": xt,
        "class_logfc": np.column_stack([xc, xn]),
        "mixed_logfc": np.column_stack([xm, xn]),
        "tree_logfc": np.column_stack([xt, xn]),
    }
    for ep in sorted(
        labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))[
            "endpoint"
        ].unique()
    ):
        t = (
            labels.filter(pl.col("endpoint") == ep)
            .join(a, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        ix = np.array([index[c] for c in t["compound_id"]])
        ov = (t["fold"].to_numpy() == outer) & ~t["is_anchor"].to_numpy()
        ot = ~ov
        anchor_groups = set(a.filter(pl.col("is_anchor"))["group"])
        groups = sorted(set(t.filter(pl.Series(ot))["group"]) - anchor_groups)
        rng = np.random.default_rng(cfg["seed"] + outer)
        vg = set(rng.permutation(groups)[: max(1, round(len(groups) * 0.2))])
        iv = ot & np.array([g in vg for g in t["group"]])
        it = ot & ~iv
        assert not (iv & t["is_anchor"].to_numpy()).any()
        assert set(t.filter(pl.Series(it))["group"]).isdisjoint(
            set(t.filter(pl.Series(iv))["group"])
        )
        d = out / f"fold-{outer}" / ep
        d.mkdir(parents=True)
        member = t.select("compound_id", "endpoint", "fold", "group", "is_anchor").with_columns(
            pl.lit(outer).alias("outer_fold"),
            pl.lit(cfg["split_run"]).alias("split_id"),
            pl.Series(
                "role",
                np.where(ov, "outer_validation", np.where(iv, "inner_validation", "inner_train")),
            ),
        )
        member.write_parquet(d / "membership.parquet")
        memberships.append(member)
        ytrain = t.filter(pl.Series(it))["value"].to_numpy()
        valid = t.filter(pl.Series(iv))
        yv, lo, hi = (valid[c].to_numpy() for c in ["value", "conf_low", "conf_high"])
        candidates, predictions = [], []
        for view, x in views.items():
            for objective, alpha in [
                ("regression_l1", None),
                ("quantile", 0.35),
                ("quantile", 0.65),
            ]:
                params = dict(
                    objective=objective,
                    n_estimators=300,
                    num_leaves=15,
                    min_child_samples=20,
                    learning_rate=0.05,
                    reg_lambda=1.0,
                    n_jobs=3,
                    random_state=cfg["seed"] + outer,
                    verbosity=-1,
                    deterministic=True,
                    force_col_wise=True,
                )
                if alpha is not None:
                    params["alpha"] = alpha
                name = f"{view}_{objective}_{alpha}"
                model = lgb.LGBMRegressor(**params)
                model.fit(x[ix[it]], ytrain)
                pv = model.predict(x[ix[iv]])
                predictions.append(pv)
                candidates.append(dict(name=name, view=view, params=params, features=x.shape[1]))
        p = np.stack(predictions)
        assert p.shape == (24, int(iv.sum())) and np.isfinite(p).all()
        inner_predictions = valid.select(
            "compound_id",
            "endpoint",
            "fold",
            "conf_low",
            "conf_high",
            pl.col("value").alias("y_true"),
        ).with_columns(*[pl.Series(c["name"], p[i]) for i, c in enumerate(candidates)])
        inner_predictions.write_parquet(d / "inner-predictions.parquet")
        counts = Counter()
        total = np.zeros(len(yv))
        history = []
        best_score = float("inf")
        best_weights = None
        for step in range(1, cfg["steps"] + 1):
            losses = [
                float(
                    scorer.rae_soft_threshold_absolute_error(
                        yv, (total + pv) / step, y_true_lower=lo, y_true_upper=hi
                    )
                )
                for pv in p
            ]
            choice = int(np.argmin(losses))
            counts[choice] += 1
            total += p[choice]
            score = losses[choice]
            history.append(
                dict(
                    step=step,
                    selected=candidates[choice]["name"],
                    score=score,
                    candidate_scores=losses,
                )
            )
            if score < best_score:
                best_score, best_weights, best_step = score, dict(counts), step
        (d / "caruana-history.json").write_text(json.dumps(history, indent=2))
        outer_prediction = np.zeros(int(ov.sum()))
        selected = []
        for i, count in best_weights.items():
            candidate = candidates[i]
            x = views[candidate["view"]]
            model = lgb.LGBMRegressor(**candidate["params"])
            model.fit(x[ix[ot]], t.filter(pl.Series(ot))["value"].to_numpy())
            prediction = model.predict(x[ix[ov]])
            weight = count / best_step
            outer_prediction += weight * prediction
            model.booster_.save_model(str(d / (candidate["name"] + ".txt")))
            selected.append({**candidate, "count": count, "weight": weight})
        assert np.isfinite(outer_prediction).all() and np.isclose(
            sum(s["weight"] for s in selected), 1
        )
        info = dict(
            endpoint=ep,
            fold=outer,
            best_step=best_step,
            inner_score=best_score,
            inner_train=int(it.sum()),
            inner_valid=int(iv.sum()),
            outer_train=int(ot.sum()),
            outer_valid=int(ov.sum()),
            selected=selected,
        )
        (d / "selected.json").write_text(json.dumps(info, indent=2))
        (d / "candidates.json").write_text(json.dumps(candidates, indent=2))
        summary.append(info)
        rows.append(
            t.filter(pl.Series(ov))
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
                pl.Series("prediction", outer_prediction),
                pl.lit(cfg["method"]).alias("method"),
                pl.lit(cfg["split_run"]).alias("split_id"),
            )
        )
        pl.concat(rows).write_parquet(out / "external-oof-progress.parquet")
        print(json.dumps({k: v for k, v in info.items() if k != "selected"}), flush=True)
result = pl.concat(rows)
assert result.height == result.select("compound_id", "endpoint", "fold").n_unique() == 6384
result.write_parquet(out / "external-oof.parquet")
pl.concat(memberships).write_parquet(out / "membership.parquet")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            selections=summary,
            feature_columns={"class": classcols, "mixed": mixedcols, "tree": treecols},
        ),
        indent=2,
    )
)
