"""Does adding frozen CYPMol D6 vectors to the deployed L4 rank-mean D6 score help?

CYPMol (isoform-conditioned, trained on public inhibition data) is used as a plain feature set for
the same pool-trained classifier family as `pool_lgbm`: official labels plus primary-screen
pseudo-negatives, fitted on the other outer folds of the frozen umap/scaffold splits. The held-fold
score is rank-averaged with the two L4 components and compared with L4 at the historical call count
under Metric O, Metric M and the pseudo analog-expansion pool, with paired grouped bootstraps.
Nothing is uploaded and no submission is prepared.
"""

import json
import os
from pathlib import Path

import evaluate_cyp_d6_official_metric as om
import numpy as np
import polars as pl
from lightgbm import LGBMClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler


def percentile_by_fold(score, fold):
    out = np.empty(len(score))
    for f in np.unique(fold):
        rows = np.flatnonzero(fold == f)
        order = np.argsort(score[rows], kind="stable")
        ranks = np.empty(len(rows))
        ranks[order] = np.arange(len(rows))
        out[rows] = (ranks + 0.5) / len(rows)
    return out


def fit_lgbm(x, rows, y, predict, cfg):
    scores = []
    for s in range(cfg["lgbm_seeds"]):
        model = LGBMClassifier(
            class_weight="balanced",
            **cfg["lgbm"],
            verbosity=-1,
            n_jobs=cfg["threads"],
            random_state=cfg["seed"] + s,
            deterministic=True,
            force_col_wise=True,
        )
        model.fit(x[rows], y[rows].astype(int))
        scores.append(np.asarray(model.predict_proba(x[predict]))[:, 1])
    return np.mean(scores, axis=0)


def fit_lr(x, rows, y, groups, predict, cfg):
    splitter = GroupKFold(n_splits=cfg["inner_folds"])
    tuning = {}
    for c in cfg["lr_C"]:
        inner = np.full(len(rows), np.nan)
        for fit, valid in splitter.split(rows, groups=groups[rows]):
            scaler = StandardScaler().fit(x[rows[fit]])
            model = LogisticRegression(
                C=c, class_weight="balanced", max_iter=cfg["lr_max_iter"], random_state=cfg["seed"]
            ).fit(scaler.transform(x[rows[fit]]), y[rows[fit]].astype(int))
            inner[valid] = model.decision_function(scaler.transform(x[rows[valid]]))
        tuning[c] = average_precision_score(y[rows].astype(int), inner)
    best = max(cfg["lr_C"], key=lambda c: tuning[c])
    scaler = StandardScaler().fit(x[rows])
    model = LogisticRegression(
        C=best, class_weight="balanced", max_iter=cfg["lr_max_iter"], random_state=cfg["seed"]
    ).fit(scaler.transform(x[rows]), y[rows].astype(int))
    return model.decision_function(scaler.transform(x[predict]))


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    oof = pl.read_parquet(inp / "testlike-pool-oof.parquet")
    assignments = pl.read_parquet(inp / "assignments.parquet")
    embeddings = pl.read_parquet(inp / "embeddings.parquet")
    analog = set(pl.read_parquet(inp / "analog-pairs.parquet")["analog_id"].to_list())
    feature_sets = {
        name: [c for c in embeddings.columns if c.startswith(tuple(prefixes))]
        for name, prefixes in cfg["feature_sets"].items()
    }
    rate = cfg["reference_calls"] / 750
    standalone, results, count_rows = [], [], []
    for split in sorted(oof["split_id"].unique().to_list()):
        rows = oof.filter(pl.col("split_id") == split)
        wide = (
            rows.filter(pl.col("model") == "pool_lgbm")
            .select("compound_id", "fold", "pool_role", "y_official", "y_testlike")
            .join(
                rows.filter(pl.col("model") == "legacy_sep15_refit").select(
                    "compound_id", pl.col("score").alias("legacy")
                ),
                on="compound_id",
            )
            .join(
                rows.filter(pl.col("model") == "pool_lgbm").select(
                    "compound_id", pl.col("score").alias("pool")
                ),
                on="compound_id",
            )
            .join(
                assignments.filter(pl.col("strategy") == split).select("compound_id", "group"),
                on="compound_id",
            )
            .join(embeddings, on="compound_id", how="left")
            .sort("compound_id")
        )
        assert wide.height == 6146 and wide["group"].null_count() == 0
        fold, groups = wide["fold"].to_numpy(), wide["group"].to_numpy().astype(str)
        y_official = wide["y_official"].fill_null(np.nan).to_numpy()
        y_testlike = wide["y_testlike"].fill_null(np.nan).to_numpy()
        role = wide["pool_role"].to_numpy()
        legacy_rank = percentile_by_fold(wide["legacy"].to_numpy(), fold)
        pool_rank = percentile_by_fold(wide["pool"].to_numpy(), fold)
        scores = {"L4": 0.5 * (legacy_rank + pool_rank)}
        for fname, columns in feature_sets.items():
            x = wide.select(columns).to_numpy().astype(np.float32)
            for kind in ("lgbm", "lr"):
                raw = np.full(len(fold), np.nan)
                for f in np.unique(fold):
                    held = np.flatnonzero(fold == f)
                    pool_rows = np.flatnonzero((fold != f) & np.isfinite(y_testlike))
                    if kind == "lgbm":
                        raw[held] = fit_lgbm(x, pool_rows, y_testlike, held, cfg)
                    else:
                        raw[held] = fit_lr(x, pool_rows, y_testlike, groups, held, cfg)
                name = f"{cfg.get('label', 'cypmol')}_{fname}_{kind}"
                rank = percentile_by_fold(raw, fold)
                scores[name] = rank
                scores[f"L4+{name}"] = (legacy_rank + pool_rank + rank) / 3
                if cfg.get("swap"):
                    scores[f"swap_{name}"] = (legacy_rank + rank) / 2
                labelled = np.isfinite(y_official)
                standalone.append(
                    dict(
                        split=split,
                        model=name,
                        ap=float(average_precision_score(y_official[labelled], rank[labelled])),
                        auc=float(roc_auc_score(y_official[labelled], rank[labelled])),
                    )
                )
                print(split, name, standalone[-1], flush=True)
        for label, rank in (("legacy", legacy_rank), ("pool", pool_rank), ("L4", scores["L4"])):
            labelled = np.isfinite(y_official)
            standalone.append(
                dict(
                    split=split,
                    model=label,
                    ap=float(average_precision_score(y_official[labelled], rank[labelled])),
                    auc=float(roc_auc_score(y_official[labelled], rank[labelled])),
                )
            )
        official, pseudo = role == "official", role == "pseudo_negative"
        in_pool = np.array([i in analog for i in wide["compound_id"].to_list()])
        populations = [
            (
                "O",
                cfg["prior"],
                "dev",
                om.official_weights(y_official, official, cfg["prior"]),
                y_official,
                None,
            ),
            (
                "M",
                cfg["phi"],
                "dev",
                om.mixture_weights(official, pseudo, cfg["phi"]),
                y_testlike,
                None,
            ),
            (
                "O",
                cfg["prior"],
                "analog_pool",
                om.official_weights(y_official, official & in_pool, cfg["prior"]),
                y_official,
                in_pool,
            ),
        ]
        for metric, _parameter, population, weights, labels, pool in populations:

            def calls_for(score, pool=pool, fold=fold, rate=rate):
                if pool is None:
                    return om.top_calls(score, fold, rate)
                calls = np.zeros(len(score), dtype=bool)
                members = np.flatnonzero(pool)
                top = members[
                    np.argsort(-score[members], kind="stable")[: round(rate * len(members))]
                ]
                calls[top] = True
                return calls

            reference = calls_for(scores["L4"])
            for model, score in scores.items():
                if model == "L4":
                    continue
                calls = calls_for(score)
                results.append(
                    dict(
                        split=split,
                        population=population,
                        metric=metric,
                        model=model,
                        mcc_l4=om.weighted_mcc(labels, reference, weights),
                        mcc_model=om.weighted_mcc(labels, calls, weights),
                        **om.group_bootstrap(
                            labels,
                            groups,
                            weights,
                            calls,
                            reference,
                            cfg["bootstrap_samples"],
                            cfg["seed"],
                        ),
                    )
                )
                if model.startswith("swap_"):
                    for count in cfg.get("counts", []):
                        wider = calls_for(score, rate=count / 750)
                        count_rows.append(
                            dict(
                                split=split,
                                population=population,
                                metric=metric,
                                model=model,
                                calls=count,
                                mcc_49=om.weighted_mcc(labels, calls, weights),
                                mcc_calls=om.weighted_mcc(labels, wider, weights),
                                **om.group_bootstrap(
                                    labels,
                                    groups,
                                    weights,
                                    wider,
                                    calls,
                                    cfg["bootstrap_samples"],
                                    cfg["seed"],
                                ),
                            )
                        )
    pl.DataFrame(standalone).write_csv(out / "standalone.csv")
    if count_rows:
        pl.DataFrame(count_rows).write_csv(out / "calls-versus-49.csv")
    pl.DataFrame(results).write_csv(out / "versus-l4.csv")
    print(pl.DataFrame(standalone))


if __name__ == "__main__":
    main()
