"""A4 TDI pool and activity-gate variants against the submitted equal-seven ensemble.

Unlike CYP2D6, the supplied CYP3A4_is_TDI labels already cover primary-screen
non-hits: 525 labelled molecules have log2fc >= -0.5 and only one of them is
positive. The labelled set is therefore not confined to primary hits and the
official labels alone define the evaluation target here. The remaining
pseudo-negative population is small (the unlabelled screened non-hits), so model B
is expected to differ little from model A; that expectation is reported, not assumed.

All models, the activity gate, its Platt calibration, feature scaling and the
logistic regularization choice are fitted inside each outer training split; held-fold
molecules never enter any fit. Positive calls are taken from every development
molecule of a held fold and scored on the supplied labels only. The submitted
equal-seven ensemble enters as a frozen complete-pool UMAP comparator and as the base
of the gate combinations; it has no scaffold counterpart.
"""

import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl
from cyp_primary_supervision import load_primary
from joblib import Parallel, delayed
from lightgbm import LGBMClassifier, LGBMRegressor
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

PRIMARY_ORDER = ["CYP2D6", "CYP3A4", "CYP1A2", "CYP2C9"]
CONTROL = "official_lgbm"
FITTED = ["official_lgbm", "pool_lgbm", "gate_x_tdi", "chemeleon_lr"]
POOL_ROLES = ["official", "pseudo_negative", "borderline", "unlabelled_hit", "unscreened"]


def featurize(smiles):
    """Fresh Morgan r2/2048 bits and RDKit 2D descriptors."""
    RDLogger.DisableLog("rdApp.*")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = np.zeros((len(smiles), 2048), dtype=np.float32)
    desc = np.full((len(smiles), len(Descriptors._descList)), np.nan, dtype=np.float32)
    for i, text in enumerate(smiles):
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            raise ValueError(f"Unparseable canonical molecule: {text}")
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(mol), fp[i])
        for j, (_, function) in enumerate(Descriptors._descList):
            try:
                desc[i, j] = function(mol)
            except (ValueError, RuntimeError, OverflowError, ZeroDivisionError):
                pass
    desc[~np.isfinite(desc)] = np.nan
    return fp, desc


def lgbm_params(cfg, seed):
    return dict(
        **cfg["lgbm"],
        verbosity=-1,
        n_jobs=cfg["lgbm_threads"],
        random_state=seed,
        deterministic=True,
        force_col_wise=True,
    )


def fit_classifier(x, rows, y, predict, seeds, cfg):
    scores = []
    for seed in seeds:
        model = LGBMClassifier(class_weight="balanced", **lgbm_params(cfg, seed))
        model.fit(x[rows], y[rows].astype(int))
        scores.append([model.predict_proba(x[p])[:, 1] for p in predict])
    return [np.mean([s[i] for s in scores], axis=0) for i in range(len(predict))]


def fit_regressor(x, rows, target, predict, seeds, cfg):
    predictions = []
    for seed in seeds:
        model = LGBMRegressor(**lgbm_params(cfg, seed))
        model.fit(x[rows], target[rows])
        predictions.append([model.predict(x[p]) for p in predict])
    return [np.mean([s[i] for s in predictions], axis=0) for i in range(len(predict))]


def fit_gate(x, screened, log2fc, groups, predict, seeds, cfg):
    """LightGBM log2fc regressor Platt-calibrated to P(log2fc < gate_threshold)."""
    inner = np.full(len(screened), np.nan)
    for fit, valid in GroupKFold(n_splits=cfg["inner_folds"]).split(
        screened, groups=groups[screened]
    ):
        inner[valid] = fit_regressor(x, screened[fit], log2fc, [screened[valid]], seeds, cfg)[0]
    active = (log2fc[screened] < cfg["gate_threshold"]).astype(int)
    calibration = LogisticRegression(C=1e4, max_iter=1000).fit(inner[:, None], active)
    predicted = fit_regressor(x, screened, log2fc, predict, seeds, cfg)
    info = dict(
        screened_train=len(screened),
        active_train=int(active.sum()),
        calibration_coef=float(calibration.coef_[0, 0]),
        calibration_intercept=float(calibration.intercept_[0]),
        inner_spearman=float(
            pl.DataFrame({"a": inner, "b": log2fc[screened]})
            .select(pl.corr("a", "b", method="spearman"))
            .item()
        ),
    )
    return [calibration.predict_proba(p[:, None])[:, 1] for p in predicted], info


def fit_chemeleon_lr(x, rows, y, groups, predict, cfg):
    """Balanced logistic regression on standardized CheMeleon; C chosen by inner AP."""
    splitter = GroupKFold(n_splits=cfg["inner_folds"])
    tuning = {}
    for c in cfg["lr_C"]:
        inner = np.full(len(rows), np.nan)
        for fit, valid in splitter.split(rows, groups=groups[rows]):
            scaler = StandardScaler().fit(x[rows[fit]])
            model = LogisticRegression(
                C=c,
                class_weight="balanced",
                max_iter=cfg["lr_max_iter"],
                random_state=cfg["seed"],
            ).fit(scaler.transform(x[rows[fit]]), y[rows[fit]].astype(int))
            inner[valid] = model.decision_function(scaler.transform(x[rows[valid]]))
        tuning[str(c)] = float(average_precision_score(y[rows].astype(int), inner))
    best = float(max(cfg["lr_C"], key=lambda c: tuning[str(c)]))
    scaler = StandardScaler().fit(x[rows])
    model = LogisticRegression(
        C=best, class_weight="balanced", max_iter=cfg["lr_max_iter"], random_state=cfg["seed"]
    ).fit(scaler.transform(x[rows]), y[rows].astype(int))
    scores = [model.predict_proba(scaler.transform(x[p]))[:, 1] for p in predict]
    return scores, dict(selected_C=best, inner_ap=tuning)


def fit_predict(task, train, predict, groups, data, cfg):
    started = time.time()
    x1, x2 = data["x_morgan"], data["x_chemeleon"]
    y_official, y_pool, log2fc = data["y_official"], data["y_pool"], data["log2fc"]
    seeds = [cfg["seed"] + s for s in range(cfg["lgbm_seeds"])]
    official = train[np.isfinite(y_official[train])]
    pool = train[np.isfinite(y_pool[train])]
    screened = train[np.isfinite(log2fc[train])]
    info = dict(
        task=task,
        train=len(train),
        official_train=len(official),
        official_positive=int(np.nansum(y_official[official])),
        pool_train=len(pool),
        pseudo_negatives=len(pool) - len(official),
        lgbm_seeds=seeds,
    )
    result = {}
    p_tdi = fit_classifier(x1, official, y_official, predict, seeds, cfg)
    result["official_lgbm"] = [(s, None) for s in p_tdi]
    pooled = fit_classifier(x1, pool, y_pool, predict, seeds, cfg)
    result["pool_lgbm"] = [(s, None) for s in pooled]
    p_active, info["gate"] = fit_gate(x1, screened, log2fc, groups, predict, seeds, cfg)
    result["gate_x_tdi"] = [(a * t, a) for a, t in zip(p_active, p_tdi, strict=True)]
    chem, info["chemeleon_lr"] = fit_chemeleon_lr(x2, official, y_official, groups, predict, cfg)
    result["chemeleon_lr"] = [(s, None) for s in chem]
    info["seconds"] = round(time.time() - started, 1)
    print(json.dumps(info), flush=True)
    return task, result, info


def percentile_rank(values):
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values))
    ranks[order] = np.arange(len(values))
    return (ranks + 0.5) / len(values)


def top_calls(score, fold, rate):
    """Call the top round(rate * n) molecules of each held fold; fold < 0 is excluded."""
    calls = np.zeros(len(score), dtype=bool)
    for f in np.unique(fold[fold >= 0]):
        rows = np.flatnonzero(fold == f)
        k = int(round(rate * len(rows)))
        calls[rows[np.argsort(-score[rows], kind="stable")[:k]]] = True
    return calls


def mcc_from_counts(tp, fp, fn, tn):
    tp, fp, fn, tn = (np.asarray(v, dtype=float) for v in (tp, fp, fn, tn))
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return np.divide(tp * tn - fp * fn, denom, out=np.zeros_like(denom), where=denom > 0)


def confusion(y, calls):
    keep = np.isfinite(y)
    y, calls = y[keep].astype(bool), calls[keep]
    return (
        int((y & calls).sum()),
        int((~y & calls).sum()),
        int((y & ~calls).sum()),
        int((~y & ~calls).sum()),
    )


def mcc(y, calls):
    return float(mcc_from_counts(*confusion(y, calls)))


def ranking(y, score):
    keep = np.isfinite(y)
    y, score = y[keep].astype(int), score[keep]
    if len(np.unique(y)) < 2:
        return None, None
    return float(roc_auc_score(y, score)), float(average_precision_score(y, score))


def decision_row(frame, calls, **keys):
    y_official = frame["y_official"].to_numpy()
    role = frame["pool_role"].to_numpy()
    tp, fp, fn, tn = confusion(y_official, calls)
    row = dict(
        **keys,
        n_molecules=len(calls),
        n_calls=int(calls.sum()),
        mcc_official=float(mcc_from_counts(tp, fp, fn, tn)),
        tp=tp,
        fp=fp,
        fn=fn,
        precision=tp / (tp + fp) if tp + fp else None,
        recall=tp / (tp + fn) if tp + fn else None,
    )
    for name in POOL_ROLES:
        row[f"calls_{name}"] = int((calls & (role == name)).sum())
    return row


def group_bootstrap(frame, calls_model, calls_control, samples, seed):
    """Paired group bootstrap of pooled MCC(model) - MCC(control) on official labels."""
    y = frame["y_official"].to_numpy()
    keep = np.isfinite(y)
    y = y[keep].astype(bool)
    groups = frame["group"].to_numpy()[keep]
    _, index = np.unique(groups, return_inverse=True)
    count = index.max() + 1
    counts = []
    for calls in (calls_model[keep], calls_control[keep]):
        counts.append(
            [
                np.bincount(index, weights=(y & calls), minlength=count),
                np.bincount(index, weights=(~y & calls), minlength=count),
                np.bincount(index, weights=(y & ~calls), minlength=count),
                np.bincount(index, weights=(~y & ~calls), minlength=count),
            ]
        )
    draws = np.random.default_rng(seed).integers(0, count, (samples, count))
    weights = np.zeros((samples, count))
    for b in range(samples):
        weights[b] = np.bincount(draws[b], minlength=count)
    delta = mcc_from_counts(*(weights @ c for c in counts[0])) - mcc_from_counts(
        *(weights @ c for c in counts[1])
    )
    point = mcc_from_counts(*(c.sum() for c in counts[0])) - mcc_from_counts(
        *(c.sum() for c in counts[1])
    )
    return dict(
        delta=float(point),
        boot_mean=float(delta.mean()),
        ci_low=float(np.quantile(delta, 0.025)),
        ci_high=float(np.quantile(delta, 0.975)),
        p_delta_gt_0=float((delta > 0).mean()),
        groups=int(count),
        molecules=int(keep.sum()),
    )


def permutation_difference(values_removed, values_added, samples, seed):
    """Two-sided permutation p for mean(added) - mean(removed) over the pooled union."""
    pooled = np.concatenate([values_removed, values_added])
    observed = values_added.mean() - values_removed.mean()
    rng = np.random.default_rng(seed)
    n = len(values_removed)
    draws = np.empty(samples)
    for b in range(samples):
        shuffled = rng.permutation(pooled)
        draws[b] = shuffled[n:].mean() - shuffled[:n].mean()
    return float(observed), float((np.abs(draws) >= abs(observed) - 1e-12).mean())


def fmt(value, digits=3):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "–"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return f"{value:.{digits}f}" if isinstance(value, float) else str(value)


def markdown_table(rows, columns):
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(c)) for c in columns) + " |")
    return "\n".join(lines)


def main():
    RDLogger.DisableLog("rdApp.*")
    started = time.time()
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    inputs = cfg["inputs"]
    endpoint, comparator = cfg["endpoint"], cfg["comparator"]
    compounds = pl.read_parquet(source / inputs["compounds"]).sort("compound_id")
    aliases = pl.read_parquet(source / inputs["aliases"])
    test_mapping = pl.read_parquet(source / inputs["test_mapping"])
    measured = pl.read_parquet(source / inputs["measurements"])
    assert measured.select("compound_id", "endpoint").is_duplicated().sum() == 0
    assert not measured["endpoint"].str.contains("log2fc").any()
    labels = measured.filter(pl.col("endpoint") == endpoint).select(
        "compound_id", pl.col("value").alias("label")
    )
    table = compounds.join(labels, on="compound_id", how="left", maintain_order="left")
    primary, primary_rows = load_primary(source / inputs["primary_screen"], aliases, table)
    log2fc = primary[:, PRIMARY_ORDER.index(cfg["primary_enzyme"])]
    assignments = pl.read_parquet(source / inputs["assignments"])
    is_blind = table["compound_id"].is_in(test_mapping["compound_id"].unique().to_list()).to_numpy()
    blind = np.flatnonzero(is_blind)
    development = np.flatnonzero(~is_blind)
    y_raw = table["label"].to_numpy().astype(float)
    assert len(blind) == 750 and not np.isfinite(y_raw[blind]).any()
    assert not np.isfinite(log2fc[blind]).any()
    assert set(assignments["compound_id"].unique().to_list()) == set(
        table["compound_id"][development].to_list()
    )

    # Pool roles.
    role = np.full(len(table), "blind", dtype=object)
    labelled, screened = np.isfinite(y_raw), np.isfinite(log2fc)
    p, dev = cfg["pool"], ~is_blind
    role[dev & labelled] = "official"
    role[dev & ~labelled & screened & (log2fc >= p["pseudo_negative_min_log2fc"])] = (
        "pseudo_negative"
    )
    role[
        dev
        & ~labelled
        & screened
        & (log2fc >= p["borderline_min_log2fc"])
        & (log2fc < p["pseudo_negative_min_log2fc"])
    ] = "borderline"
    role[dev & ~labelled & screened & (log2fc < p["borderline_min_log2fc"])] = "unlabelled_hit"
    role[dev & ~labelled & ~screened] = "unscreened"
    y_official = np.where(role == "official", y_raw, np.nan)
    y_pool = np.where(role == "pseudo_negative", 0.0, y_official)
    official_rows = role == "official"
    pool_audit = dict(
        counts={name: int((role == name).sum()) for name in np.unique(role)},
        development=len(development),
        official_positive=int(np.nansum(y_official)),
        official_positive_rate=float(np.nanmean(y_official[official_rows])),
        official_with_log2fc=int((official_rows & screened).sum()),
        official_weak_primary=int((official_rows & screened & (log2fc >= -0.5)).sum()),
        official_weak_primary_positive=int(
            (official_rows & screened & (log2fc >= -0.5) & (y_raw == 1)).sum()
        ),
        official_positive_log2fc_max=float(np.nanmax(log2fc[official_rows & (y_raw == 1)])),
        official_positive_log2fc_ge_gate=int(
            (official_rows & (y_raw == 1) & (log2fc >= cfg["gate_threshold"])).sum()
        ),
        pool_n=int(np.isfinite(y_pool).sum()),
        pool_positive_rate=float(np.nanmean(y_pool)),
        screened_development=int((dev & screened).sum()),
        gate_active_development=int((dev & (log2fc < cfg["gate_threshold"])).sum()),
        primary_source_rows=int(primary_rows["source_rows"].sum()),
    )
    print(json.dumps(pool_audit), flush=True)

    # Features.
    smiles = table["standardized_smiles"].to_list()
    chunks = np.array_split(np.arange(len(smiles)), cfg["feature_jobs"])
    parts = Parallel(n_jobs=cfg["feature_jobs"])(
        delayed(featurize)([smiles[i] for i in chunk]) for chunk in chunks
    )
    fp = np.concatenate([a for a, _ in parts])
    desc = np.concatenate([b for _, b in parts])
    x_morgan = np.concatenate([fp, desc], axis=1)
    embeddings = pl.read_parquet(source / inputs["chemeleon"])
    columns = [c for c in embeddings.columns if c.startswith("chemeleon_")]
    assert len(columns) == 2048
    x_chemeleon = (
        table.select("compound_id")
        .join(embeddings, on="compound_id", how="left", validate="1:1", maintain_order="left")
        .select(columns)
        .to_numpy()
        .astype(np.float32)
    )
    assert np.isfinite(x_chemeleon).all()
    print(
        json.dumps(
            dict(
                features_seconds=round(time.time() - started, 1),
                morgan_rdkit=x_morgan.shape[1],
                chemeleon=x_chemeleon.shape[1],
            )
        ),
        flush=True,
    )
    data = dict(
        x_morgan=x_morgan,
        x_chemeleon=x_chemeleon,
        y_official=y_official,
        y_pool=y_pool,
        log2fc=log2fc,
    )

    # Outer tasks.
    family_groups, family_folds, tasks = {}, {}, []
    for split in cfg["split_strategies"]:
        ordered = table.select("compound_id").join(
            assignments.filter(pl.col("strategy") == split).select("compound_id", "fold", "group"),
            on="compound_id",
            how="left",
            maintain_order="left",
        )
        family_folds[split] = ordered["fold"].fill_null(-1).to_numpy()
        family_groups[split] = ordered["group"].fill_null("blind").to_numpy().astype(str)
        for f in range(cfg["outer_folds"]):
            train = development[family_folds[split][development] != f]
            held = development[family_folds[split][development] == f]
            assert set(family_groups[split][train]).isdisjoint(family_groups[split][held])
            tasks.append(((split, int(f)), train, [held], family_groups[split]))
    tasks.append((("blind", -1), development, [blind], family_groups[cfg["refit_inner_groups"]]))
    outputs = Parallel(n_jobs=cfg["task_jobs"])(
        delayed(fit_predict)(key, train, predict, groups, data, cfg)
        for key, train, predict, groups in tasks
    )
    fit_info = [info for _, _, info in outputs]

    # Keyed OOF and blind tables for the fitted models.
    held_rows = {key: predict[0] for key, _, predict, _ in tasks}
    names = test_mapping.select("compound_id", "molecule_name").unique()
    oof_parts, blind_parts = [], []
    for (split, f), result, _ in outputs:
        rows = held_rows[(split, f)]
        for model, ((score, active),) in result.items():
            part = pl.DataFrame(
                dict(
                    compound_id=table["compound_id"].gather(rows),
                    model=[model] * len(rows),
                    score=np.asarray(score, dtype=np.float64),
                    p_active=np.asarray(
                        active if active is not None else np.full(len(rows), np.nan),
                        dtype=np.float64,
                    ),
                )
            )
            if split == "blind":
                blind_parts.append(part)
                continue
            oof_parts.append(
                part.with_columns(
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit(split).alias("split_id"),
                    pl.lit(f).alias("fold"),
                    pl.Series("pool_role", role[rows].astype(str)),
                    pl.Series("y_official", y_official[rows]),
                    pl.Series("group", family_groups[split][rows]),
                )
            )
    oof = pl.concat(oof_parts)
    blind_scores = pl.concat(blind_parts)

    # Frozen equal-seven comparator on the umap folds.
    equal7 = pl.read_parquet(source / inputs["equal7_oof"])
    assert equal7["split_id"].unique().to_list() == [cfg["comparator_split"]]
    template = oof.filter(
        (pl.col("model") == CONTROL) & (pl.col("split_id") == cfg["comparator_split"])
    ).drop("score", "p_active")
    joined = template.join(
        equal7.select("compound_id", pl.col("fold").alias("equal7_fold"), "y_true", "score"),
        on="compound_id",
        how="left",
    )
    equal7_audit = dict(
        rows=joined.height,
        missing=int(joined["score"].null_count()),
        fold_mismatch=int((joined["fold"] != joined["equal7_fold"]).sum()),
        label_mismatch=int(
            joined.filter(
                (pl.col("y_true").is_null() != pl.col("y_official").is_nan())
                | (pl.col("y_true").is_not_null() & (pl.col("y_true") != pl.col("y_official")))
            ).height
        ),
    )
    assert equal7_audit["missing"] == 0 and equal7_audit["fold_mismatch"] == 0, equal7_audit
    assert equal7_audit["label_mismatch"] == 0, equal7_audit
    oof = pl.concat(
        [
            oof,
            joined.drop("equal7_fold", "y_true")
            .with_columns(
                pl.lit(comparator).alias("model"), pl.lit(None, dtype=pl.Float64).alias("p_active")
            )
            .select(oof.columns),
        ]
    )
    blind_equal7 = pl.read_parquet(source / inputs["equal7_blind"]).select(
        "compound_id", pl.col("a4_prediction").alias("score")
    )
    blind_scores = pl.concat(
        [
            blind_scores,
            pl.DataFrame({"compound_id": table["compound_id"].gather(blind)})
            .join(blind_equal7, on="compound_id", how="left")
            .with_columns(
                pl.lit(comparator).alias("model"), pl.lit(None, dtype=pl.Float64).alias("p_active")
            )
            .select(blind_scores.columns),
        ]
    )
    assert blind_scores["score"].null_count() == 0
    print(json.dumps(dict(equal7_audit=equal7_audit)), flush=True)

    # Combination variants, built per held fold (and once over the blind set).
    def add_variants(frame, keys):
        wide = (
            frame.filter(pl.col("model").is_in([comparator, CONTROL, "chemeleon_lr", "gate_x_tdi"]))
            .select([*keys, "compound_id", "model", "score", "p_active"])
            .pivot(on="model", index=[*keys, "compound_id"], values=["score", "p_active"])
        )
        score_of = {m: f"score_{m}" for m in [comparator, CONTROL, "chemeleon_lr", "gate_x_tdi"]}
        wide = wide.with_columns(
            pl.col("p_active_gate_x_tdi").alias("p_active"),
        )
        extra = []
        for key_values, part in wide.group_by(keys, maintain_order=True):
            e = percentile_rank(part[score_of[comparator]].to_numpy())
            a = percentile_rank(part[score_of[CONTROL]].to_numpy())
            d = percentile_rank(part[score_of["chemeleon_lr"]].to_numpy())
            gate = part["p_active"].to_numpy()
            variants = {
                "equal7_x_gate": e * gate,
                "equal7_x_sqrt_gate": e * np.sqrt(gate),
                "rank_mean_equal7_official": 0.5 * (e + a),
                "rank_mean_equal7_chemeleon": 0.5 * (e + d),
            }
            for tau in cfg["gate_filters"]:
                variants[f"equal7_filter_{tau:g}"] = np.where(gate >= tau, e, e - 1.0)
            for model, values in variants.items():
                extra.append(
                    part.select([*keys, "compound_id"]).with_columns(
                        pl.lit(model).alias("model"),
                        pl.Series("score", values),
                        pl.Series("p_active", gate),
                    )
                )
            del key_values
        return pl.concat(extra)

    variant_oof = add_variants(
        oof.filter(pl.col("split_id") == cfg["comparator_split"]), ["split_id", "fold"]
    )
    base = oof.filter(
        (pl.col("model") == CONTROL) & (pl.col("split_id") == cfg["comparator_split"])
    )
    variant_oof = variant_oof.join(
        base.select("compound_id", "endpoint", "pool_role", "y_official", "group"),
        on="compound_id",
        how="left",
    )
    oof = pl.concat([oof, variant_oof.select(oof.columns)])
    variant_blind = add_variants(
        blind_scores.with_columns(pl.lit("blind").alias("split")), ["split"]
    ).drop("split")
    blind_scores = pl.concat([blind_scores, variant_blind.select(blind_scores.columns)])

    oof = oof.with_columns(pl.col(["score", "p_active", "y_official"]).fill_nan(None)).sort(
        "split_id", "model", "fold", "compound_id"
    )
    oof.drop("group").write_parquet(out / "oof_scores.parquet")
    blind_scores = (
        blind_scores.join(names, on="compound_id", how="left")
        .select("molecule_name", "compound_id", "model", "score", "p_active")
        .with_columns(pl.col(["score", "p_active"]).fill_nan(None))
        .sort("model", "compound_id")
    )
    assert blind_scores["molecule_name"].null_count() == 0
    blind_scores.write_parquet(out / "test_scores.parquet")

    # Evaluation.
    rates, reference = cfg["rates"], cfg["reference_rate"]
    assert reference in rates
    ranking_rows, curve_rows, nested_rows, boot_rows, deploy = [], [], [], [], {}
    for split in cfg["split_strategies"]:
        present = (
            oof.filter(pl.col("split_id") == split)["model"].unique(maintain_order=True).to_list()
        )
        fixed_calls, nested_calls = {}, {}
        for model in present:
            frame = oof.filter((pl.col("split_id") == split) & (pl.col("model") == model)).sort(
                "fold", "compound_id"
            )
            frame = frame.with_columns(pl.col("y_official").fill_null(np.nan))
            score = frame["score"].to_numpy()
            fold = frame["fold"].to_numpy()
            y_off = frame["y_official"].to_numpy()
            auc, ap = ranking(y_off, score)
            fold_rank = [ranking(y_off[fold == f], score[fold == f]) for f in np.unique(fold)]
            ranking_rows.append(
                dict(
                    split_id=split,
                    model=model,
                    auroc=auc,
                    ap=ap,
                    fold_mean_auroc=float(np.mean([m[0] for m in fold_rank])),
                    fold_mean_ap=float(np.mean([m[1] for m in fold_rank])),
                    n_official=int(np.isfinite(y_off).sum()),
                    n_pool=len(score),
                )
            )
            calls_by_rate = {}
            for rate in rates:
                calls = top_calls(score, fold, rate)
                calls_by_rate[rate] = calls
                curve_rows.append(
                    decision_row(frame, calls, split_id=split, model=model, rate=rate, fold=-1)
                )
                for f in np.unique(fold):
                    mask = fold == f
                    curve_rows.append(
                        decision_row(
                            frame.filter(pl.Series(mask)),
                            calls[mask],
                            split_id=split,
                            model=model,
                            rate=rate,
                            fold=int(f),
                        )
                    )
            fixed_calls[model] = (frame, calls_by_rate[reference])
            selected = np.zeros(len(score), dtype=bool)
            for f in np.unique(fold):
                other = fold != f
                inner = {r: mcc(y_off[other], calls_by_rate[r][other]) for r in rates}
                rate = max(rates, key=lambda r: (inner[r], -r))
                held = fold == f
                selected[held] = calls_by_rate[rate][held]
                nested_rows.append(
                    decision_row(
                        frame.filter(pl.Series(held)),
                        calls_by_rate[rate][held],
                        split_id=split,
                        model=model,
                        fold=int(f),
                        selected_rate=rate,
                        selection_mcc=inner[rate],
                    )
                )
            nested_rows.append(
                decision_row(
                    frame, selected, split_id=split, model=model, fold=-1, selected_rate=None
                )
            )
            nested_calls[model] = (frame, selected)
            full = {r: mcc(y_off, calls_by_rate[r]) for r in rates}
            deploy.setdefault(split, {})[model] = max(rates, key=lambda r: (full[r], -r))
        references = [comparator, CONTROL] if comparator in present else [CONTROL]
        for control in references:
            control_frame = nested_calls[control][0]
            for model in present:
                if model == control:
                    continue
                for operating, calls_map in [("nested", nested_calls), ("reference", fixed_calls)]:
                    frame, calls = calls_map[model]
                    _, control_calls = calls_map[control]
                    assert frame["compound_id"].equals(control_frame["compound_id"])
                    boot_rows.append(
                        dict(
                            split_id=split,
                            model=model,
                            control=control,
                            operating_point=operating,
                            **group_bootstrap(
                                frame,
                                calls,
                                control_calls,
                                cfg["bootstrap_samples"],
                                cfg["seed"],
                            ),
                        )
                    )
    ranking_table = pl.DataFrame(ranking_rows)
    curves = pl.DataFrame(curve_rows, infer_schema_length=None)
    nested = pl.DataFrame(nested_rows, infer_schema_length=None)
    bootstrap = pl.DataFrame(boot_rows)
    ranking_table.write_csv(out / "ranking_metrics.csv")
    curves.write_csv(out / "rate_curves.csv")
    nested.write_csv(out / "nested_selection.csv")
    bootstrap.write_csv(out / "bootstrap_deltas.csv")

    # Blind diagnostics against the submitted 249 A4 positives.
    submission = pl.read_parquet(source / inputs["submission"]).select(
        pl.col("Molecule_Name").alias("molecule_name"), pl.col(endpoint).alias("call")
    )
    assert submission.height == 750
    submitted = set(
        submission.join(names, on="molecule_name", how="left")
        .filter(pl.col("call"))["compound_id"]
        .to_list()
    )
    count = cfg["reference_count"]
    assert len(submitted) == count
    comparator_blind = blind_scores.filter(pl.col("model") == comparator).sort(
        ["score", "compound_id"], descending=[True, False]
    )
    assert set(comparator_blind["compound_id"].to_list()[:count]) == submitted
    blind_rows = []
    for model in blind_scores["model"].unique(maintain_order=True).to_list():
        frame = blind_scores.filter(pl.col("model") == model).sort(
            ["score", "compound_id"], descending=[True, False]
        )
        ids = frame["compound_id"].to_list()
        gate = frame["p_active"].to_numpy()
        blind_rows.append(
            dict(
                model=model,
                overlap_at_249=len(set(ids[:count]) & submitted),
                added_vs_submitted=count - len(set(ids[:count]) & submitted),
                mean_p_active=float(np.nanmean(gate)) if np.isfinite(gate).any() else None,
                p_active_gt_0_5=(int((gate > 0.5).sum()) if np.isfinite(gate).all() else None),
                **{
                    f"{split}_nested_rates": ",".join(
                        f"{r:.3f}"
                        for r in nested.filter(
                            (pl.col("split_id") == split)
                            & (pl.col("model") == model)
                            & (pl.col("fold") >= 0)
                        )["selected_rate"].to_list()
                    )
                    for split in cfg["split_strategies"]
                },
            )
        )
    blind_table = pl.DataFrame(blind_rows, infer_schema_length=None)
    blind_table.write_csv(out / "blind_summary.csv")
    gate_blind = blind_scores.filter(pl.col("model") == "gate_x_tdi")["p_active"].to_numpy()
    blind_gate = dict(
        n=len(gate_blind),
        mean=float(gate_blind.mean()),
        quantiles={str(q): float(np.quantile(gate_blind, q)) for q in [0.1, 0.25, 0.5, 0.75, 0.9]},
        gt_0_5=int((gate_blind > 0.5).sum()),
        **{f"lt_{tau:g}": int((gate_blind < tau).sum()) for tau in cfg["gate_filters"]},
        development_active_rate=pool_audit["gate_active_development"]
        / pool_audit["screened_development"],
    )

    # Historical transition concordance on the blind scores.
    payloads = {}
    for name, alias in inputs["payloads"]:
        frame = pl.read_parquet(source / alias).select(
            pl.col("Molecule_Name").alias("molecule_name"), pl.col(endpoint).alias("call")
        )
        assert frame.height == 750
        payloads[name] = set(
            frame.join(names, on="molecule_name", how="left")
            .filter(pl.col("call"))["compound_id"]
            .to_list()
        )
    score_columns = {}
    for model in sorted(blind_scores["model"].unique().to_list()):
        frame = blind_scores.filter(pl.col("model") == model)
        score_columns[model] = dict(
            zip(frame["compound_id"].to_list(), frame["score"].to_list(), strict=True)
        )
    gate_frame = blind_scores.filter(pl.col("model") == "gate_x_tdi")
    score_columns["gate_p_active"] = dict(
        zip(gate_frame["compound_id"].to_list(), gate_frame["p_active"].to_list(), strict=True)
    )
    concordance_rows = []
    for earlier, later, lb_delta in cfg["transitions"]:
        removed = sorted(payloads[earlier] - payloads[later])
        added = sorted(payloads[later] - payloads[earlier])
        for column, mapping in score_columns.items():
            a = np.array([mapping[i] for i in removed])
            b = np.array([mapping[i] for i in added])
            difference, p_value = permutation_difference(
                a, b, cfg["permutation_samples"], cfg["seed"]
            )
            concordance_rows.append(
                dict(
                    transition=f"{earlier}->{later}",
                    lb_delta=lb_delta,
                    column=column,
                    n_removed=len(removed),
                    n_added=len(added),
                    mean_removed=float(a.mean()),
                    mean_added=float(b.mean()),
                    difference=difference,
                    permutation_p=p_value,
                    sign_agrees=bool(np.sign(difference) == np.sign(lb_delta)),
                )
            )
    concordance = pl.DataFrame(concordance_rows)
    concordance.write_csv(out / "transition_concordance.csv")
    (out / "fit_info.json").write_text(json.dumps(fit_info, indent=2) + "\n")

    decision_main = []
    for r in nested_rows:
        if r["fold"] != -1:
            continue
        fixed = next(
            c
            for c in curve_rows
            if c["fold"] == -1
            and c["split_id"] == r["split_id"]
            and c["model"] == r["model"]
            and c["rate"] == reference
        )
        decision_main.append(
            dict(
                split_id=r["split_id"],
                model=r["model"],
                nested_mcc=r["mcc_official"],
                nested_calls=r["n_calls"],
                reference_mcc=fixed["mcc_official"],
                reference_precision=fixed["precision"],
                reference_recall=fixed["recall"],
                reference_calls_official=fixed["calls_official"],
                reference_calls_unscreened=fixed["calls_unscreened"],
            )
        )
    curve_wide = (
        curves.filter(pl.col("fold") == -1)
        .select("split_id", "model", "rate", "mcc_official")
        .with_columns(
            pl.col("rate").map_elements(lambda r: f"{100 * r:.1f}%", return_dtype=pl.String)
        )
        .pivot(on="rate", index=["split_id", "model"], values="mcc_official")
    )
    summary = dict(
        config=cfg,
        endpoint=endpoint,
        pool=pool_audit,
        equal7_comparator=equal7_audit,
        gate_choice=(
            "LightGBM regressor on CYP3A4 primary log2fc (all screened training molecules), "
            "Platt-calibrated to P(log2fc < gate_threshold) on inner group-fold OOF predictions"
        ),
        ranking=ranking_rows,
        decision_pooled=[r for r in curve_rows if r["fold"] == -1],
        nested=[r for r in nested_rows if r["fold"] == -1],
        nested_per_fold=[r for r in nested_rows if r["fold"] >= 0],
        bootstrap=boot_rows,
        deployment_rate=deploy,
        blind=blind_rows,
        blind_gate=blind_gate,
        transition_concordance=concordance_rows,
        runtime_seconds=round(time.time() - started, 1),
        caveats=[
            "The equal-seven comparator exists on umap folds only; scaffold comparisons "
            "use official_lgbm as the reference instead",
            "Only 161 unlabelled screened non-hits exist for A4, so model B barely differs "
            "from model A by construction",
            "Blind rank-based calls say nothing about the blind scored-label mask",
            "Two historical A4 transitions are few and share one hidden test set",
        ],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    report = [
        "# A4 pool and activity-gate variants",
        "",
        "## Pool counts",
        "",
        "```json",
        json.dumps(pool_audit, indent=2),
        "```",
        "",
        "## Ranking (pooled held-fold, official labels)",
        "",
        markdown_table(ranking_rows, ["split_id", "model", "auroc", "ap", "n_official", "n_pool"]),
        "",
        "## Decision",
        "",
        markdown_table(decision_main, list(decision_main[0].keys())),
        "",
        "## Official-label MCC by call rate (pooled)",
        "",
        markdown_table(curve_wide.to_dicts(), curve_wide.columns),
        "",
        "## Paired group bootstrap",
        "",
        markdown_table(
            boot_rows,
            [
                "split_id",
                "model",
                "control",
                "operating_point",
                "delta",
                "ci_low",
                "ci_high",
                "p_delta_gt_0",
            ],
        ),
        "",
        "## Blind diagnostics",
        "",
        markdown_table(blind_rows, list(blind_rows[0].keys())),
        "",
        "```json",
        json.dumps(blind_gate, indent=2),
        "```",
        "",
        "## Historical A4 transition concordance",
        "",
        markdown_table(
            concordance_rows,
            [
                "transition",
                "lb_delta",
                "column",
                "n_removed",
                "n_added",
                "mean_removed",
                "mean_added",
                "difference",
                "permutation_p",
                "sign_agrees",
            ],
        ),
        "",
    ]
    (out / "report.md").write_text("\n".join(report))
    print(json.dumps(dict(runtime_seconds=summary["runtime_seconds"])), flush=True)


if __name__ == "__main__":
    main()
