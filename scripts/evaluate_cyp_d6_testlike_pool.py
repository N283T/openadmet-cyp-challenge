"""D6 TDI evaluation on a test-like pool with primary-screen pseudo-negatives.

Every supplied CYP2D6_is_TDI label belongs to a single-concentration primary-screen
hit (TDI arm, 50 uM, log2fc < -0.5). Screened molecules with log2fc >= -0.25 and no
label are treated as near-certain "assigned negatives" (pseudo-negatives). The blind
test was measured in dose response for every molecule, so official + pseudo-negative
labels form a pool whose positive rate is close to the estimated blind rate.

All models, the activity gate, its calibration, feature scalers and the logistic
regularization choice are fitted inside each outer training split. Held-fold
molecules never enter any fit. Positive calls are made from all development
molecules of a held fold (complete pool); metrics use only the scored subsets.
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

ENDPOINT = "CYP2D6_is_TDI"
CONTROL = "official_lgbm"
LEGACY = "legacy_sep15_refit"
HISTORICAL_RATE = 49 / 750


def featurize(smiles):
    """Fresh Morgan r2/2048 bits and RDKit 2D descriptors (same as zero-start)."""
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
    """Balanced LightGBM classifier averaged over seeds; returns probabilities."""
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
    """LightGBM log2fc regressor + Platt calibration to P(log2fc < gate threshold).

    The calibration is fitted on inner group-fold OOF regressor predictions of the
    screened training molecules only.
    """
    inner = np.full(len(screened), np.nan)
    splitter = GroupKFold(n_splits=cfg["inner_folds"])
    for fit, valid in splitter.split(screened, groups=groups[screened]):
        inner[valid] = fit_regressor(x, screened[fit], log2fc, [screened[valid]], seeds, cfg)[0]
    active = (log2fc[screened] < cfg["gate_threshold"]).astype(int)
    calibration = LogisticRegression(C=1e4, max_iter=1000).fit(inner[:, None], active)
    predicted = fit_regressor(x, screened, log2fc, predict, seeds, cfg)
    p_active = [calibration.predict_proba(p[:, None])[:, 1] for p in predicted]
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
    return p_active, info


def fit_chemeleon_lr(x, pool, y, groups, predict, cfg, seed):
    """Balanced logistic regression on standardized CheMeleon; C chosen by inner AP."""
    splitter = GroupKFold(n_splits=cfg["inner_folds"])
    tuning = {}
    for c in cfg["lr_C"]:
        inner = np.full(len(pool), np.nan)
        for fit, valid in splitter.split(pool, groups=groups[pool]):
            scaler = StandardScaler().fit(x[pool[fit]])
            model = LogisticRegression(
                C=c, class_weight="balanced", max_iter=cfg["lr_max_iter"], random_state=seed
            ).fit(scaler.transform(x[pool[fit]]), y[pool[fit]].astype(int))
            inner[valid] = model.decision_function(scaler.transform(x[pool[valid]]))
        tuning[str(c)] = float(average_precision_score(y[pool].astype(int), inner))
    best = float(max(cfg["lr_C"], key=lambda c: tuning[str(c)]))
    scaler = StandardScaler().fit(x[pool])
    model = LogisticRegression(
        C=best, class_weight="balanced", max_iter=cfg["lr_max_iter"], random_state=seed
    ).fit(scaler.transform(x[pool]), y[pool].astype(int))
    scores = [model.predict_proba(scaler.transform(x[p]))[:, 1] for p in predict]
    return scores, dict(selected_C=best, inner_ap=tuning)


def fit_predict(task, train, predict, groups, data, cfg):
    """Fit every recipe on `train` rows and score each array in `predict`."""
    started = time.time()
    x1, x2 = data["x_morgan"], data["x_chemeleon"]
    y_official, y_testlike, log2fc = data["y_official"], data["y_testlike"], data["log2fc"]
    seeds = [cfg["seed"] + s for s in range(cfg["lgbm_seeds"])]
    official = train[np.isfinite(y_official[train])]
    pool = train[np.isfinite(y_testlike[train])]
    screened = train[np.isfinite(log2fc[train])]
    result, info = (
        {},
        dict(
            task=task,
            train=len(train),
            official_train=len(official),
            official_positive=int(np.nansum(y_official[official])),
            pool_train=len(pool),
            lgbm_seeds=seeds,
        ),
    )
    p_tdi = fit_classifier(x1, official, y_official, predict, seeds, cfg)
    result["official_lgbm"] = [(s, None, None) for s in p_tdi]
    pooled = fit_classifier(x1, pool, y_testlike, predict, seeds, cfg)
    result["pool_lgbm"] = [(s, None, None) for s in pooled]
    p_active, info["gate_morgan"] = fit_gate(x1, screened, log2fc, groups, predict, seeds, cfg)
    result["gate_x_tdi"] = [(a * t, a, t) for a, t in zip(p_active, p_tdi, strict=True)]
    chem_pool, info["pool_chemeleon"] = fit_chemeleon_lr(
        x2, pool, y_testlike, groups, predict, cfg, cfg["seed"]
    )
    result["pool_chemeleon"] = [(s, None, None) for s in chem_pool]
    chem_tdi = fit_classifier(x2, official, y_official, predict, seeds, cfg)
    chem_active, info["gate_chemeleon"] = fit_gate(
        x2, screened, log2fc, groups, predict, seeds, cfg
    )
    result["gate_x_tdi_chemeleon"] = [
        (a * t, a, t) for a, t in zip(chem_active, chem_tdi, strict=True)
    ]
    info["seconds"] = round(time.time() - started, 1)
    print(json.dumps(info), flush=True)
    return task, result, info


def percentile_rank(values):
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values))
    ranks[order] = np.arange(len(values))
    return (ranks + 0.5) / len(values)


def top_calls(score, fold, rate):
    """Call the top round(rate * n) molecules of each fold (stable compound order)."""
    calls = np.zeros(len(score), dtype=bool)
    for f in np.unique(fold):
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
    y_testlike = frame["y_testlike"].to_numpy()
    role = frame["pool_role"].to_numpy()
    tp, fp, fn, tn = confusion(y_testlike, calls)
    row = dict(
        **keys,
        n_molecules=len(calls),
        n_calls=int(calls.sum()),
        mcc_official=mcc(y_official, calls),
        mcc_testlike=float(mcc_from_counts(tp, fp, fn, tn)),
        tp_testlike=tp,
        fp_testlike=fp,
        fn_testlike=fn,
        precision_testlike=tp / (tp + fp) if tp + fp else None,
        recall_testlike=tp / (tp + fn) if tp + fn else None,
    )
    for name in ["official", "pseudo_negative", "borderline", "unscreened", "unlabelled_hit"]:
        row[f"calls_{name}"] = int((calls & (role == name)).sum())
    return row


def group_bootstrap(frame, calls_model, calls_control, label, samples, seed):
    """Paired group bootstrap of pooled MCC(model) - MCC(control) on a label set."""
    y = frame[label].to_numpy()
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
    model = mcc_from_counts(*(weights @ c for c in counts[0]))
    control = mcc_from_counts(*(weights @ c for c in counts[1]))
    delta = model - control
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


def fmt(value, digits=3):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "–"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


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
    compounds = pl.read_parquet(source / inputs["compounds"]).sort("compound_id")
    aliases = pl.read_parquet(source / inputs["aliases"])
    test_mapping = pl.read_parquet(source / inputs["test_mapping"])
    measured = pl.read_parquet(source / inputs["measurements"])
    assert measured.select("compound_id", "endpoint").is_duplicated().sum() == 0
    assert not measured["endpoint"].str.contains("log2fc").any()
    labels = measured.filter(pl.col("endpoint") == ENDPOINT).select(
        "compound_id", pl.col("value").alias("label")
    )
    table = compounds.join(labels, on="compound_id", how="left", maintain_order="left")
    primary, primary_rows = load_primary(source / inputs["primary_screen"], aliases, table)
    log2fc = primary[:, 0]
    assignments = pl.read_parquet(source / inputs["assignments"])
    blind_ids = set(test_mapping["compound_id"].to_list())
    is_blind = table["compound_id"].is_in(list(blind_ids)).to_numpy()
    blind = np.flatnonzero(is_blind)
    development = np.flatnonzero(~is_blind)
    y_raw = table["label"].to_numpy().astype(float)
    assert len(blind) == 750 and not np.isfinite(y_raw[blind]).any()
    assert not np.isfinite(log2fc[blind]).any()
    assert set(assignments["compound_id"].unique().to_list()) == set(
        table["compound_id"][development].to_list()
    )

    # Pool roles for development molecules.
    role = np.full(len(table), "blind", dtype=object)
    labelled = np.isfinite(y_raw)
    screened = np.isfinite(log2fc)
    p = cfg["pool"]
    dev = ~is_blind
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
    y_testlike = np.where(role == "pseudo_negative", 0.0, y_official)
    counts = {name: int((role == name).sum()) for name in np.unique(role)}
    official_rows = role == "official"
    pool_audit = dict(
        counts=counts,
        development=len(development),
        official_positive=int(np.nansum(y_official)),
        official_with_log2fc=int((official_rows & screened).sum()),
        official_log2fc_max=float(np.nanmax(log2fc[official_rows])),
        official_positive_log2fc_max=float(np.nanmax(log2fc[official_rows & (y_raw == 1)])),
        official_positive_log2fc_ge_gate=int(
            (official_rows & (y_raw == 1) & (log2fc >= cfg["gate_threshold"])).sum()
        ),
        testlike_n=int(np.isfinite(y_testlike).sum()),
        testlike_positive_rate=float(np.nanmean(y_testlike)),
        screened_development=int((dev & screened).sum()),
        gate_active_development=int((dev & (log2fc < cfg["gate_threshold"])).sum()),
        primary_source_rows=int(primary_rows["source_rows"].sum()),
    )
    print(json.dumps(pool_audit), flush=True)

    # Optional smoke subsample (development only; blind is always complete).
    smoke = cfg.get("smoke")
    if smoke:
        rng = np.random.default_rng(cfg["seed"])
        development = np.sort(rng.choice(development, smoke["development_sample"], replace=False))

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
                descriptor_nan_fraction=float(np.isnan(desc).mean()),
                fingerprint_on_bits_mean=float(fp.sum(axis=1).mean()),
                chemeleon=x_chemeleon.shape[1],
            )
        ),
        flush=True,
    )
    data = dict(
        x_morgan=x_morgan,
        x_chemeleon=x_chemeleon,
        y_official=y_official,
        y_testlike=y_testlike,
        log2fc=log2fc,
    )

    # Outer tasks: fixed umap/scaffold folds plus a full-development blind refit.
    family_groups, family_folds, tasks = {}, {}, []
    for split in cfg["split_strategies"]:
        ordered = table.select("compound_id").join(
            assignments.filter(pl.col("strategy") == split).select("compound_id", "fold", "group"),
            on="compound_id",
            how="left",
            maintain_order="left",
        )
        fold = ordered["fold"].fill_null(-1).to_numpy()
        groups = ordered["group"].fill_null("blind").to_numpy().astype(str)
        family_groups[split], family_folds[split] = groups, fold
        for f in cfg.get("evaluated_folds", range(cfg["outer_folds"])):
            train = development[fold[development] != f]
            held = development[fold[development] == f]
            assert set(groups[train]).isdisjoint(groups[held])
            tasks.append(((split, int(f)), train, [held], groups))
    refit_groups = family_groups[cfg["refit_inner_groups"]]
    tasks.append((("blind", -1), development, [blind], refit_groups))
    outputs = Parallel(n_jobs=cfg["task_jobs"])(
        delayed(fit_predict)(key, train, predict, groups, data, cfg)
        for key, train, predict, groups in tasks
    )
    fit_info = [info for _, _, info in outputs]

    # Assemble keyed OOF and blind tables.
    oof_parts, blind_part = [], None
    held_rows = {key: predict[0] for key, _, predict, _ in tasks}
    names = test_mapping.select("compound_id", "molecule_name").unique()
    for (split, f), result, _ in outputs:
        rows = held_rows[(split, f)]
        for model, ((score, active, tdi),) in result.items():
            part = pl.DataFrame(
                dict(
                    compound_id=table["compound_id"].gather(rows),
                    model=[model] * len(rows),
                    score=np.asarray(score, dtype=np.float64),
                    p_active=np.asarray(
                        active if active is not None else np.full(len(rows), np.nan),
                        dtype=np.float64,
                    ),
                    p_tdi=np.asarray(
                        tdi if tdi is not None else np.full(len(rows), np.nan),
                        dtype=np.float64,
                    ),
                )
            )
            if split == "blind":
                blind_part = part if blind_part is None else pl.concat([blind_part, part])
                continue
            oof_parts.append(
                part.with_columns(
                    pl.lit(ENDPOINT).alias("endpoint"),
                    pl.lit(split).alias("split_id"),
                    pl.lit(f).alias("fold"),
                    pl.Series("pool_role", role[rows].astype(str)),
                    pl.Series("y_official", y_official[rows]),
                    pl.Series("y_testlike", y_testlike[rows]),
                    pl.Series("group", family_groups[split][rows]),
                )
            )
    oof = pl.concat(oof_parts)

    # Legacy comparator: Sep 15 refit nested OOF, seed-averaged score.
    legacy_audit = dict(used=False)
    legacy = pl.read_parquet(source / inputs["legacy_oof"])
    legacy_mean = legacy.group_by("compound_id", "split_id").agg(
        pl.col("score").mean(),
        pl.col("fold").n_unique().alias("fold_variants"),
        pl.col("fold").first(),
        pl.col("seed").n_unique().alias("seeds"),
        pl.col("label").first().alias("legacy_label"),
    )
    reference = oof.filter(pl.col("model") == CONTROL).select(
        "compound_id", "split_id", "fold", "y_official"
    )
    joined = reference.join(
        legacy_mean, on=["compound_id", "split_id"], how="left", suffix="_legacy"
    )
    label_mismatch = joined.filter(
        (pl.col("legacy_label").is_not_nan() != pl.col("y_official").is_not_nan())
        | (pl.col("legacy_label").is_not_nan() & (pl.col("legacy_label") != pl.col("y_official")))
    ).height
    legacy_audit.update(
        missing=int(joined["score"].null_count()),
        fold_mismatch=int((joined["fold"] != joined["fold_legacy"]).sum()),
        seed_counts=sorted(joined["seeds"].drop_nulls().unique().to_list()),
        fold_inconsistent_across_seeds=int((joined["fold_variants"] > 1).sum()),
        label_mismatch=int(label_mismatch),
        rows=joined.height,
    )
    if (
        legacy_audit["missing"] == 0
        and legacy_audit["fold_mismatch"] == 0
        and legacy_audit["fold_inconsistent_across_seeds"] == 0
        and legacy_audit["label_mismatch"] == 0
    ):
        legacy_audit["used"] = True
        template = oof.filter(pl.col("model") == CONTROL).drop("score", "p_active", "p_tdi")
        legacy_rows = template.join(
            legacy_mean.select("compound_id", "split_id", "score"),
            on=["compound_id", "split_id"],
            how="left",
        ).with_columns(
            pl.lit(LEGACY).alias("model"),
            pl.lit(None, dtype=pl.Float64).alias("p_active"),
            pl.lit(None, dtype=pl.Float64).alias("p_tdi"),
        )
        oof = pl.concat([oof, legacy_rows.select(oof.columns)])
    print(json.dumps(dict(legacy_audit=legacy_audit)), flush=True)

    # Equal rank-mean of the two gate models (within each held fold / within blind).
    def add_rank_mean(frame, keys):
        wide = (
            frame.filter(pl.col("model").is_in(["gate_x_tdi", "gate_x_tdi_chemeleon"]))
            .with_columns(
                (
                    (pl.col("score").rank("ordinal").over([*keys, "model"]) - 0.5)
                    / pl.len().over([*keys, "model"])
                ).alias("pct")
            )
            .group_by([*keys, "compound_id"])
            .agg(pl.col("pct").mean().alias("score"))
        )
        base = frame.filter(pl.col("model") == "gate_x_tdi").drop("score", "model")
        extra = base.join(wide, on=[*keys, "compound_id"]).with_columns(
            pl.lit("rank_mean_gates").alias("model"),
            pl.lit(None, dtype=pl.Float64).alias("p_active"),
            pl.lit(None, dtype=pl.Float64).alias("p_tdi"),
        )
        return pl.concat([frame, extra.select(frame.columns)])

    oof = add_rank_mean(oof, ["split_id", "fold"])
    oof = oof.with_columns(pl.col(["score", "p_active", "p_tdi"]).fill_nan(None))
    oof = oof.with_columns(pl.col(["y_official", "y_testlike"]).fill_nan(None))
    oof = oof.sort("split_id", "model", "fold", "compound_id")
    oof.drop("group").write_parquet(out / "oof_scores.parquet")
    blind_scores = add_rank_mean(blind_part.with_columns(pl.lit("blind").alias("split")), ["split"])
    blind_scores = (
        blind_scores.drop("split")
        .join(names, on="compound_id", how="left")
        .select("molecule_name", "compound_id", "model", "score", "p_active", "p_tdi")
        .with_columns(pl.col(["score", "p_active", "p_tdi"]).fill_nan(None))
        .sort("model", "compound_id")
    )
    assert blind_scores["molecule_name"].null_count() == 0
    blind_scores.write_parquet(out / "test_scores.parquet")

    # Evaluation.
    models = oof["model"].unique(maintain_order=True).to_list()
    rates = cfg["rates"]
    ranking_rows, curve_rows, nested_rows, boot_rows, deploy = [], [], [], [], {}
    for split in cfg["split_strategies"]:
        nested_calls = {}
        fixed_calls = {}
        for model in models:
            frame = oof.filter((pl.col("split_id") == split) & (pl.col("model") == model)).sort(
                "fold", "compound_id"
            )
            frame = frame.with_columns(pl.col(["y_official", "y_testlike"]).fill_null(np.nan))
            score = frame["score"].to_numpy()
            fold = frame["fold"].to_numpy()
            y_off = frame["y_official"].to_numpy()
            y_tl = frame["y_testlike"].to_numpy()
            auc_o, ap_o = ranking(y_off, score)
            auc_t, ap_t = ranking(y_tl, score)
            fold_metrics = [
                (
                    ranking(y_off[fold == f], score[fold == f]),
                    ranking(y_tl[fold == f], score[fold == f]),
                )
                for f in np.unique(fold)
            ]

            def mean_of(values):
                values = [v for v in values if v is not None]
                return float(np.mean(values)) if values else None

            ranking_rows.append(
                dict(
                    split_id=split,
                    model=model,
                    auroc_official=auc_o,
                    ap_official=ap_o,
                    auroc_testlike=auc_t,
                    ap_testlike=ap_t,
                    fold_mean_auroc_official=mean_of([m[0][0] for m in fold_metrics]),
                    fold_mean_ap_official=mean_of([m[0][1] for m in fold_metrics]),
                    fold_mean_auroc_testlike=mean_of([m[1][0] for m in fold_metrics]),
                    fold_mean_ap_testlike=mean_of([m[1][1] for m in fold_metrics]),
                    n_official=int(np.isfinite(y_off).sum()),
                    n_testlike=int(np.isfinite(y_tl).sum()),
                    n_scored_pool=len(score),
                )
            )
            calls_by_rate = {}
            for rate in [*rates, HISTORICAL_RATE]:
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
            fixed_calls[model] = (frame, calls_by_rate[HISTORICAL_RATE])
            # Nested rate selection on the other folds (pooled test-like MCC).
            selected = np.zeros(len(score), dtype=bool)
            folds_present = np.unique(fold)
            for f in folds_present:
                other = fold != f
                if not other.any():
                    continue
                inner = {r: mcc(y_tl[other], calls_by_rate[r][other]) for r in rates}
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
                        selection_mcc_testlike=inner[rate],
                    )
                )
            nested_rows.append(
                decision_row(
                    frame, selected, split_id=split, model=model, fold=-1, selected_rate=None
                )
            )
            nested_calls[model] = (frame, selected)
            full = {r: mcc(y_tl, calls_by_rate[r]) for r in rates}
            deploy.setdefault(split, {})[model] = max(rates, key=lambda r: (full[r], -r))
        control_frame = nested_calls[CONTROL][0]
        for model in models:
            if model == CONTROL:
                continue
            for operating, calls_map in [("nested", nested_calls), ("fixed_6.53pct", fixed_calls)]:
                frame, calls = calls_map[model]
                _, control_calls = calls_map[CONTROL]
                assert frame["compound_id"].equals(control_frame["compound_id"])
                for label in ["y_testlike", "y_official"]:
                    boot_rows.append(
                        dict(
                            split_id=split,
                            model=model,
                            control=CONTROL,
                            operating_point=operating,
                            label_set=label.removeprefix("y_"),
                            **group_bootstrap(
                                frame,
                                calls,
                                control_calls,
                                label,
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

    # Blind diagnostics.
    sep15 = pl.read_parquet(source / inputs["sep15_payload"]).select(
        pl.col("Molecule_Name").alias("molecule_name"), pl.col(ENDPOINT).alias("sep15_call")
    )
    sep15 = sep15.join(names, on="molecule_name", how="left")
    assert sep15.height == 750 and sep15["compound_id"].null_count() == 0
    sep15_positive = set(sep15.filter(pl.col("sep15_call"))["compound_id"].to_list())
    blind_rows = []
    for model in blind_scores["model"].unique(maintain_order=True).to_list():
        frame = blind_scores.filter(pl.col("model") == model).sort(
            ["score", "compound_id"], descending=[True, False]
        )
        ids = frame["compound_id"].to_list()
        row = dict(
            model=model,
            sep15_positives=len(sep15_positive),
            top49_overlap_sep15=len(set(ids[:49]) & sep15_positive),
            p_active_gt_0_5=(
                int((frame["p_active"] > 0.5).sum())
                if frame["p_active"].null_count() == 0
                else None
            ),
            mean_p_active=frame["p_active"].mean(),
        )
        for split in cfg["split_strategies"]:
            if model not in deploy.get(split, {}):
                continue
            rate = deploy[split][model]
            k = int(round(rate * 750))
            row[f"{split}_rate"] = rate
            row[f"{split}_calls"] = k
            row[f"{split}_overlap_sep15"] = len(set(ids[:k]) & sep15_positive)
            # Score-threshold transfer: OOF k-th score per held fold, averaged.
            thresholds = []
            held = oof.filter((pl.col("split_id") == split) & (pl.col("model") == model))
            for f in held["fold"].unique().to_list():
                s = np.sort(held.filter(pl.col("fold") == f)["score"].to_numpy())[::-1]
                kk = int(round(rate * len(s)))
                if kk:
                    thresholds.append(s[kk - 1])
            threshold = float(np.mean(thresholds))
            row[f"{split}_oof_threshold"] = threshold
            row[f"{split}_threshold_calls"] = int((frame["score"] >= threshold).sum())
            nested_rates = nested.filter(
                (pl.col("split_id") == split) & (pl.col("model") == model) & (pl.col("fold") >= 0)
            )["selected_rate"].to_list()
            row[f"{split}_nested_rates"] = ",".join(f"{r:.2f}" for r in nested_rates)
        blind_rows.append(row)
    blind_table = pl.DataFrame(blind_rows, infer_schema_length=None)
    blind_table.write_csv(out / "blind_summary.csv")
    (out / "fit_info.json").write_text(json.dumps(fit_info, indent=2) + "\n")

    # Summary and report.
    summary = dict(
        config=cfg,
        endpoint=ENDPOINT,
        pool=pool_audit,
        smoke=bool(smoke),
        evaluated_development=len(development),
        legacy_comparator=legacy_audit,
        gate_choice=(
            "LightGBM regressor on D6 primary log2fc (all screened training molecules), "
            "Platt-calibrated to P(log2fc < gate_threshold) on inner group-fold OOF predictions"
        ),
        ranking=ranking_rows,
        decision_pooled=[r for r in curve_rows if r["fold"] == -1],
        nested=[r for r in nested_rows if r["fold"] == -1],
        nested_per_fold=[r for r in nested_rows if r["fold"] >= 0],
        bootstrap=boot_rows,
        deployment_rate=deploy,
        blind=blind_rows,
        runtime_seconds=round(time.time() - started, 1),
        caveats=[
            "Pseudo-negatives assume primary non-hits would be assigned negative",
            "Blind molecules were not primary screened; the test-like pool is an analogy",
            "Rates and models are evaluated on umap/scaffold families only",
            "Blind rank-based calls do not establish the blind scored-label mask",
        ],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")

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
            and c["rate"] == HISTORICAL_RATE
        )
        decision_main.append(
            dict(
                split_id=r["split_id"],
                model=r["model"],
                nested_mcc_testlike=r["mcc_testlike"],
                nested_mcc_official=r["mcc_official"],
                nested_calls=r["n_calls"],
                fixed653_mcc_testlike=fixed["mcc_testlike"],
                fixed653_mcc_official=fixed["mcc_official"],
                fixed653_calls_borderline=fixed["calls_borderline"],
                fixed653_calls_unscreened=fixed["calls_unscreened"],
            )
        )
    curve_wide = (
        curves.filter(pl.col("fold") == -1)
        .select("split_id", "model", "rate", "mcc_testlike")
        .with_columns(
            pl.col("rate").map_elements(lambda r: f"{100 * r:.2f}%", return_dtype=pl.String)
        )
        .pivot(on="rate", index=["split_id", "model"], values="mcc_testlike")
    )
    report = [
        "# D6 test-like pool evaluation",
        "",
        f"Development molecules evaluated: {len(development)}; smoke: {bool(smoke)}.",
        "",
        "## Pool counts",
        "",
        "```json",
        json.dumps(pool_audit, indent=2),
        "```",
        "",
        "## Ranking (pooled OOF)",
        "",
        markdown_table(
            ranking_rows,
            ["split_id", "model", "auroc_official", "ap_official", "auroc_testlike", "ap_testlike"],
        ),
        "",
        "## Decision (pooled held-fold calls)",
        "",
        markdown_table(decision_main, list(decision_main[0].keys())),
        "",
        "## Test-like MCC by call rate (pooled)",
        "",
        markdown_table(curve_wide.to_dicts(), curve_wide.columns),
        "",
        "## Paired group bootstrap versus official_lgbm (test-like labels)",
        "",
        markdown_table(
            [r for r in boot_rows if r["label_set"] == "testlike"],
            ["split_id", "model", "operating_point", "delta", "ci_low", "ci_high", "p_delta_gt_0"],
        ),
        "",
        "## Blind diagnostics",
        "",
        markdown_table(blind_rows, list(blind_rows[0].keys())),
        "",
        f"Legacy comparator audit: `{json.dumps(legacy_audit)}`",
        "",
    ]
    (out / "report.md").write_text("\n".join(report))
    print(json.dumps(dict(runtime_seconds=summary["runtime_seconds"])), flush=True)


if __name__ == "__main__":
    main()
