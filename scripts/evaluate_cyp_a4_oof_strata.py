"""Where does the submitted A4 equal-seven ensemble get its score?

Stratifies the complete-pool UMAP out-of-fold predictions of the submitted CYP3A4 TDI
equal-seven ensemble by primary-screen activity, label provenance and chemical
proximity to the blinded set, then reweights the labelled rows towards the blind
population with two independent estimators of the blind-side stratum mix:

1. `nn_transfer`  - each blind molecule inherits the stratum of its nearest
   development molecule by Morgan Tanimoto (k = 1 and a similarity-weighted k = 5).
2. `gate`         - a fold-excluded LightGBM regressor of CYP3A4 primary log2fc,
   Platt-calibrated to P(log2fc < tau) for tau in {-1, -2}, applied to the blind
   molecules; both sides are binned with the *predicted* log2fc so no stratum is
   empty on either side.

No new classification model is fitted: the ensemble scores are consumed frozen. The
log2fc regressor is fitted strictly inside the training part of every outer fold.
"""

import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl
from cyp_primary_supervision import load_primary
from joblib import Parallel, delayed
from lightgbm import LGBMRegressor
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

PRIMARY_ORDER = ["CYP2D6", "CYP3A4", "CYP1A2", "CYP2C9"]


def featurize(smiles):
    """Morgan r2/2048 bits plus RDKit 2D descriptors, and the bit vectors themselves."""
    RDLogger.DisableLog("rdApp.*")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = np.zeros((len(smiles), 2048), dtype=np.float32)
    desc = np.full((len(smiles), len(Descriptors._descList)), np.nan, dtype=np.float32)
    bits = []
    for i, text in enumerate(smiles):
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            raise ValueError(f"Unparseable canonical molecule: {text}")
        vector = generator.GetFingerprint(mol)
        bits.append(vector)
        DataStructs.ConvertToNumpyArray(vector, fp[i])
        for j, (_, function) in enumerate(Descriptors._descList):
            try:
                desc[i, j] = function(mol)
            except (ValueError, RuntimeError, OverflowError, ZeroDivisionError):
                pass
    desc[~np.isfinite(desc)] = np.nan
    return fp, desc, bits


def lgbm_params(cfg, seed):
    return dict(
        **cfg["lgbm"],
        verbosity=-1,
        n_jobs=cfg["lgbm_threads"],
        random_state=seed,
        deterministic=True,
        force_col_wise=True,
    )


def fit_regressor(x, rows, target, predict, seeds, cfg):
    predictions = []
    for seed in seeds:
        model = LGBMRegressor(**lgbm_params(cfg, seed))
        model.fit(x[rows], target[rows])
        predictions.append([model.predict(x[p]) for p in predict])
    return [np.mean([s[i] for s in predictions], axis=0) for i in range(len(predict))]


def fit_gate(task, x, screened, log2fc, groups, predict, cfg):
    """Predicted log2fc plus Platt-calibrated P(log2fc < tau) for each threshold."""
    started = time.time()
    seeds = [cfg["seed"] + s for s in range(cfg["lgbm_seeds"])]
    inner = np.full(len(screened), np.nan)
    for fit, valid in GroupKFold(n_splits=cfg["inner_folds"]).split(
        screened, groups=groups[screened]
    ):
        inner[valid] = fit_regressor(x, screened[fit], log2fc, [screened[valid]], seeds, cfg)[0]
    predicted = fit_regressor(x, screened, log2fc, predict, seeds, cfg)
    calibrated, info = (
        {},
        dict(
            task=task,
            screened_train=len(screened),
            inner_spearman=float(
                pl.DataFrame({"a": inner, "b": log2fc[screened]})
                .select(pl.corr("a", "b", method="spearman"))
                .item()
            ),
            seconds=None,
        ),
    )
    for tau in cfg["gate_thresholds"]:
        active = (log2fc[screened] < tau).astype(int)
        model = LogisticRegression(C=1e4, max_iter=1000).fit(inner[:, None], active)
        calibrated[tau] = [model.predict_proba(p[:, None])[:, 1] for p in predicted]
        info[f"active_train_{tau}"] = int(active.sum())
        info[f"calibration_coef_{tau}"] = float(model.coef_[0, 0])
    info["seconds"] = round(time.time() - started, 1)
    print(json.dumps(info), flush=True)
    return task, predicted, calibrated, info


def bin_log2fc(values, cuts):
    """Ordered primary-activity bins; NaN becomes the explicit `none` level.

    With cuts [-3, -2, -1, -0.5] the levels are `<= -3`, `(-3, -2]`, `(-2, -1]`,
    `(-1, -0.5)`, `>= -0.5` and `none`, so the last level reproduces the
    "weak or no primary inhibition" population used in the pool audit.
    """
    strong, mid, weak, edge = cuts
    labels = [
        f"<= {strong:g}",
        f"({strong:g}, {mid:g}]",
        f"({mid:g}, {weak:g}]",
        f"({weak:g}, {edge:g})",
        f">= {edge:g}",
        "none",
    ]
    out = np.full(len(values), "none", dtype=object)
    finite = np.isfinite(values)
    out[finite & (values <= strong)] = labels[0]
    out[finite & (values > strong) & (values <= mid)] = labels[1]
    out[finite & (values > mid) & (values <= weak)] = labels[2]
    out[finite & (values > weak) & (values < edge)] = labels[3]
    out[finite & (values >= edge)] = labels[4]
    assert set(out) <= set(labels)
    return out, labels


def quartile_labels(values, edges):
    out = np.full(len(values), "Q4", dtype=object)
    out[values <= edges[2]] = "Q3"
    out[values <= edges[1]] = "Q2"
    out[values <= edges[0]] = "Q1"
    return out


def mcc_from_rates(tpr, fpr, prevalence):
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator > 0 else 0.0


def weighted_rates(y, calls, weights):
    positive, negative = y == 1, y == 0
    tpr = float(np.average(calls[positive], weights=weights[positive]))
    fpr = float(np.average(calls[negative], weights=weights[negative]))
    return tpr, fpr


def top_calls(score, fold, count_per_750):
    """Call the top round(rate * n) of each held fold; rows with fold < 0 are excluded."""
    calls = np.zeros(len(score), dtype=bool)
    rate = count_per_750 / 750.0
    for f in np.unique(fold[fold >= 0]):
        rows = np.flatnonzero(fold == f)
        k = int(round(rate * len(rows)))
        calls[rows[np.argsort(-score[rows], kind="stable")[:k]]] = True
    return calls


def fmt(value, digits=3):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "–"
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
    compounds = pl.read_parquet(source / inputs["compounds"]).sort("compound_id")
    aliases = pl.read_parquet(source / inputs["aliases"])
    test_mapping = pl.read_parquet(source / inputs["test_mapping"])
    measured = pl.read_parquet(source / inputs["measurements"])
    assignments = pl.read_parquet(source / inputs["assignments"])

    labels = measured.filter(pl.col("endpoint") == cfg["endpoint"]).select(
        "compound_id", pl.col("value").alias("label")
    )
    pic50 = measured.filter(pl.col("endpoint") == cfg["pic50_endpoint"]).select(
        "compound_id", pl.col("value").alias("pic50")
    )
    table = (
        compounds.join(labels, on="compound_id", how="left", maintain_order="left")
        .join(pic50, on="compound_id", how="left", maintain_order="left")
        .select("compound_id", "standardized_smiles", "label", "pic50")
    )
    primary, _ = load_primary(source / inputs["primary_screen"], aliases, table)
    log2fc = primary[:, PRIMARY_ORDER.index(cfg["primary_enzyme"])]

    blind_ids = test_mapping["compound_id"].unique().to_list()
    is_blind = table["compound_id"].is_in(blind_ids).to_numpy()
    blind = np.flatnonzero(is_blind)
    development = np.flatnonzero(~is_blind)
    y = table["label"].to_numpy().astype(float)
    assert len(blind) == 750 and not np.isfinite(y[blind]).any()
    assert not np.isfinite(log2fc[blind]).any()

    # Frozen submitted ensemble OOF on the complete UMAP pool.
    oof = pl.read_parquet(source / inputs["equal7_oof"])
    assert oof["model"].unique().to_list() == [cfg["oof_model"]]
    assert oof["split_id"].unique().to_list() == [cfg["split_strategy"]]
    assert oof.height == len(development), (oof.height, len(development))
    aligned = table.select("compound_id").join(
        oof.select("compound_id", "fold", "group", "y_true", "score"),
        on="compound_id",
        how="left",
        validate="1:1",
        maintain_order="left",
    )
    score = aligned["score"].fill_null(np.nan).to_numpy().astype(float)
    fold = aligned["fold"].fill_null(-1).to_numpy().astype(int)
    groups = aligned["group"].fill_null("blind").to_numpy().astype(str)
    oof_y = aligned["y_true"].fill_null(np.nan).to_numpy().astype(float)
    assert np.array_equal(np.isfinite(oof_y), np.isfinite(y))
    assert np.array_equal(oof_y[np.isfinite(oof_y)], y[np.isfinite(y)])
    assert np.isfinite(score[development]).all()
    assert (fold[blind] == -1).all() and (fold[development] >= 0).all()
    declared = table.select("compound_id").join(
        assignments.filter(pl.col("strategy") == cfg["split_strategy"]).select(
            "compound_id", "fold", "group"
        ),
        on="compound_id",
        how="left",
        maintain_order="left",
    )
    assert np.array_equal(declared["fold"].fill_null(-1).to_numpy().astype(int), fold)
    assert np.array_equal(declared["group"].fill_null("blind").to_numpy().astype(str), groups)
    labelled = np.flatnonzero(np.isfinite(y))
    assert set(labelled) <= set(development)

    # Features and nearest-neighbour similarities.
    smiles = table["standardized_smiles"].to_list()
    chunks = np.array_split(np.arange(len(smiles)), cfg["feature_jobs"])
    parts = Parallel(n_jobs=cfg["feature_jobs"])(
        delayed(featurize)([smiles[i] for i in chunk]) for chunk in chunks
    )
    fp = np.concatenate([a for a, _, _ in parts])
    desc = np.concatenate([b for _, b, _ in parts])
    bits = [v for _, _, part in parts for v in part]
    x_morgan = np.concatenate([fp, desc], axis=1)
    blind_bits = [bits[i] for i in blind]
    dev_bits = [bits[i] for i in development]
    nearest_blind = np.zeros(len(table))
    for i in range(len(table)):
        nearest_blind[i] = max(DataStructs.BulkTanimotoSimilarity(bits[i], blind_bits))
    nn_top, nearest_dev_of_blind = [], np.zeros(len(blind))
    for j, b in enumerate(blind):
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(bits[b], dev_bits))
        order = np.argsort(-sims, kind="stable")[: max(cfg["nn_transfer_k"])]
        nn_top.append((development[order], sims[order]))
        nearest_dev_of_blind[j] = float(sims[order[0]])
    print(
        json.dumps(
            dict(
                feature_seconds=round(time.time() - started, 1),
                nearest_blind_mean_labelled=float(nearest_blind[labelled].mean()),
                nearest_dev_of_blind_mean=float(nearest_dev_of_blind.mean()),
            )
        ),
        flush=True,
    )

    # Fold-excluded log2fc gate: five outer folds plus a development-wide blind refit.
    tasks = []
    for f in range(cfg["outer_folds"]):
        train = development[fold[development] != f]
        held = development[fold[development] == f]
        assert set(groups[train]).isdisjoint(groups[held])
        tasks.append((("fold", f), train, [held]))
    tasks.append((("blind", -1), development, [blind]))
    gate_outputs = Parallel(n_jobs=cfg["task_jobs"])(
        delayed(fit_gate)(
            key,
            x_morgan,
            train[np.isfinite(log2fc[train])],
            log2fc,
            groups,
            predict,
            cfg,
        )
        for key, train, predict in tasks
    )
    predicted_log2fc = np.full(len(table), np.nan)
    p_active = {tau: np.full(len(table), np.nan) for tau in cfg["gate_thresholds"]}
    held_rows = {key: predict[0] for key, _, predict in tasks}
    gate_info = []
    for key, predicted, calibrated, info in gate_outputs:
        rows = held_rows[key]
        predicted_log2fc[rows] = predicted[0]
        for tau in cfg["gate_thresholds"]:
            p_active[tau][rows] = calibrated[tau][0]
        gate_info.append(info)
    assert np.isfinite(predicted_log2fc).all()

    # Stratifications.
    observed_bin, bin_labels = bin_log2fc(log2fc, cfg["log2fc_bins"])
    predicted_bin, _ = bin_log2fc(predicted_log2fc, cfg["log2fc_bins"])
    provenance = np.where(
        np.isfinite(table["pic50"].to_numpy().astype(float)),
        "direct_pic50",
        "assigned_negative",
    ).astype(object)
    provenance[~np.isfinite(y)] = "unlabelled"
    edges = np.quantile(nearest_blind[labelled], [0.25, 0.5, 0.75])
    similarity_bin = quartile_labels(nearest_blind, edges)
    strata = {
        "log2fc_bin": (observed_bin, [*bin_labels]),
        "pic50_provenance": (provenance, ["direct_pic50", "assigned_negative"]),
        "nearest_blind_quartile": (similarity_bin, ["Q1", "Q2", "Q3", "Q4"]),
        "predicted_log2fc_bin": (predicted_bin, [b for b in bin_labels if b != "none"]),
    }
    for tau in cfg["gate_thresholds"]:
        strata[f"gate_{tau:g}"] = (
            np.where(p_active[tau] >= 0.5, "active", "inactive").astype(object),
            ["active", "inactive"],
        )

    # Per-stratum ranking and 249/750 decision contribution.
    reference = cfg["reference_count"]
    calls_by_count = {c: top_calls(score, fold, c) for c in cfg["counts_per_750"]}
    stratum_rows = []
    for name, (values, levels) in strata.items():
        for level in levels:
            rows = labelled[values[labelled] == level]
            if len(rows) == 0:
                stratum_rows.append(dict(stratification=name, stratum=level, n=0))
                continue
            yy, ss = y[rows].astype(int), score[rows]
            calls = calls_by_count[reference][rows]
            tp = int(((yy == 1) & calls).sum())
            stratum_rows.append(
                dict(
                    stratification=name,
                    stratum=level,
                    n=len(rows),
                    positives=int(yy.sum()),
                    positive_rate=float(yy.mean()),
                    auroc=float(roc_auc_score(yy, ss)) if len(np.unique(yy)) > 1 else None,
                    ap=float(average_precision_score(yy, ss)) if len(np.unique(yy)) > 1 else None,
                    mean_score=float(ss.mean()),
                    calls_at_249=int(calls.sum()),
                    tp_at_249=tp,
                    precision_at_249=tp / calls.sum() if calls.sum() else None,
                    recall_at_249=tp / yy.sum() if yy.sum() else None,
                    share_of_labelled=len(rows) / len(labelled),
                    share_of_calls_at_249=(
                        calls.sum() / calls_by_count[reference][labelled].sum()
                        if calls_by_count[reference][labelled].sum()
                        else None
                    ),
                    mean_nearest_blind=float(nearest_blind[rows].mean()),
                )
            )
    strata_table = pl.DataFrame(stratum_rows, infer_schema_length=None)
    strata_table.write_csv(out / "strata.csv")

    # Blind-side stratum mixes.
    def nn_mix(values, levels, k):
        counts = {level: 0.0 for level in levels}
        total = 0.0
        for neighbours, sims in nn_top:
            take = neighbours[:k]
            weights = np.asarray(sims[:k], dtype=float)
            weights = weights / weights.sum() if weights.sum() > 0 else np.ones(k) / k
            for n, w in zip(take, weights, strict=True):
                level = values[n]
                if level in counts:
                    counts[level] += float(w)
                    total += float(w)
        return {level: (counts[level] / total if total else 0.0) for level in levels}

    weightings, mix_rows = {}, []
    for name, (values, levels) in strata.items():
        labelled_mix = {level: float((values[labelled] == level).mean()) for level in levels}
        assert set(values[labelled]) <= set(levels), (name, set(values[labelled]) - set(levels))
        if name.startswith("predicted_log2fc") or name.startswith("gate_"):
            blind_mix = {level: float((values[blind] == level).mean()) for level in levels}
            weightings[f"{name}:direct"] = (name, blind_mix)
            mix_rows.append(dict(estimator=f"{name}:direct", stratification=name, **blind_mix))
        for k in cfg["nn_transfer_k"]:
            blind_mix = nn_mix(values, levels, k)
            weightings[f"{name}:nn{k}"] = (name, blind_mix)
            mix_rows.append(dict(estimator=f"{name}:nn{k}", stratification=name, **blind_mix))
        mix_rows.append(dict(estimator=f"{name}:labelled", stratification=name, **labelled_mix))
    pl.DataFrame(mix_rows, infer_schema_length=None).write_csv(out / "stratum_mix.csv")

    # Reweighted population metrics and rate curves.
    yy = y[labelled].astype(int)
    ss = score[labelled]
    reweight_rows, curve_rows = [], []
    base = np.ones(len(labelled))
    cases = [("none", None, base)]
    for key, (name, blind_mix) in weightings.items():
        values = strata[name][0][labelled]
        labelled_mix = {level: float((values == level).mean()) for level in strata[name][1]}
        weights = np.array(
            [
                blind_mix.get(level, 0.0) / labelled_mix[level] if labelled_mix[level] > 0 else 0.0
                for level in values
            ]
        )
        cases.append((key, name, weights))
    for key, name, weights in cases:
        if weights.sum() <= 0 or weights[yy == 1].sum() <= 0 or weights[yy == 0].sum() <= 0:
            reweight_rows.append(dict(weighting=key, stratification=name, usable=False))
            continue
        prevalence = float(np.average(yy, weights=weights))
        row = dict(
            weighting=key,
            stratification=name,
            usable=True,
            effective_n=float(weights.sum() ** 2 / np.square(weights).sum()),
            projected_prevalence=prevalence,
            auroc=float(roc_auc_score(yy, ss, sample_weight=weights)),
            ap=float(average_precision_score(yy, ss, sample_weight=weights)),
        )
        for count in cfg["counts_per_750"]:
            calls = calls_by_count[count][labelled]
            tpr, fpr = weighted_rates(yy, calls, weights)
            projected = mcc_from_rates(tpr, fpr, prevalence)
            fixed = mcc_from_rates(tpr, fpr, cfg["fixed_prior"])
            row[f"mcc_{count}"] = projected
            curve_rows.append(
                dict(
                    weighting=key,
                    stratification=name,
                    count_per_750=count,
                    all_pool_calls=int(calls_by_count[count].sum()),
                    scored_calls=int(calls.sum()),
                    projected_prevalence=prevalence,
                    tpr=tpr,
                    fpr=fpr,
                    projected_mcc=projected,
                    mcc_at_fixed_prior=fixed,
                )
            )
        reweight_rows.append(row)
    reweighted = pl.DataFrame(reweight_rows, infer_schema_length=None)
    curves = pl.DataFrame(curve_rows, infer_schema_length=None)
    reweighted.write_csv(out / "reweighted_metrics.csv")
    curves.write_csv(out / "reweighted_rate_curves.csv")
    best = (
        curves.sort("projected_mcc", descending=True)
        .group_by("weighting", maintain_order=True)
        .first()
        .select("weighting", "count_per_750", "projected_mcc", "projected_prevalence")
    )

    # Blind-side gate distribution and the submitted vector check.
    submission = pl.read_parquet(source / inputs["submission"])
    names = test_mapping.select("compound_id", "molecule_name").unique()
    submitted = (
        submission.select(
            pl.col("Molecule_Name").alias("molecule_name"),
            pl.col(cfg["endpoint"]).alias("call"),
        )
        .join(names, on="molecule_name", how="left")
        .filter(pl.col("call"))
    )
    blind_scores = pl.read_parquet(source / inputs["equal7_blind"])
    blind_frame = (
        pl.DataFrame({"compound_id": table["compound_id"].gather(blind)})
        .join(blind_scores.select("compound_id", "a4_prediction"), on="compound_id", how="left")
        .with_columns(
            pl.Series("predicted_log2fc", predicted_log2fc[blind]),
            *[
                pl.Series(f"p_active_{tau:g}", p_active[tau][blind])
                for tau in cfg["gate_thresholds"]
            ],
            pl.Series("nearest_dev_tanimoto", nearest_dev_of_blind),
        )
    )
    assert blind_frame["a4_prediction"].null_count() == 0
    submitted_ids = set(submitted["compound_id"].to_list())
    top_reference = set(
        blind_frame.sort("a4_prediction", descending=True).head(reference)["compound_id"].to_list()
    )
    blind_check = dict(
        submitted_positives=len(submitted_ids),
        top_reference_matches_submission=bool(top_reference == submitted_ids),
        overlap=len(top_reference & submitted_ids),
    )
    assert blind_check["top_reference_matches_submission"], blind_check
    blind_frame = blind_frame.with_columns(
        pl.col("compound_id").is_in(list(submitted_ids)).alias("submitted_call")
    )
    blind_frame.write_csv(out / "blind_gate.csv")
    blind_gate_summary = (
        blind_frame.group_by("submitted_call")
        .agg(
            pl.len().alias("n"),
            pl.col("predicted_log2fc").mean().alias("mean_predicted_log2fc"),
            *[
                pl.col(f"p_active_{tau:g}").mean().alias(f"mean_p_active_{tau:g}")
                for tau in cfg["gate_thresholds"]
            ],
            pl.col("nearest_dev_tanimoto").mean().alias("mean_nearest_dev_tanimoto"),
        )
        .sort("submitted_call", descending=True)
    )

    summary = dict(
        config=cfg,
        pool=dict(
            molecules=len(table),
            development=len(development),
            blind=len(blind),
            labelled=len(labelled),
            labelled_positive=int(np.nansum(y)),
            labelled_positive_rate=float(np.nanmean(y[labelled])),
            labelled_with_direct_pic50=int((provenance[labelled] == "direct_pic50").sum()),
            labelled_assigned_negative=int((provenance[labelled] == "assigned_negative").sum()),
            labelled_with_primary=int(np.isfinite(log2fc[labelled]).sum()),
            labelled_positive_max_log2fc=float(np.nanmax(log2fc[labelled][y[labelled] == 1])),
            screened_development=int(np.isfinite(log2fc[development]).sum()),
            unlabelled_primary_hits=int(
                (
                    ~np.isfinite(y[development])
                    & np.isfinite(log2fc[development])
                    & (log2fc[development] < -0.5)
                ).sum()
            ),
            unlabelled_primary_non_hits=int(
                (
                    ~np.isfinite(y[development])
                    & np.isfinite(log2fc[development])
                    & (log2fc[development] >= -0.5)
                ).sum()
            ),
            nearest_blind_quartile_edges=[float(e) for e in edges],
        ),
        gate=gate_info,
        strata=stratum_rows,
        stratum_mix=mix_rows,
        reweighted=reweight_rows,
        rate_curves=curve_rows,
        best_count_by_weighting=best.to_dicts(),
        blind_check=blind_check,
        blind_gate_summary=blind_gate_summary.to_dicts(),
        runtime_seconds=round(time.time() - started, 1),
        caveats=[
            "The ensemble OOF is consumed frozen; small TabICL/TabPFN reload differences "
            "documented for run e74748cd are inherited",
            "Blind stratum mixes are estimates from training chemistry, not measurements",
            "The gate reweighting bins both sides by predicted log2fc, so it measures "
            "predicted-activity transfer rather than true activity transfer",
            "Importance weights cannot correct for a stratum the labelled pool does not cover",
        ],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")

    report = [
        "# A4 equal-seven OOF stratification and blind reweighting",
        "",
        "## Pool",
        "",
        "```json",
        json.dumps(summary["pool"], indent=2),
        "```",
        "",
        "## Strata (labelled rows, 249/750 call rate on the complete pool)",
        "",
        markdown_table(
            stratum_rows,
            [
                "stratification",
                "stratum",
                "n",
                "positive_rate",
                "auroc",
                "ap",
                "calls_at_249",
                "tp_at_249",
                "precision_at_249",
                "mean_nearest_blind",
            ],
        ),
        "",
        "## Reweighted population metrics",
        "",
        markdown_table(
            reweight_rows,
            [
                "weighting",
                "usable",
                "effective_n",
                "projected_prevalence",
                "auroc",
                "ap",
                *[f"mcc_{c}" for c in cfg["counts_per_750"]],
            ],
        ),
        "",
        "## Best count per weighting",
        "",
        markdown_table(best.to_dicts(), best.columns),
        "",
        "## Blind gate by submitted call",
        "",
        markdown_table(blind_gate_summary.to_dicts(), blind_gate_summary.columns),
        "",
        f"Submitted-vector check: `{json.dumps(blind_check)}`",
        "",
    ]
    (out / "report.md").write_text("\n".join(report))
    print(
        json.dumps(
            dict(
                runtime_seconds=summary["runtime_seconds"],
                blind_check=blind_check,
                best=best.to_dicts(),
            )
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
