"""Induce the D6 label-coverage defect on 1A2/2C9/3A4, then test the pseudo-negative fix.

D6 cannot test the pseudo-negative device, because D6 has no labels in the region the
device relabels. CYP1A2, CYP2C9 and CYP3A4 do. This script therefore manufactures the
D6 situation on those three enzymes and judges the repair against labels that neither
the defect arm nor the fix arm ever saw.

Three matched arms share features, folds and seeds:

- `defect`   fits on primary hits only (log2fc < borderline_min), so weak-primary and
             missing-primary labels are hidden from fitting and from rate selection,
             mimicking D6's coverage.
- `fixed`    adds the hidden weak-primary rows back as pseudo-negatives (label False,
             never their real value).
- `fixed_full` is the D6 device verbatim: every weak-primary row, labelled or not, enters
             as a pseudo-negative, which matches the dose D6 receives (2,693 added rows
             against 1,495 visible ones) rather than only the labelled share.
- `ceiling_screened` fits on the real labels of every screened row, so it bounds the part
             of the gap that a weak-primary repair could possibly close.
- `ceiling`  fits on the real, complete label set, including missing-primary rows.
- `fixed_missing` (only where a missing-primary stratum exists) additionally relabels
             the missing-primary rows False, a deliberately wrong variant whose failure
             checks that the metric responds to label quality.

Every arm is judged on the real, complete labels of the held-out fold. Calls are taken
from all development molecules of the held fold; metrics use the labelled subsets.
Nothing here touches the blind set. No candidate, prepared submission or upload.
"""

import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl
from cyp_primary_supervision import load_primary
from joblib import Parallel, delayed
from lightgbm import LGBMClassifier
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sklearn.metrics import average_precision_score, roc_auc_score

DEFECT = "defect"
FIXED = "fixed"
FIXED_FULL = "fixed_full"
CEILING_SCREENED = "ceiling_screened"
CEILING = "ceiling"
FIXED_MISSING = "fixed_missing"
PRIMARY_ORDER = ["CYP2D6", "CYP3A4", "CYP1A2", "CYP2C9"]


def featurize(smiles):
    """Fresh Morgan r2/2048 bits and RDKit 2D descriptors (same recipe as the D6 device)."""
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


def mixture_mcc(y, calls, strata_masks, weights):
    """MCC after reweighting the scored strata to a target composition.

    The labelled strong-primary and weak-primary rows are samples of two populations
    whose mixture the blind scored set is believed to be. Reweighting them to a target
    mixture fraction keeps every label real while matching the composition that matters.
    """
    tp = fp = fn = tn = 0.0
    effective_numerator = effective_denominator = 0.0
    for name, weight in weights.items():
        mask = strata_masks[name]
        n = int(mask.sum())
        if n == 0 or weight == 0:
            continue
        scale = weight / n
        labels = y[mask].astype(bool)
        called = calls[mask]
        tp += scale * float((labels & called).sum())
        fp += scale * float((~labels & called).sum())
        fn += scale * float((labels & ~called).sum())
        tn += scale * float((~labels & ~called).sum())
        effective_numerator += weight
        effective_denominator += weight**2 / n
    effective = effective_numerator**2 / effective_denominator if effective_denominator else 0.0
    return float(mcc_from_counts(tp, fp, fn, tn)), effective


def top_calls(score, rate):
    calls = np.zeros(len(score), dtype=bool)
    k = int(round(rate * len(score)))
    calls[np.argsort(-score, kind="stable")[:k]] = True
    return calls


def fit_arm(x, rows, y, predict, seeds, cfg):
    """Balanced LightGBM averaged over seeds; `y` holds the arm's visible labels."""
    scores = []
    for seed in seeds:
        model = LGBMClassifier(
            class_weight="balanced",
            **cfg["lgbm"],
            verbosity=-1,
            n_jobs=cfg["lgbm_threads"],
            random_state=seed,
            deterministic=True,
            force_col_wise=True,
        )
        model.fit(x[rows], y[rows].astype(int))
        scores.append(model.predict_proba(x[predict])[:, 1])
    return np.mean(scores, axis=0)


def fit_fold(enzyme, split, fold, train, held, arms, x, cfg):
    """Fit every arm of one (enzyme, split, fold) cell and score the held fold."""
    started = time.time()
    seeds = [cfg["seed"] + s for s in range(cfg["lgbm_seeds"])]
    scores, info = {}, dict(enzyme=enzyme, split_id=split, fold=fold, train=len(train))
    for name, labels in arms.items():
        rows = train[np.isfinite(labels[train])]
        scores[name] = fit_arm(x, rows, labels, held, seeds, cfg)
        info[f"{name}_rows"] = len(rows)
        info[f"{name}_positive"] = int(np.nansum(labels[rows]))
    info["seconds"] = round(time.time() - started, 1)
    print(json.dumps(info), flush=True)
    return enzyme, split, fold, scores, info


def group_bootstrap(y, groups, calls_a, calls_b, samples, seed, row_weights=None):
    """Paired group bootstrap of pooled MCC(a) - MCC(b), optionally stratum-reweighted."""
    keep = np.isfinite(y)
    if row_weights is not None:
        keep = keep & (row_weights > 0)
    y = y[keep].astype(bool)
    groups = groups[keep]
    weight = np.ones(int(keep.sum())) if row_weights is None else row_weights[keep]
    _, index = np.unique(groups, return_inverse=True)
    count = index.max() + 1
    counts = []
    for calls in (calls_a[keep], calls_b[keep]):
        counts.append(
            [
                np.bincount(index, weights=(y & calls) * weight, minlength=count),
                np.bincount(index, weights=(~y & calls) * weight, minlength=count),
                np.bincount(index, weights=(y & ~calls) * weight, minlength=count),
                np.bincount(index, weights=(~y & ~calls) * weight, minlength=count),
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


def fmt(value, digits=3):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "–"
    return f"{value:.{digits}f}" if isinstance(value, float) else str(value)


def markdown_table(rows, columns):
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(c)) for c in columns) + " |")
    return "\n".join(lines)


def build_labels(y_real, log2fc, development, cfg):
    """Visible label vectors for each arm, plus the scoring strata of the real labels."""
    labelled = np.isfinite(y_real)
    screened = np.isfinite(log2fc)
    in_development = np.zeros(len(y_real), dtype=bool)
    in_development[development] = True
    cut, border = cfg["pseudo_negative_min_log2fc"], cfg["borderline_min_log2fc"]
    strata = dict(
        strong=labelled & screened & (log2fc < border),
        borderline=labelled & screened & (log2fc >= border) & (log2fc < cut),
        weak=labelled & screened & (log2fc >= cut),
        missing=labelled & ~screened,
    )
    # The D6 device relabels every weak-primary row, including the unlabelled ones.
    weak_any = in_development & screened & (log2fc >= cut)
    visible_strong = np.where(strata["strong"], y_real, np.nan)
    arms = {
        DEFECT: visible_strong,
        FIXED: np.where(strata["weak"], 0.0, visible_strong),
        FIXED_FULL: np.where(weak_any, 0.0, visible_strong),
        CEILING_SCREENED: np.where(labelled & screened, y_real, np.nan),
        CEILING: np.where(labelled, y_real, np.nan),
    }
    if strata["missing"].any():
        arms[FIXED_MISSING] = np.where((in_development & ~screened) | weak_any, 0.0, visible_strong)
    strata["weak_any_unlabelled"] = weak_any & ~labelled
    return arms, strata


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
    primary, _ = load_primary(source / inputs["primary_screen"], aliases, compounds)
    is_blind = compounds["compound_id"].is_in(test_mapping["compound_id"].to_list()).to_numpy()
    assert is_blind.sum() == 750
    development = np.flatnonzero(~is_blind)
    assert set(assignments["compound_id"].unique().to_list()) == set(
        compounds["compound_id"][development].to_list()
    )

    # Labels, primary screen and induced-defect masks per enzyme.
    enzyme_data, census = {}, []
    for enzyme in cfg["enzymes"]:
        y_real = (
            compounds.select("compound_id")
            .join(
                measured.filter(pl.col("endpoint") == f"{enzyme}_is_TDI").select(
                    "compound_id", pl.col("value").alias("y")
                ),
                on="compound_id",
                how="left",
                maintain_order="left",
            )["y"]
            .to_numpy()
            .astype(float)
        )
        log2fc = primary[:, PRIMARY_ORDER.index(enzyme)]
        assert not np.isfinite(y_real[is_blind]).any()
        assert not np.isfinite(log2fc[is_blind]).any()
        arms, strata = build_labels(y_real, log2fc, development, cfg)
        enzyme_data[enzyme] = dict(y_real=y_real, log2fc=log2fc, arms=arms, strata=strata)
        row = dict(
            enzyme=enzyme,
            labelled=int(np.isfinite(y_real).sum()),
            positives=int(np.nansum(y_real)),
            arms=",".join(arms),
        )
        for name, mask in strata.items():
            row[f"{name}_n"] = int(mask.sum())
            row[f"{name}_positives"] = int(np.nansum(y_real[mask]))
        row["device_dose"] = (row["weak_n"] + row["weak_any_unlabelled_n"]) / max(
            row["strong_n"], 1
        )
        row["hidden_from_defect"] = int(
            (strata["weak"] | strata["borderline"] | strata["missing"]).sum()
        )
        row["hidden_positives"] = int(
            np.nansum(y_real[strata["weak"] | strata["borderline"] | strata["missing"]])
        )
        census.append(row)
    print(json.dumps(census), flush=True)
    pl.DataFrame(census).write_csv(out / "label_census.csv")

    # Features.
    smiles = compounds["standardized_smiles"].to_list()
    chunks = np.array_split(np.arange(len(smiles)), cfg["feature_jobs"])
    parts = Parallel(n_jobs=cfg["feature_jobs"])(
        delayed(featurize)([smiles[i] for i in chunk]) for chunk in chunks
    )
    x = np.concatenate(
        [
            np.concatenate([a for a, _ in parts]),
            np.concatenate([b for _, b in parts]),
        ],
        axis=1,
    )
    print(
        json.dumps(dict(features=x.shape[1], seconds=round(time.time() - started, 1))), flush=True
    )

    # Outer tasks.
    family = {}
    tasks = []
    for split in cfg["split_strategies"]:
        ordered = compounds.select("compound_id").join(
            assignments.filter(pl.col("strategy") == split).select("compound_id", "fold", "group"),
            on="compound_id",
            how="left",
            maintain_order="left",
        )
        fold = ordered["fold"].fill_null(-1).to_numpy()
        groups = ordered["group"].fill_null("blind").to_numpy().astype(str)
        family[split] = (fold, groups)
        for f in range(cfg["outer_folds"]):
            train = development[fold[development] != f]
            held = development[fold[development] == f]
            assert set(groups[train]).isdisjoint(groups[held])
            for enzyme in cfg["enzymes"]:
                tasks.append((enzyme, split, int(f), train, held))
    outputs = Parallel(n_jobs=cfg["task_jobs"])(
        delayed(fit_fold)(enzyme, split, f, train, held, enzyme_data[enzyme]["arms"], x, cfg)
        for enzyme, split, f, train, held in tasks
    )
    fit_info = [info for *_, info in outputs]
    (out / "fit_info.json").write_text(json.dumps(fit_info, indent=2) + "\n")

    # Assemble held-fold scores keyed by molecule.
    held_rows = {(e, s, f): held for e, s, f, _, held in tasks}
    oof_parts = []
    for enzyme, split, f, scores, _ in outputs:
        rows = held_rows[(enzyme, split, f)]
        data = enzyme_data[enzyme]
        for arm, score in scores.items():
            oof_parts.append(
                pl.DataFrame(
                    dict(
                        compound_id=compounds["compound_id"].gather(rows),
                        endpoint=[f"{enzyme}_is_TDI"] * len(rows),
                        enzyme=[enzyme] * len(rows),
                        split_id=[split] * len(rows),
                        fold=[f] * len(rows),
                        arm=[arm] * len(rows),
                        score=np.asarray(score, dtype=np.float64),
                        y_real=data["y_real"][rows],
                        y_visible=data["arms"][arm][rows],
                        log2fc=data["log2fc"][rows],
                        stratum=np.select(
                            [data["strata"][k][rows] for k in ["strong", "borderline", "weak"]],
                            ["strong", "borderline", "weak"],
                            default=np.where(
                                data["strata"]["missing"][rows], "missing", "unlabelled"
                            ),
                        ),
                        group=family[split][1][rows],
                    )
                )
            )
    oof = pl.concat(oof_parts).sort("enzyme", "split_id", "arm", "fold", "compound_id")
    oof.with_columns(pl.col(["y_real", "y_visible", "log2fc"]).fill_nan(None)).write_parquet(
        out / "oof_scores.parquet"
    )

    # Evaluation.
    rates = cfg["rates"]
    ranking_rows, curve_rows, nested_rows, boot_rows = [], [], [], []
    mixture_rows = []
    for enzyme in cfg["enzymes"]:
        data = enzyme_data[enzyme]
        arm_names = list(data["arms"])
        for split in cfg["split_strategies"]:
            frames, calls_grid, nested_calls = {}, {}, {}
            for arm in arm_names:
                frame = oof.filter(
                    (pl.col("enzyme") == enzyme)
                    & (pl.col("split_id") == split)
                    & (pl.col("arm") == arm)
                ).sort("fold", "compound_id")
                frame = frame.with_columns(pl.col(["y_real", "y_visible"]).fill_null(np.nan))
                frames[arm] = frame
                score = frame["score"].to_numpy()
                fold = frame["fold"].to_numpy()
                y_real = frame["y_real"].to_numpy()
                stratum = frame["stratum"].to_numpy()
                for name, mask in [
                    ("broad", np.isfinite(y_real)),
                    ("strong", stratum == "strong"),
                    ("weak", stratum == "weak"),
                    ("missing", stratum == "missing"),
                ]:
                    if mask.sum() == 0:
                        continue
                    auroc, ap = ranking(np.where(mask, y_real, np.nan), score)
                    ranking_rows.append(
                        dict(
                            enzyme=enzyme,
                            split_id=split,
                            arm=arm,
                            population=name,
                            n=int(mask.sum()),
                            positives=int(np.nansum(y_real[mask])),
                            auroc=auroc,
                            ap=ap,
                        )
                    )
                # Calls come from every development molecule of the held fold.
                per_rate = {}
                for rate in rates:
                    calls = np.zeros(len(score), dtype=bool)
                    for f in np.unique(fold):
                        rows = np.flatnonzero(fold == f)
                        calls[rows] = top_calls(score[rows], rate)
                    per_rate[rate] = calls
                    curve_rows.append(
                        dict(
                            enzyme=enzyme,
                            split_id=split,
                            arm=arm,
                            rate=rate,
                            n_calls=int(calls.sum()),
                            mcc_real_broad=mcc(y_real, calls),
                            mcc_real_strong=mcc(
                                np.where(stratum == "strong", y_real, np.nan), calls
                            ),
                            mcc_visible=mcc(frame["y_visible"].to_numpy(), calls),
                            calls_labelled=int((calls & np.isfinite(y_real)).sum()),
                        )
                    )
                calls_grid[arm] = per_rate
                # Composition-matched metric: reweight the real-label strata to a target
                # strong-primary fraction, so the judge matches the blind scored mixture.
                strata_masks = dict(
                    strong=(stratum == "strong") & np.isfinite(y_real),
                    weak=(stratum == "weak") & np.isfinite(y_real),
                )
                for phi in cfg["mixture_phi"]:
                    for rate in rates:
                        value, effective = mixture_mcc(
                            y_real,
                            per_rate[rate],
                            strata_masks,
                            dict(strong=phi, weak=1.0 - phi),
                        )
                        mixture_rows.append(
                            dict(
                                enzyme=enzyme,
                                split_id=split,
                                arm=arm,
                                phi_strong=phi,
                                rate=rate,
                                mcc_mixture=value,
                                effective_n=effective,
                                n_strong=int(strata_masks["strong"].sum()),
                                n_weak=int(strata_masks["weak"].sum()),
                            )
                        )
                # Nested rate selection on the other folds, using only this arm's visible labels.
                y_visible = frame["y_visible"].to_numpy()
                selected = np.zeros(len(score), dtype=bool)
                for f in np.unique(fold):
                    other = fold != f
                    inner = {r: mcc(y_visible[other], per_rate[r][other]) for r in rates}
                    rate = max(rates, key=lambda r: (inner[r], -r))
                    held_mask = fold == f
                    selected[held_mask] = per_rate[rate][held_mask]
                    nested_rows.append(
                        dict(
                            enzyme=enzyme,
                            split_id=split,
                            arm=arm,
                            fold=int(f),
                            selected_rate=rate,
                            selection_mcc_visible=inner[rate],
                            mcc_real_broad=mcc(y_real[held_mask], selected[held_mask]),
                        )
                    )
                nested_rows.append(
                    dict(
                        enzyme=enzyme,
                        split_id=split,
                        arm=arm,
                        fold=-1,
                        selected_rate=None,
                        selection_mcc_visible=None,
                        mcc_real_broad=mcc(y_real, selected),
                        n_calls=int(selected.sum()),
                    )
                )
                nested_calls[arm] = selected
            # Paired group bootstrap against the induced-defect arm, on real labels.
            base = frames[DEFECT]
            y_real = base["y_real"].to_numpy()
            groups = base["group"].to_numpy()
            stratum = base["stratum"].to_numpy()
            for arm in arm_names:
                if arm == DEFECT:
                    continue
                assert frames[arm]["compound_id"].equals(base["compound_id"])
                for operating, calls_map in [
                    ("nested", {a: nested_calls[a] for a in arm_names}),
                    *[
                        (f"rate_{r}", {a: calls_grid[a][r] for a in arm_names})
                        for r in cfg["bootstrap_rates"]
                    ],
                ]:
                    populations = [
                        ("broad", np.isfinite(y_real), None),
                        ("strong", stratum == "strong", None),
                    ]
                    strong_mask = (stratum == "strong") & np.isfinite(y_real)
                    weak_mask = (stratum == "weak") & np.isfinite(y_real)
                    for phi in cfg["bootstrap_phi"]:
                        if weak_mask.sum() == 0 or strong_mask.sum() == 0:
                            continue
                        weights = np.zeros(len(y_real))
                        weights[strong_mask] = phi / strong_mask.sum()
                        weights[weak_mask] = (1.0 - phi) / weak_mask.sum()
                        populations.append(
                            (f"mixture_phi_{phi:.2f}", strong_mask | weak_mask, weights)
                        )
                    for population, mask, weights in populations:
                        boot_rows.append(
                            dict(
                                enzyme=enzyme,
                                split_id=split,
                                arm=arm,
                                control=DEFECT,
                                operating_point=operating,
                                population=population,
                                label_set="real",
                                **group_bootstrap(
                                    np.where(mask, y_real, np.nan),
                                    groups,
                                    calls_map[arm],
                                    calls_map[DEFECT],
                                    cfg["bootstrap_samples"],
                                    cfg["seed"],
                                    weights,
                                ),
                            )
                        )
    pl.DataFrame(ranking_rows).write_csv(out / "ranking_metrics.csv")
    pl.DataFrame(mixture_rows).write_csv(out / "mixture_curves.csv")
    pl.DataFrame(curve_rows).write_csv(out / "rate_curves.csv")
    pl.DataFrame(nested_rows, infer_schema_length=None).write_csv(out / "nested_selection.csv")
    pl.DataFrame(boot_rows).write_csv(out / "bootstrap_deltas.csv")

    # Recovery fraction: how much of the defect-to-ceiling gap the fix recovers.
    recovery = []
    pooled_nested = {
        (r["enzyme"], r["split_id"], r["arm"]): r for r in nested_rows if r["fold"] == -1
    }
    for (enzyme, split, arm), row in pooled_nested.items():
        if arm == DEFECT:
            continue
        defect = pooled_nested[(enzyme, split, DEFECT)]["mcc_real_broad"]
        entry = dict(
            enzyme=enzyme,
            split_id=split,
            arm=arm,
            mcc_defect=defect,
            mcc_arm=row["mcc_real_broad"],
            recovered=row["mcc_real_broad"] - defect,
        )
        for name, key in [("screened", CEILING_SCREENED), ("complete", CEILING)]:
            ceiling = pooled_nested[(enzyme, split, key)]["mcc_real_broad"]
            gap = ceiling - defect
            entry[f"mcc_ceiling_{name}"] = ceiling
            entry[f"gap_{name}"] = gap
            entry[f"recovered_fraction_{name}"] = (
                (row["mcc_real_broad"] - defect) / gap if gap > 0 else None
            )
        recovery.append(entry)
    pl.DataFrame(recovery).write_csv(out / "recovery.csv")

    summary = dict(
        config=cfg,
        census=census,
        ranking=ranking_rows,
        nested_pooled=[r for r in nested_rows if r["fold"] == -1],
        nested_per_fold=[r for r in nested_rows if r["fold"] >= 0],
        recovery=recovery,
        mixture=mixture_rows,
        bootstrap=boot_rows,
        runtime_seconds=round(time.time() - started, 1),
        caveats=[
            "The induced defect mimics D6 coverage; it does not reproduce D6 assay chemistry",
            "One LightGBM recipe is used for every arm; the deployed A4 ensemble is not refitted",
            "1A2 and 2C9 carry 23 and 38 positives, so their intervals are wide by construction",
            "fixed_full trains on unlabelled weak rows whose truth is bounded, not observed",
            "A rate is a fraction of all 6,146 development rows, because calls use complete pools",
        ],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    report = [
        "# Induced D6 coverage defect and the pseudo-negative fix",
        "",
        "## Label census",
        "",
        markdown_table(census, list(census[0].keys())),
        "",
        "## Nested-rate MCC on real complete labels (pooled held folds)",
        "",
        markdown_table(
            recovery,
            [
                "enzyme",
                "split_id",
                "arm",
                "mcc_defect",
                "mcc_arm",
                "mcc_ceiling_screened",
                "mcc_ceiling_complete",
                "recovered",
                "recovered_fraction_screened",
            ],
        ),
        "",
        "## Ranking on real labels",
        "",
        markdown_table(
            ranking_rows,
            ["enzyme", "split_id", "arm", "population", "n", "positives", "auroc", "ap"],
        ),
        "",
        "## Paired group bootstrap versus the induced-defect arm (real labels)",
        "",
        markdown_table(
            boot_rows,
            [
                "enzyme",
                "split_id",
                "arm",
                "operating_point",
                "population",
                "delta",
                "ci_low",
                "ci_high",
                "p_delta_gt_0",
            ],
        ),
        "",
    ]
    (out / "report.md").write_text("\n".join(report))
    print(json.dumps(dict(runtime_seconds=summary["runtime_seconds"])), flush=True)


if __name__ == "__main__":
    main()
