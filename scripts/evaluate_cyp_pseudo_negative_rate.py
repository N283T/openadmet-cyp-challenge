"""Direct, non-circular test of the pseudo-negative assumption on three enzymes.

Every supplied CYP2D6_is_TDI label belongs to a primary-screen hit, so the D6
"pseudo-negative" device (unlabelled rows with primary log2fc >= -0.25 relabelled
False) cannot be checked against D6 labels. CYP1A2, CYP2C9 and CYP3A4 do carry
labels in the weak-primary region, so the assumption is directly measurable there.

Part A measures the positive rate among labelled weak-primary rows, pooled and per
enzyme, with exact Clopper-Pearson upper bounds, and sweeps the log2fc cut so that
-0.25 is a point on a curve rather than a premise. It also reports the
missing-primary stratum, which the D6 device never relabels.

Part B asks how many truly positive rows would have to hide inside the 2,693 D6
pseudo-negatives to overturn the pool-versus-official ranking, using the frozen
held-fold scores of run 20261002T215345688167Z-7920b49a. The required count is
converted to a required positive rate and compared with the Part A bound.

Nothing is fitted. No candidate, prepared submission or upload results.
"""

import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl
from cyp_primary_supervision import load_primary
from scipy.stats import beta

ENZYMES = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
VALIDATING = ["CYP1A2", "CYP2C9", "CYP3A4"]
CONTROL = "official_lgbm"
POOL = "pool_lgbm"
HISTORICAL_RATE = 49 / 750


def clopper_pearson(successes, trials, alpha):
    """One-sided exact upper bound and the two-sided interval for a binomial rate."""
    if trials == 0:
        return dict(rate=None, upper_one_sided=None, ci_low=None, ci_high=None)
    low = 0.0 if successes == 0 else float(beta.ppf(alpha / 2, successes, trials - successes + 1))
    high = (
        1.0
        if successes == trials
        else float(beta.ppf(1 - alpha / 2, successes + 1, trials - successes))
    )
    upper = (
        1.0
        if successes == trials
        else float(beta.ppf(1 - alpha, successes + 1, trials - successes))
    )
    return dict(rate=successes / trials, upper_one_sided=upper, ci_low=low, ci_high=high)


def mcc_from_counts(tp, fp, fn, tn):
    tp, fp, fn, tn = (np.asarray(v, dtype=float) for v in (tp, fp, fn, tn))
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return np.divide(tp * tn - fp * fn, denom, out=np.zeros_like(denom), where=denom > 0)


def top_calls(score, fold, rate):
    calls = np.zeros(len(score), dtype=bool)
    for f in np.unique(fold):
        rows = np.flatnonzero(fold == f)
        k = int(round(rate * len(rows)))
        calls[rows[np.argsort(-score[rows], kind="stable")[:k]]] = True
    return calls


def cell_counts(y, calls):
    """Confusion counts on finite labels only."""
    keep = np.isfinite(y)
    y, calls = y[keep].astype(bool), calls[keep]
    return (
        float((y & calls).sum()),
        float((~y & calls).sum()),
        float((y & ~calls).sum()),
        float((~y & ~calls).sum()),
    )


def flip_delta(base, cells, allocation):
    """MCC delta after relabelling `allocation` pseudo-negatives per call cell as positive.

    `cells` orders the pseudo-negative population as (both, control only, pool only,
    neither). A relabelled row that a model calls moves FP -> TP for that model; a
    relabelled row it does not call moves TN -> FN.
    """
    out = []
    for name, (tp, fp, fn, tn) in base.items():
        called = allocation[0] + (allocation[1] if name == CONTROL else allocation[2])
        uncalled = allocation[3] + (allocation[2] if name == CONTROL else allocation[1])
        out.append(mcc_from_counts(tp + called, fp - called, fn + uncalled, tn - uncalled))
    return float(out[1] - out[0]), cells


def adversarial_minimum(base, cells, k, step):
    """Smallest MCC(pool) - MCC(control) over every way of placing k true positives."""
    best, best_allocation = np.inf, None
    limits = [min(int(c), k) for c in cells]
    grid = [range(0, limit + 1, step) for limit in limits[:3]]
    for a in grid[0]:
        for b in grid[1]:
            if a + b > k:
                break
            for c in grid[2]:
                if a + b + c > k:
                    break
                d = k - a - b - c
                if d > limits[3]:
                    continue
                delta, _ = flip_delta(base, cells, (a, b, c, d))
                if delta < best:
                    best, best_allocation = delta, (a, b, c, d)
    return best, best_allocation


def proportional_delta(base, cells, k):
    """MCC delta when k true positives fall across the call cells in proportion to size."""
    total = sum(cells)
    allocation = tuple(k * c / total for c in cells)
    return flip_delta(base, cells, allocation)[0]


def smallest_flip(base, cells, kind, step, limit):
    """Smallest k whose best (adversarial) or expected (proportional) delta is <= 0."""
    for k in range(0, limit + 1):
        delta = (
            adversarial_minimum(base, cells, k, step)[0]
            if kind == "adversarial"
            else proportional_delta(base, cells, k)
        )
        if delta <= 0:
            return k, delta
    return None, None


def main():
    started = time.time()
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    inputs = cfg["inputs"]
    alpha = cfg["alpha"]
    compounds = pl.read_parquet(source / inputs["compounds"]).sort("compound_id")
    aliases = pl.read_parquet(source / inputs["aliases"])
    test_mapping = pl.read_parquet(source / inputs["test_mapping"])
    measured = pl.read_parquet(source / inputs["measurements"])
    primary, primary_rows = load_primary(source / inputs["primary_screen"], aliases, compounds)
    order = ["CYP2D6", "CYP3A4", "CYP1A2", "CYP2C9"]
    log2fc = {enzyme: primary[:, order.index(enzyme)] for enzyme in ENZYMES}
    blind = set(test_mapping["compound_id"].to_list())
    is_blind = compounds["compound_id"].is_in(list(blind)).to_numpy()
    assert is_blind.sum() == 750
    for enzyme in ENZYMES:
        assert not np.isfinite(log2fc[enzyme][is_blind]).any()

    # Raw Emax labels, used only to reconcile the canonical table's count difference.
    mapping = aliases.select("molecule_name", "compound_id").unique()
    emax = (
        pl.read_csv(source / inputs["emax"])
        .rename({"Molecule_Name": "molecule_name"})
        .join(mapping, on="molecule_name", how="left", validate="m:1")
    )
    canonical = {
        enzyme: compounds.select("compound_id")
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
        for enzyme in ENZYMES
    }

    # Label-source reconciliation.
    provenance = []
    for enzyme in ENZYMES:
        column = f"{enzyme}_is_TDI"
        raw = emax.select("compound_id", pl.col(column).alias("y")).filter(
            pl.col("y").is_not_null()
        )
        raw_screened = raw.join(
            primary_rows.filter(pl.col("enzyme") == enzyme).select("compound_id", "log2fc"),
            on="compound_id",
            how="left",
        )
        y = canonical[enzyme]
        provenance.append(
            dict(
                enzyme=enzyme,
                raw_emax_labelled=raw.height,
                raw_emax_positive=int(raw["y"].sum()),
                raw_emax_unlinked=int(raw["compound_id"].null_count()),
                raw_emax_without_primary=int(raw_screened["log2fc"].null_count()),
                canonical_labelled=int(np.isfinite(y).sum()),
                canonical_positive=int(np.nansum(y)),
                canonical_without_primary=int(
                    (np.isfinite(y) & ~np.isfinite(log2fc[enzyme])).sum()
                ),
                labelled_max_log2fc=float(
                    np.nanmax(np.where(np.isfinite(y), log2fc[enzyme], -np.inf))
                ),
                screened_molecules=int(np.isfinite(log2fc[enzyme]).sum()),
            )
        )
    print(json.dumps(provenance), flush=True)

    # Part A: the positive rate among labelled weak-primary rows.
    cut = cfg["pseudo_negative_min_log2fc"]
    band_low = cfg["borderline_min_log2fc"]
    strata_rows = []
    for enzyme in ENZYMES:
        y, primary_values = canonical[enzyme], log2fc[enzyme]
        labelled = np.isfinite(y)
        screened = np.isfinite(primary_values)
        for name, mask in [
            ("weak_primary", labelled & screened & (primary_values >= cut)),
            (
                "borderline",
                labelled & screened & (primary_values >= band_low) & (primary_values < cut),
            ),
            ("primary_hit", labelled & screened & (primary_values < band_low)),
            ("no_primary", labelled & ~screened),
        ]:
            successes = int(np.nansum(y[mask]))
            strata_rows.append(
                dict(
                    enzyme=enzyme,
                    stratum=name,
                    n=int(mask.sum()),
                    positives=successes,
                    **clopper_pearson(successes, int(mask.sum()), alpha),
                )
            )
    for name in ["weak_primary", "borderline", "primary_hit", "no_primary"]:
        subset = [r for r in strata_rows if r["stratum"] == name and r["enzyme"] in VALIDATING]
        n = sum(r["n"] for r in subset)
        successes = sum(r["positives"] for r in subset)
        strata_rows.append(
            dict(
                enzyme="pooled_1A2_2C9_3A4",
                stratum=name,
                n=n,
                positives=successes,
                **clopper_pearson(successes, n, alpha),
            )
        )
    strata = pl.DataFrame(strata_rows)
    strata.write_csv(out / "label_strata.csv")

    # Part A continued: the cut curve.
    curve_rows = []
    grid = np.round(np.arange(cfg["cut_from"], cfg["cut_to"] - 1e-9, -cfg["cut_step"]), 4)
    for value in grid:
        pooled_n, pooled_pos = 0, 0
        for enzyme in ENZYMES:
            y, primary_values = canonical[enzyme], log2fc[enzyme]
            mask = np.isfinite(y) & np.isfinite(primary_values) & (primary_values >= value)
            successes = int(np.nansum(y[mask]))
            curve_rows.append(
                dict(
                    scope=enzyme,
                    cut=float(value),
                    n=int(mask.sum()),
                    positives=successes,
                    **clopper_pearson(successes, int(mask.sum()), alpha),
                )
            )
            if enzyme in VALIDATING:
                pooled_n += int(mask.sum())
                pooled_pos += successes
        curve_rows.append(
            dict(
                scope="pooled_1A2_2C9_3A4",
                cut=float(value),
                n=pooled_n,
                positives=pooled_pos,
                **clopper_pearson(pooled_pos, pooled_n, alpha),
            )
        )
    curve = pl.DataFrame(curve_rows)
    curve.write_csv(out / "cut_curve.csv")

    # Part B: how many hidden D6 positives would overturn the pool-versus-official ranking.
    oof = pl.read_parquet(source / inputs["testlike_oof"])
    assert set(oof["pool_role"].unique().to_list()) >= {"official", "pseudo_negative"}
    overturn_rows = []
    for split in oof["split_id"].unique().to_list():
        frames = {}
        for model in (CONTROL, POOL):
            frames[model] = (
                oof.filter((pl.col("split_id") == split) & (pl.col("model") == model))
                .sort("compound_id")
                .with_columns(pl.col(["y_official", "y_testlike"]).fill_null(np.nan))
            )
        assert frames[CONTROL]["compound_id"].equals(frames[POOL]["compound_id"])
        role = frames[CONTROL]["pool_role"].to_numpy()
        y_testlike = frames[CONTROL]["y_testlike"].to_numpy()
        pseudo = role == "pseudo_negative"
        scored = np.isfinite(y_testlike)
        for rate in cfg["overturn_rates"]:
            calls = {
                model: top_calls(frame["score"].to_numpy(), frame["fold"].to_numpy(), rate)
                for model, frame in frames.items()
            }
            base = {model: cell_counts(y_testlike, calls[model]) for model in (CONTROL, POOL)}
            point = float(mcc_from_counts(*base[POOL]) - mcc_from_counts(*base[CONTROL]))
            cells = (
                float((pseudo & calls[CONTROL] & calls[POOL]).sum()),
                float((pseudo & calls[CONTROL] & ~calls[POOL]).sum()),
                float((pseudo & ~calls[CONTROL] & calls[POOL]).sum()),
                float((pseudo & ~calls[CONTROL] & ~calls[POOL]).sum()),
            )
            row = dict(
                split_id=split,
                rate=rate,
                scored_pool=int(scored.sum()),
                pseudo_negatives=int(pseudo.sum()),
                delta_mcc=point,
                calls_control=int(calls[CONTROL].sum()),
                calls_pool=int(calls[POOL].sum()),
                pseudo_both=int(cells[0]),
                pseudo_control_only=int(cells[1]),
                pseudo_pool_only=int(cells[2]),
                pseudo_neither=int(cells[3]),
            )
            for kind in ("adversarial", "proportional"):
                k, delta = smallest_flip(
                    base, cells, kind, cfg["overturn_step"], cfg["overturn_limit"]
                )
                row[f"{kind}_k"] = k
                row[f"{kind}_delta_at_k"] = delta
                row[f"{kind}_required_rate"] = None if k is None else k / int(pseudo.sum())
            overturn_rows.append(row)
    overturn = pl.DataFrame(overturn_rows, infer_schema_length=None)
    overturn.write_csv(out / "overturn_counts.csv")

    pooled_weak = next(
        r
        for r in strata_rows
        if r["enzyme"] == "pooled_1A2_2C9_3A4" and r["stratum"] == "weak_primary"
    )
    summary = dict(
        config=cfg,
        provenance=provenance,
        primary_source_rows=int(primary_rows["source_rows"].sum()),
        pooled_weak_primary=pooled_weak,
        strata=strata_rows,
        cut_curve=curve_rows,
        overturn=overturn_rows,
        runtime_seconds=round(time.time() - started, 1),
        caveats=[
            "The weak-primary rate is measured on 1A2/2C9/3A4 and transferred to D6 by assumption",
            "The D6 assay differs in readout, so a shared mechanism is argued, not demonstrated",
            "Part B conditions on the frozen held-fold scores of one run and two recipes",
            "The adversarial minimum is a grid search over allocations, not a closed-form proof",
        ],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(dict(pooled_weak_primary=pooled_weak, seconds=summary["runtime_seconds"])))


if __name__ == "__main__":
    main()
