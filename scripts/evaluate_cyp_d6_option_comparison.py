"""Compare three D6 options in public-metric units on the recorded sizing basis.

The options are the gate-filtered vector at an unchanged 49 calls, the pure monotone
extension of the deployed vector to the top 94 rows that run
20261003T025354186349Z-621b6ad3 pre-registered, and holding the current vector.

Every option is evaluated against the fitted confusion matrix of the deployed payload from
run 20261003T012000305555Z-bca66716: 374 scored rows of 750, 26 scored positives, 28 scored
calls with 11 true positives and 17 false positives, point MCC 0.36165. A swap of k calls
at a fixed count is simulated by drawing the true positives among the scored removed calls
from the deployed realised precision and those among the scored added calls from each
marginal-precision scenario, then recomputing MCC exactly rather than linearising it. The
reported interval is the simulated 10th to 90th percentile, so it is a spread over outcomes
at a fixed assumed precision and not a posterior over the precision itself.

Read-only arithmetic. Nothing is fitted and no candidate, prepared submission or upload is
created.
"""

import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl


def mcc_from_counts(tp, fp, fn, tn):
    tp, fp, fn, tn = (np.asarray(v, dtype=float) for v in (tp, fp, fn, tn))
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return np.divide(tp * tn - fp * fn, denom, out=np.zeros_like(denom), where=denom > 0)


def simulate_swap(sizing, scored_swapped, removed_precision, added_precision, rng, draws):
    """Swap `scored_swapped` scored calls at a fixed total call count."""
    tp0 = sizing["deployed_true_positives"]
    fp0 = sizing["deployed_false_positives"]
    fn0 = sizing["deployed_false_negatives"]
    scored = sizing["scored_rows"]
    calls = sizing["deployed_scored_calls"]
    tn0 = scored - calls - fn0
    k = int(round(scored_swapped))
    k_removed = min(k, tp0 + fp0)
    tp_out = rng.binomial(k_removed, removed_precision, draws)
    tp_out = np.minimum(tp_out, tp0)
    tp_in = rng.binomial(k, added_precision, draws)
    tp = tp0 - tp_out + tp_in
    fp = fp0 + tp_out - tp_in
    fn = fn0 + tp_out - tp_in
    tn = tn0 - tp_out + tp_in
    return mcc_from_counts(tp, fp, fn, tn)


def simulate_extension(sizing, scored_added, added_precision, rng, draws):
    """Pure monotone addition: nothing is removed, so only the call count grows."""
    tp0 = sizing["deployed_true_positives"]
    fp0 = sizing["deployed_false_positives"]
    fn0 = sizing["deployed_false_negatives"]
    scored = sizing["scored_rows"]
    calls = sizing["deployed_scored_calls"]
    tn0 = scored - calls - fn0
    k = int(round(scored_added))
    tp_in = rng.binomial(k, added_precision, draws)
    tp_in = np.minimum(tp_in, fn0)
    tp = tp0 + tp_in
    fp = fp0 + k - tp_in
    fn = fn0 - tp_in
    tn = tn0 - (k - tp_in)
    return mcc_from_counts(tp, fp, fn, tn)


def quantiles(values, interval):
    low = (1 - interval) / 2
    return float(np.quantile(values, low)), float(np.quantile(values, 1 - low))


def fmt(value, digits=4):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "–"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:+.{digits}f}" if abs(value) < 1 else f"{value:.{digits}f}"
    return str(value)


def markdown_table(rows, columns, digits=4):
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(c), digits) for c in columns) + " |")
    return "\n".join(lines)


def main():
    started = time.time()
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    sizing = cfg["sizing"]
    base = float(
        mcc_from_counts(
            sizing["deployed_true_positives"],
            sizing["deployed_false_positives"],
            sizing["deployed_false_negatives"],
            sizing["scored_rows"]
            - sizing["deployed_scored_calls"]
            - sizing["deployed_false_negatives"],
        )
    )
    assert abs(base - sizing["deployed_point_mcc"]) < 5e-5, base
    rng = np.random.default_rng(cfg["seed"])
    draws = cfg["simulations"]
    interval = cfg["interval"]

    swaps = pl.read_csv(source / cfg["inputs"]["blind_swaps"]).filter(
        (pl.col("reference") == "deployed_oct03") & (pl.col("n_calls") == 49)
    )
    bands = pl.read_csv(source / cfg["inputs"]["marginal_precision"])

    rows = []
    # Option 1: the gate filter at an unchanged call count.
    for record in swaps.to_dicts():
        changed = int(record["removed"])
        if changed == 0:
            continue
        neighbour_ratio = (
            record["mean_neighbour_rate_added"] / record["mean_neighbour_rate_removed"]
            if record["mean_neighbour_rate_removed"]
            else None
        )
        scenarios = dict(cfg["added_precision_scenarios"])
        if neighbour_ratio is not None:
            scenarios["neighbour_rate_ratio"] = cfg["removed_precision"] * float(neighbour_ratio)
        for name, added_precision in scenarios.items():
            scored_swapped = changed * cfg["scored_fraction_swap"]
            values = simulate_swap(
                sizing,
                scored_swapped,
                cfg["removed_precision"],
                added_precision,
                rng,
                draws,
            )
            low, high = quantiles(values - base, interval)
            rows.append(
                dict(
                    option="gate_filter",
                    variant=record["model"],
                    calls=49,
                    changed_calls=changed,
                    scored_changed=round(scored_swapped, 1),
                    scenario=name,
                    added_precision=float(added_precision),
                    neighbour_rate_ratio=neighbour_ratio,
                    expected_delta=float((values - base).mean()),
                    sd_delta=float((values - base).std()),
                    interval_low=low,
                    interval_high=high,
                    p_delta_positive=float((values > base).mean()),
                )
            )
    # Option 2: the pre-registered monotone extension to 94 calls.
    for name, added_precision in cfg["extension_precision_scenarios"].items():
        scored_added = cfg["extension_added_calls"] * cfg["scored_fraction_extension"]
        values = simulate_extension(sizing, scored_added, added_precision, rng, draws)
        low, high = quantiles(values - base, interval)
        rows.append(
            dict(
                option="extension_94",
                variant="l4_top94",
                calls=94,
                changed_calls=cfg["extension_added_calls"],
                scored_changed=round(scored_added, 1),
                scenario=name,
                added_precision=float(added_precision),
                expected_delta=float((values - base).mean()),
                sd_delta=float((values - base).std()),
                interval_low=low,
                interval_high=high,
                p_delta_positive=float((values > base).mean()),
            )
        )
    # Option 3: hold.
    rows.append(
        dict(
            option="hold",
            variant="deployed_oct03",
            calls=49,
            changed_calls=0,
            scored_changed=0.0,
            scenario="unchanged",
            added_precision=None,
            expected_delta=0.0,
            sd_delta=0.0,
            interval_low=0.0,
            interval_high=0.0,
            p_delta_positive=0.0,
        )
    )
    options = pl.DataFrame(rows, infer_schema_length=None)
    options.write_csv(out / "option_comparison.csv")

    # Sensitivity of every option to the scored fraction of the changed ranks.
    sensitivity = []
    for fraction in cfg["scored_fraction_sensitivity"]:
        for record in swaps.filter(pl.col("model") == "l4_filter_tau0.2").to_dicts():
            values = simulate_swap(
                sizing,
                record["removed"] * fraction,
                cfg["removed_precision"],
                cfg["added_precision_scenarios"]["offline_mixture_ratio"],
                rng,
                draws,
            )
            low, high = quantiles(values - base, interval)
            sensitivity.append(
                dict(
                    option="gate_filter",
                    variant=record["model"],
                    scored_fraction=fraction,
                    expected_delta=float((values - base).mean()),
                    interval_low=low,
                    interval_high=high,
                )
            )
        values = simulate_extension(
            sizing,
            cfg["extension_added_calls"] * fraction,
            cfg["extension_precision_scenarios"]["uncalibrated_mixture"],
            rng,
            draws,
        )
        low, high = quantiles(values - base, interval)
        sensitivity.append(
            dict(
                option="extension_94",
                variant="l4_top94",
                scored_fraction=fraction,
                expected_delta=float((values - base).mean()),
                interval_low=low,
                interval_high=high,
            )
        )
    pl.DataFrame(sensitivity).write_csv(out / "scored_fraction_sensitivity.csv")

    # How often the pre-registered extension produces a decisive outcome.
    decisive = []
    for name, added_precision in cfg["extension_precision_scenarios"].items():
        values = simulate_extension(
            sizing,
            cfg["extension_added_calls"] * cfg["scored_fraction_extension"],
            added_precision,
            rng,
            draws,
        )
        delta = values - base
        decisive.append(
            dict(
                scenario=name,
                added_precision=float(added_precision),
                p_above_plus_080=float((delta > 0.080).mean()),
                p_below_minus_080=float((delta < -0.080).mean()),
                p_decisive=float(((delta > 0.080) | (delta < -0.080)).mean()),
                p_inside=float(((delta <= 0.080) & (delta >= -0.080)).mean()),
            )
        )
    pl.DataFrame(decisive).write_csv(out / "extension_decisiveness.csv")

    summary = dict(
        config=cfg,
        deployed_point_mcc=base,
        band_precision=bands.to_dicts(),
        options=rows,
        scored_fraction_sensitivity=sensitivity,
        extension_decisiveness=decisive,
        runtime_seconds=round(time.time() - started, 1),
        caveats=[
            "The confusion matrices are reconstructions from published bootstrap aggregates, "
            "not disclosed labels",
            "Each interval is a spread over outcomes at one assumed marginal precision; it is "
            "not a posterior over that precision",
            "The scored fraction of the changed ranks is assumed, and the fitted 49-call "
            "payloads are scored at 0.57 to 0.65 while the overall fraction is 0.499",
            "The neighbour-rate scenario transfers a labelled-neighbour ratio to a Live "
            "precision ratio, which is an assumption of proportionality",
        ],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    report = [
        "# Three D6 options in public-metric units",
        "",
        f"Deployed point MCC {base:.5f} against a reported Live {sizing['deployed_live_mcc']:.6f}.",
        "",
        "## Expected change and spread",
        "",
        markdown_table(
            rows,
            [
                "option",
                "variant",
                "calls",
                "changed_calls",
                "scored_changed",
                "scenario",
                "added_precision",
                "expected_delta",
                "sd_delta",
                "interval_low",
                "interval_high",
                "p_delta_positive",
            ],
        ),
        "",
        "## Sensitivity to the scored fraction of the changed ranks",
        "",
        markdown_table(
            sensitivity,
            [
                "option",
                "variant",
                "scored_fraction",
                "expected_delta",
                "interval_low",
                "interval_high",
            ],
        ),
        "",
        "## Probability that the pre-registered extension gives a decisive outcome",
        "",
        markdown_table(
            decisive,
            [
                "scenario",
                "added_precision",
                "p_above_plus_080",
                "p_below_minus_080",
                "p_decisive",
                "p_inside",
            ],
        ),
        "",
    ]
    (out / "report.md").write_text("\n".join(report) + "\n")
    print(json.dumps(dict(runtime_seconds=summary["runtime_seconds"]), indent=2), flush=True)


if __name__ == "__main__":
    main()
