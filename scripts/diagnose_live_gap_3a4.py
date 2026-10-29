"""Read-only calibration of the 3A4 live-board result against development OOF.

Reconstructs the deployed fixed-weight composition on the 2,300-row incumbent cohort for the sent
(baseline) and TDI-proxy variants, then asks how large and how variable the proxy-minus-baseline
ST-RAE difference is on random series-grouped subsets of roughly the live-board size (about 375
compounds), and how large the per-compound prediction change is compared with the blind change.
"""

import importlib.util
import json
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "artifacts/runs/regression/20261029T054536454025Z-220f26c8"
ENDPOINT = "CYP3A4_pIC50_direct_inhibition"


def main():
    inputs, cfg = RUN / "inputs", json.loads((RUN / "config.json").read_text())
    spec = importlib.util.spec_from_file_location("official", inputs / "custom_scoring_functions.py")
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    weights = cfg["weights"]
    cohort = (
        pl.read_parquet(inputs / "control-3a4.parquet")
        .filter((pl.col("method") == "caruana") & (pl.col("endpoint") == ENDPOINT))
        .sort("compound_id")
        .select("compound_id", "fold", "y_true", "conf_low", "conf_high")
    )

    def column(file, method=None):
        frame = pl.read_parquet(inputs / file).filter(pl.col("endpoint") == ENDPOINT)
        if method is not None:
            frame = frame.filter(pl.col("method") == method)
        return cohort.select("compound_id").join(
            frame.select("compound_id", "prediction"), on="compound_id", how="left"
        )["prediction"].to_numpy()

    members = {name: column(file) for name, file in cfg["member_oof"].items()}
    arms = {}
    for arm in ("baseline", "proxy"):
        swapped = {
            name: column(s["file"], f"{s['method']}--{arm}")
            for name, s in cfg["swapped_members"].items()
        }
        arms[arm] = sum(weights[n] * swapped.get(n, members[n]) for n in weights)
    y, lo, hi = (cohort[c].to_numpy() for c in ("y_true", "conf_low", "conf_high"))

    def st(rows, p):
        return float(
            scorer.rae_soft_threshold_absolute_error(
                y[rows], p[rows], y_true_upper=hi[rows], y_true_lower=lo[rows]
            )
        )

    groups = (
        cohort.join(
            pl.read_parquet(inputs / "parent-assignments.parquet").select("compound_id", "group"),
            on="compound_id",
            how="left",
        )["group"].to_numpy()
    )
    unique = np.unique(groups)
    index = {g: np.flatnonzero(groups == g) for g in unique}
    rng = np.random.default_rng(1)
    delta_pred = arms["proxy"] - arms["baseline"]
    print("dev per-compound |proxy - baseline| mean", float(np.abs(delta_pred).mean()))
    print("dev corr(proxy-baseline, y - baseline)", float(np.corrcoef(delta_pred, y - arms["baseline"])[0, 1]))
    full = st(np.arange(len(y)), arms["proxy"]) - st(np.arange(len(y)), arms["baseline"])
    print("dev full delta", full)
    for size in (375, 750):
        deltas = []
        for _ in range(4000):
            order = rng.permutation(unique)
            picked, count = [], 0
            for g in order:
                picked.append(index[g])
                count += len(index[g])
                if count >= size:
                    break
            rows = np.concatenate(picked)
            deltas.append(st(rows, arms["proxy"]) - st(rows, arms["baseline"]))
        d = np.array(deltas)
        print(
            f"subset ~{size}: sd {d.std():.4f}, mean {d.mean():.4f}, "
            f"P(delta >= +0.0065) {np.mean(d >= 0.0065):.4f}, 5-95% {np.quantile(d, [.05, .95]).round(4)}"
        )


if __name__ == "__main__":
    main()
