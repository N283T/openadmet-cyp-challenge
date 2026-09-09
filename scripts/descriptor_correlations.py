"""Descriptive label correlations on measured training compounds; no model or split."""

import json
import os
import warnings
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import false_discovery_control, pearsonr, spearmanr

from cyp_challenge.storage import write_json


def main():
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    inputs = Path(os.environ["CYP_INPUT_DIR"])
    output = Path(os.environ["CYP_OUTPUT_DIR"])
    measurements = pl.read_parquet(next(inputs.glob("*measurements.parquet")))
    families = []
    for family in ["rdkit", "mordred"]:
        frame = pl.concat(
            [pl.read_parquet(p) for p in sorted(inputs.glob(f"*-{family}-*.parquet"))]
        )
        assert frame["compound_id"].n_unique() == frame.height
        families.append(
            frame.rename({c: f"{family}:{c}" for c in frame.columns if c != "compound_id"})
        )
    features = families[0].join(families[1], on="compound_id", how="inner", validate="1:1")
    names = [c for c in features.columns if c != "compound_id"]
    all_results, endpoint_info, tops = [], [], []
    rng = np.random.default_rng(config["seed"])
    for endpoint in config["endpoints"]:
        target = measurements.filter(pl.col("endpoint") == endpoint).select("compound_id", "value")
        assert target["compound_id"].n_unique() == target.height, (
            "Resolve conflicting observations first"
        )
        joined = target.join(features, on="compound_id", how="left", validate="1:1").sort(
            "compound_id"
        )
        assert joined.height == target.height
        y = joined["value"].to_numpy()
        matrix = joined.select(names).to_numpy()
        info = {
            "endpoint": endpoint,
            "n": len(y),
            "min": float(y.min()),
            "max": float(y.max()),
            "positive_n": int(np.sum(y == 1)) if endpoint.endswith("is_TDI") else None,
            "below_4_n": int(np.sum(y < 4)) if not endpoint.endswith("is_TDI") else None,
        }
        endpoint_info.append(info)
        candidates = []
        for i, name in enumerate(names):
            x = matrix[:, i]
            mask = np.isfinite(x) & np.isfinite(y)
            n = int(mask.sum())
            record = {
                "endpoint": endpoint,
                "descriptor": name,
                "n": n,
                "coverage": n / len(y),
                "pearson": None,
                "spearman": None,
                "p": None,
                "status": "insufficient",
            }
            if n >= config["minimum_n"] and np.unique(x[mask]).size > 1:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    r = pearsonr(x[mask], y[mask])
                    rho = spearmanr(x[mask], y[mask])
                if np.isfinite(rho.statistic):
                    record.update(
                        pearson=float(r.statistic),
                        spearman=float(rho.statistic),
                        p=float(rho.pvalue),
                        status="ok",
                    )
                    if record["coverage"] >= config["minimum_coverage"]:
                        candidates.append((abs(record["spearman"]), i, record))
            all_results.append(record)
        candidates.sort(key=lambda x: x[0], reverse=True)
        top = dict(candidates[0][2])
        x = matrix[:, candidates[0][1]]
        mask = np.isfinite(x) & np.isfinite(y)
        x, yy = x[mask], y[mask]
        boot = []
        for _ in range(config["bootstrap_repetitions"]):
            idx = rng.integers(0, len(yy), len(yy))
            boot.append(float(spearmanr(x[idx], yy[idx]).statistic))
        top["bootstrap_95_low"], top["bootstrap_95_high"] = np.nanquantile(
            boot, [0.025, 0.975]
        ).tolist()
        if not endpoint.endswith("is_TDI"):
            reliable = yy >= 4
            top["sensitivity_ge4_n"] = int(reliable.sum())
            top["sensitivity_ge4_rho"] = float(spearmanr(x[reliable], yy[reliable]).statistic)
        # Independently recompute Pearson using NumPy covariance on the same paired rows.
        assert np.isclose(np.corrcoef(x, yy)[0, 1], top["pearson"], atol=1e-12)
        tops.append(top)
        print(endpoint, info["n"], top["descriptor"], top["spearman"], flush=True)
    valid = [r for r in all_results if r["status"] == "ok"]
    for row, q in zip(valid, false_discovery_control([r["p"] for r in valid]), strict=True):
        row["q_bh"] = float(q)
    table = pl.DataFrame(all_results, infer_schema_length=None)
    table.write_csv(output / "correlations.csv")
    pl.DataFrame(endpoint_info).write_csv(output / "endpoints.csv")
    write_json(output / "top_correlations.json", tops)
    # Simple, interpretable descriptors on every endpoint for comparison with the maxima.
    common = [
        "MolWt",
        "MolLogP",
        "TPSA",
        "NumHDonors",
        "NumHAcceptors",
        "NumRotatableBonds",
        "RingCount",
        "FractionCSP3",
    ]
    table.filter(pl.col("descriptor").is_in(["rdkit:" + c for c in common])).write_csv(
        output / "common_descriptors.csv"
    )
    lines = [
        "# Descriptor–target correlation screening",
        "",
        "Source dataset: `"
        + config["dataset_run"]
        + "`; feature set: `"
        + config["feature_run"]
        + "`.",
        "Measured training labels only, one row per standardized compound and endpoint.",
        "Pairwise complete cases, no imputation, models, splits, or feature removal.",
        "Rankings require n >= 100 and >= 80% coverage. Full results: correlations.csv.",
        "BH q-values span all pairs; molecular dependence limits inference.",
        "This is exploratory selection on the same data, not out-of-sample performance.",
        "Bootstrap intervals condition on feature selection and assume independent molecules.",
        "",
    ]
    for endpoint in config["endpoints"]:
        lines += [
            "## " + endpoint,
            "",
            "| Descriptor | n | Coverage | Pearson | Spearman |",
            "|---|---:|---:|---:|---:|",
        ]
        selected = sorted(
            [
                r
                for r in valid
                if r["endpoint"] == endpoint and r["coverage"] >= config["minimum_coverage"]
            ],
            key=lambda r: abs(r["spearman"]),
            reverse=True,
        )[:10]
        for row in selected:
            lines.append(
                f"| {row['descriptor']} | {row['n']} | {row['coverage']:.1%} | "
                f"{row['pearson']:.3f} | {row['spearman']:.3f} |"
            )
        lines.append("")
    (output / "report.md").write_text("\n".join(lines))
    write_json(
        output / "summary.json",
        {
            "config": config,
            "endpoints": endpoint_info,
            "total_pairs": len(all_results),
            "valid_pairs": len(valid),
            "qualification": "Exploratory associations; no causal or predictive validation claims",
            "qa": "Unique keys and join size checked; top Pearson recomputed with NumPy",
        },
    )


if __name__ == "__main__":
    main()
