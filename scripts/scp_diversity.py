"""Describe fixed D6 predictions and residual diversity without fitting weights."""

import itertools
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr, spearmanr


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    frame = pl.read_parquet(inp / cfg["oof"]).filter(
        pl.col("endpoint") == "CYP2D6_pIC50_direct_inhibition"
    )
    keys = ["compound_id", "endpoint", "fold"]
    arrays, errors, metrics = {}, {}, []
    reference = frame.filter(pl.col("method") == "previous_global").sort(keys)
    assert reference.height == 1457
    for (method,), table in frame.group_by("method"):
        table = table.sort(keys)
        assert table.height == 1457 and table.unique(keys).height == 1457
        assert table.select(keys).equals(reference.select(keys))
        np.testing.assert_array_equal(table["y_true"], reference["y_true"])
        y, pred = table["y_true"].to_numpy(), table["prediction"].to_numpy()
        assert np.isfinite(pred).all()
        residual = pred - y
        arrays[method], errors[method] = pred, residual
        rho = float(spearmanr(y, pred).statistic)
        metrics.append(
            dict(
                method=method,
                n=len(y),
                prediction_mean=float(pred.mean()),
                prediction_sd=float(pred.std()),
                label_mean=float(y.mean()),
                bias=float(residual.mean()),
                rmse=float(np.sqrt(np.mean(residual**2))),
                mae=float(np.abs(residual).mean()),
                residual_sd=float(residual.std()),
                spearman=rho,
                pearson=float(pearsonr(y, pred).statistic),
                spearman_minus_blog_reported_0503=rho - 0.503,
            )
        )
    pairs = []
    for a, b in itertools.combinations(sorted(arrays), 2):
        ea, eb = errors[a], errors[b]
        mean_error = (ea + eb) / 2
        component_mse = (np.mean(ea**2) + np.mean(eb**2)) / 2
        pair_mse = np.mean(mean_error**2)
        disagreement_gain = np.mean((ea - eb) ** 2) / 4
        assert np.isclose(component_mse - pair_mse, disagreement_gain)
        pairs.append(
            dict(
                first=a,
                second=b,
                n=1457,
                prediction_pearson=float(pearsonr(arrays[a], arrays[b]).statistic),
                residual_pearson=float(pearsonr(ea, eb).statistic),
                error_cross_moment=float(np.mean(ea * eb)),
                mean_component_mse=float(component_mse),
                equal_pair_mse=float(pair_mse),
                disagreement_mse_reduction=float(disagreement_gain),
            )
        )
    pl.DataFrame(metrics).sort("rmse").write_csv(out / "d6-model-errors.csv")
    pl.DataFrame(pairs).write_csv(out / "d6-pairwise-diversity.csv")
    for prefix, values in [("prediction", arrays), ("residual", errors)]:
        names = sorted(values)
        correlation = np.corrcoef(np.stack([values[name] for name in names]))
        pl.DataFrame(correlation, schema=names).insert_column(
            0, pl.Series("method", names)
        ).write_csv(out / f"d6-{prefix}-pearson-matrix.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                source_evaluation=cfg["evaluation_run"],
                methods=sorted(arrays),
                n=1457,
                residual_definition="prediction minus measured pIC50",
                weights_fitted=False,
                pair_decomposition="MSE(mean(a,b)) = mean(MSE(a),MSE(b)) - mean((a-b)^2)/4",
                blog_comparison="Published reference 0.503; evaluation cohorts are not matched.",
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
