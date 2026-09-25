"""Compare depth members against exact submitted 1A2, 2C9 and 3A4 OOF recipes."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr


def load_scorer(path):
    spec = importlib.util.spec_from_file_location("official_scoring", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    scorer = load_scorer(inp / cfg["scorer"])
    keys = ["compound_id", "endpoint", "fold"]
    all_metrics, all_predictions, weight_rows = [], [], []
    for target, setup in cfg["targets"].items():
        endpoint = setup["endpoint"]
        summary = json.loads((inp / setup["bank_summary"]).read_text())
        component_specs = summary["config"]["components"]
        names = [item["name"] for item in component_specs]
        bank = (
            pl.read_parquet(inp / setup["bank_oof"])
            .filter(pl.col("method") == setup["method"])
            .sort(keys)
        )
        assert bank.height == setup["expected_rows"]
        table = bank.select(*keys, "y_true", "conf_low", "conf_high")
        for item in component_specs:
            alias = f"{target}-{item['name']}.parquet"
            component = pl.read_parquet(inp / alias).filter(pl.col("endpoint") == endpoint)
            assert component.height == table.height
            table = table.join(
                component.select(*keys, pl.col("prediction").alias(item["name"])),
                on=keys,
                validate="1:1",
            ).sort(keys)
        folds = table["fold"].to_numpy()
        weight_matrix = np.empty((table.height, len(names)), dtype=float)
        if setup["method"] == "equal_all":
            weight_matrix[:] = 1 / len(names)
        else:
            fits = pl.read_csv(inp / setup["bank_fits"]).filter(pl.col("method") == setup["method"])
            assert fits.height == 5
            for row in fits.iter_rows(named=True):
                weights = np.asarray(json.loads(row["weights"]), dtype=float)
                assert len(weights) == len(names) and np.isclose(weights.sum(), 1)
                weight_matrix[folds == row["outer_fold"]] = weights
        x = table.select(names).to_numpy()
        control = np.sum(x * weight_matrix, axis=1)
        saved = bank["prediction"].to_numpy()
        max_difference = float(np.max(np.abs(control - saved)))
        assert max_difference < 1e-9, (target, max_difference)
        for fold in range(5):
            mask = folds == fold
            weights = weight_matrix[mask][0]
            for name, value in zip(names, weights, strict=True):
                weight_rows.append(
                    {"target": target, "fold": fold, "member": name, "weight": float(value)}
                )

        def aligned(filename, name, depth=None, endpoint=endpoint, table=table, target=target):
            frame = pl.read_parquet(inp / filename).filter(pl.col("endpoint") == endpoint)
            if depth is not None:
                frame = frame.filter(pl.col("depth") == depth)
            assert frame.height == table.height, (target, name, frame.height)
            joined = (
                table.select(*keys)
                .join(
                    frame.select(*keys, pl.col("prediction").alias(name)), on=keys, validate="1:1"
                )
                .sort(keys)
            )
            assert joined.height == table.height
            return joined[name].to_numpy()

        candidates = {}
        for depth in (1, 3, 6):
            candidates[f"chemeleon_d{depth}_lgbm"] = aligned(
                cfg["chemeleon_depths"], f"chemeleon_d{depth}_lgbm", depth
            )
        candidates["chemeleon_d1_tabpfn"] = aligned(cfg["chemeleon_tabpfn"], "chemeleon_d1_tabpfn")
        candidates["chemeleon_d1_tabicl"] = aligned(cfg["chemeleon_tabicl"], "chemeleon_d1_tabicl")
        for depth in (1, 2, 4):
            candidates[f"jeremy_d{depth}_lgbm"] = aligned(
                setup["jeremy_oof"], f"jeremy_d{depth}_lgbm", depth
            )
        assert (
            np.max(np.abs(candidates["jeremy_d4_lgbm"] - table[setup["ft_member"]].to_numpy())) == 0
        )
        candidates["chemeleon_1_6"] = (
            candidates["chemeleon_d1_lgbm"] + candidates["chemeleon_d6_lgbm"]
        ) / 2
        candidates["chemeleon_1_3_6"] = np.mean(
            [candidates[f"chemeleon_d{depth}_lgbm"] for depth in (1, 3, 6)], axis=0
        )
        candidates["jeremy_1_4"] = (candidates["jeremy_d1_lgbm"] + candidates["jeremy_d4_lgbm"]) / 2
        candidates["jeremy_1_2_4"] = np.mean(
            [candidates[f"jeremy_d{depth}_lgbm"] for depth in (1, 2, 4)], axis=0
        )
        recipes = {"current": control}
        for name in (
            "chemeleon_d1_lgbm",
            "chemeleon_d1_tabpfn",
            "chemeleon_d1_tabicl",
            "chemeleon_1_6",
            "chemeleon_1_3_6",
            "jeremy_d1_lgbm",
            "jeremy_d2_lgbm",
            "jeremy_1_4",
            "jeremy_1_2_4",
        ):
            candidate = candidates[name]
            if setup["method"] != "caruana":
                active = (weight_matrix > 0).sum(axis=1)
                assert np.all(active == active[0])
                recipes[f"add_{name}"] = (active * control + candidate) / (active + 1)
            replaced = (
                setup["ft_member"] if name.startswith("jeremy") else setup["chemeleon_member"]
            )
            index = names.index(replaced)
            recipes[f"replace_{replaced}__{name}"] = control + weight_matrix[:, index] * (
                candidate - x[:, index]
            )
        y = table["y_true"].to_numpy()
        low = table["conf_low"].to_numpy()
        high = table["conf_high"].to_numpy()
        for name, prediction in recipes.items():
            all_predictions.append(
                table.select(*keys).with_columns(
                    pl.lit(target).alias("target"),
                    pl.lit(name).alias("scenario"),
                    pl.Series("prediction", prediction),
                )
            )
            for fold in [-1, 0, 1, 2, 3, 4]:
                mask = np.ones(len(y), dtype=bool) if fold == -1 else folds == fold
                all_metrics.append(
                    {
                        "target": target,
                        "scenario": name,
                        "fold": fold,
                        "n": int(mask.sum()),
                        "ST_RAE": float(
                            scorer.rae_soft_threshold_absolute_error(
                                y[mask],
                                prediction[mask],
                                y_true_lower=low[mask],
                                y_true_upper=high[mask],
                            )
                        ),
                        "MAE": float(np.mean(np.abs(y[mask] - prediction[mask]))),
                        "Spearman": float(spearmanr(y[mask], prediction[mask]).statistic),
                    }
                )
    pl.DataFrame(all_metrics).write_csv(out / "metrics.csv")
    pl.concat(all_predictions).write_parquet(out / "oof.parquet")
    pl.DataFrame(weight_rows).write_csv(out / "current-member-weights.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "rule": (
                    "Reproduce current bank exactly; replace an atomic member at its saved "
                    "fold weight, or add one equal member to an equal-mean bank. "
                    "Caruana additions require refitting selection separately."
                ),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
