"""Apply source placement as explicitly assumed test-distribution scenarios."""

import __future__

import ast
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr, spearmanr


class SourceSeries:
    def __init__(self, values):
        self.values = values

    def mean(self):
        return self.values.mean()

    def std(self):
        return self.values.std(ddof=1)


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    syntax = ast.parse((inp / cfg["source"]).read_text())
    names = {"FLOOR", "BLIND_MOMENTS", "SOLVED_PEARSON", "STRAE_MOMENTS", "OOF_TO_BLIND"}
    selected = [
        node
        for node in syntax.body
        if (
            isinstance(node, ast.Assign)
            and any(isinstance(t, ast.Name) and t.id in names for t in node.targets)
        )
        or (isinstance(node, ast.FunctionDef) and node.name == "place")
    ]
    namespace = {}
    exec(
        compile(
            ast.Module(body=selected, type_ignores=[]),
            "frozen_source_place",
            "exec",
            flags=__future__.annotations.compiler_flag,
        ),
        namespace,
    )
    namespace["VALUE_COLUMNS"] = {
        iso: iso + "_pIC50_direct_inhibition" for iso in namespace["BLIND_MOMENTS"]
    }
    test = pl.read_parquet(inp / cfg["test"])
    oof = pl.read_parquet(inp / cfg["oof"]) if "oof" in cfg else None
    if "submitted" in cfg:
        submitted = pl.read_parquet(inp / cfg["submitted"]).unpivot(
            index=["Molecule_Name", "SMILES"], variable_name="endpoint", value_name="sent"
        )
        paired = test.join(submitted, on=["Molecule_Name", "SMILES", "endpoint"], validate="1:1")
        assert paired.height == 3000
        np.testing.assert_array_equal(paired["prediction"], paired["sent"])
    assert test.height == 3000 and test.unique(["compound_id", "endpoint"]).height == 3000
    arrays, source_frame, rho = {}, {}, {}
    for iso, endpoint in namespace["VALUE_COLUMNS"].items():
        table = test.filter(pl.col("endpoint") == endpoint).sort("compound_id")
        assert table.height == 750
        arrays[iso] = table
        source_frame[endpoint] = SourceSeries(table["prediction"].to_numpy().astype(float))
        if oof is not None:
            labeled = oof.filter(pl.col("endpoint") == endpoint)
            rho[iso] = float(pearsonr(labeled["prediction"], labeled["y_true"]).statistic)
    scaled = {iso: r * namespace["OOF_TO_BLIND"][iso] for iso, r in rho.items()}
    scenarios = {"raw": {iso: dict(scale=1.0, offset=0.0) for iso in arrays}}
    for name, correlations, strae in [
        ("author_solved_r2", namespace["SOLVED_PEARSON"], False),
        ("author_solved_strae", namespace["SOLVED_PEARSON"], True),
        ("our_oof_times_author_ratio_r2", scaled, False),
        ("our_oof_times_author_ratio_strae", scaled, True),
    ]:
        if not correlations:
            continue
        scenarios[name] = namespace["place"](source_frame, correlations, strae)
    base = scenarios.get("our_oof_times_author_ratio_strae", scenarios["author_solved_strae"])
    scenarios["d6_only_author_strae"] = {iso: dict(c) for iso, c in scenarios["raw"].items()}
    scenarios["d6_only_author_strae"]["CYP2D6"] = dict(base["CYP2D6"])
    for fraction in [0.25, 0.5, 0.75]:
        scenarios[f"partial_{fraction}"] = {
            iso: dict(scale=1 + fraction * (c["scale"] - 1), offset=fraction * c["offset"])
            for iso, c in base.items()
        }
    d6 = source_frame["CYP2D6_pIC50_direct_inhibition"]
    for center in [3.32, 3.57, 3.82]:
        for spread in [0.675, 0.9, 1.125]:
            scenario = {iso: dict(c) for iso, c in base.items()}
            scenario["CYP2D6"] = dict(offset=center - d6.mean(), scale=spread / d6.std())
            scenarios[f"d6_center{center}_sd{spread}"] = scenario
    predictions, metrics = [], []
    for name, parameters in scenarios.items():
        for iso, table in arrays.items():
            values = table["prediction"].to_numpy().astype(float)
            cal = parameters[iso]
            assert cal["scale"] > 0
            affine = values.mean() + cal["scale"] * (values - values.mean()) + cal["offset"]
            untouched = name == "raw" or (name == "d6_only_author_strae" and iso != "CYP2D6")
            placed = values.copy() if untouched else np.maximum(namespace["FLOOR"], affine)
            if untouched:
                np.testing.assert_array_equal(placed, values)
            assert np.isfinite(placed).all() and (np.diff(placed[np.argsort(values)]) >= 0).all()
            predictions.append(
                table.with_columns(pl.Series("prediction", placed), pl.lit(name).alias("scenario"))
            )
            metrics.append(
                dict(
                    scenario=name,
                    endpoint=iso,
                    n=750,
                    raw_mean=float(values.mean()),
                    raw_sd=float(values.std(ddof=1)),
                    offset=cal["offset"],
                    scale=cal["scale"],
                    intended_mean=float(affine.mean()),
                    intended_sd=float(affine.std(ddof=1)),
                    actual_mean=float(placed.mean()),
                    actual_sd=float(placed.std(ddof=1)),
                    minimum=float(placed.min()),
                    maximum=float(placed.max()),
                    floored=int((affine < namespace["FLOOR"]).sum()) if name != "raw" else 0,
                    rank_correlation_to_raw=float(spearmanr(values, placed).statistic),
                    mean_absolute_change=float(np.abs(placed - values).mean()),
                )
            )
    pl.concat(predictions).write_parquet(out / "scenario-test-predictions.parquet")
    pl.DataFrame(metrics).write_csv(out / "scenario-distributions.csv")
    (out / "assumptions.json").write_text(
        json.dumps(
            dict(
                source_commit=cfg["source_commit"],
                author_constants={name: namespace[name] for name in names},
                measured_oof_pearson=rho,
                assumed_blind_pearson_from_author_ratio=scaled,
                scenario_count=len(scenarios),
                test_labels_used=False,
                submitted=False,
                interpretation="Assumed author moments/ratios; no measured test score.",
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
