"""Refit the CYPMol-ensemble readouts on all labelled rows and blend them into a payload.

For each configured endpoint the member is the mean of the prespecified recipes, each fitted on
every row of the incumbent OOF cohort (no folds) and predicted on the 750 blind compounds. The
payload column is `(1 - w) * base + w * member` at a fixed per-endpoint weight, moment-matched
back to the base column's mean and SD so only the ranking changes (the form evaluated on the OOF
with per-fold moment matching); all other columns are copied from the base payload unchanged. The
fold-averaged screen predictions are kept as a diagnostic of the full-refit blind column.
"""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import regression_embedding_member as screen
import torch
from threadpoolctl import threadpool_limits


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    threadpool_limits(limits=cfg["threads"])
    torch.set_num_threads(cfg["threads"])
    frame = pl.read_parquet(inp / "embeddings.parquet")
    blocks = {}
    for name, block in cfg["feature_blocks"].items():
        cols = [c for c in frame.columns if c.startswith(tuple(block["prefixes"]))]
        blocks[name] = frame.select("compound_id", *cols).rename({c: f"{name}:{c}" for c in cols})
    blind_ids = pl.read_parquet(inp / cfg["blind_compounds"])["compound_id"].sort()
    assert blind_ids.len() == 750
    screen_blind = pl.read_parquet(inp / "screen-blind-mean.parquet")

    members, diagnostics = [], {}
    for endpoint, weight in cfg["weights"].items():
        control = cfg["controls"][endpoint]
        oof = pl.read_parquet(inp / control["file"])
        for column, value in control["filter"].items():
            oof = oof.filter(pl.col(column) == value)
        oof = oof.filter(pl.col("endpoint") == endpoint).sort("compound_id")
        assert oof.height == control["rows"]
        y = oof["y_true"].to_numpy()
        n = oof.height
        ids = pl.concat([oof.select("compound_id"), pl.DataFrame({"compound_id": blind_ids})])
        per_recipe = []
        for recipe in cfg["recipes"]:
            name = f"{recipe['blocks']}--{recipe['selection']}--{recipe['readout']}"
            table = ids
            for block in recipe["blocks"].split("+"):
                table = table.join(blocks[block], on="compound_id", how="left", validate="1:1")
            assert table["compound_id"].equals(ids["compound_id"])
            x = table.drop("compound_id").to_numpy().astype(np.float32)
            x[~np.isfinite(x)] = np.nan
            medians = np.nanmedian(x[:n], axis=0)
            medians[~np.isfinite(medians)] = 0.0
            xx = np.where(np.isfinite(x), x, medians).astype(np.float32)
            seed = cfg["seed"]
            columns = screen.select_columns(xx[:n], y, recipe["selection"], seed, cfg["threads"])
            model = screen.build_model(recipe["readout"], seed, cfg, inp)
            model.fit(xx[:n][:, columns], y)
            prediction = np.asarray(model.predict(xx[n:][:, columns]), dtype=float)
            assert np.isfinite(prediction).all()
            per_recipe.append(prediction)
            print(endpoint, name, float(prediction.mean()), float(prediction.std()), flush=True)
            del model
            torch.cuda.empty_cache()
        member = np.mean(per_recipe, axis=0)
        reference = (
            screen_blind.filter(pl.col("endpoint") == endpoint)
            .group_by("compound_id")
            .agg(pl.col("prediction").mean())
            .join(pl.DataFrame({"compound_id": blind_ids}), on="compound_id", how="semi")
            .sort("compound_id")["prediction"]
            .to_numpy()
        )
        diagnostics[endpoint] = {
            "weight": weight,
            "member_mean": float(member.mean()),
            "member_sd": float(member.std(ddof=1)),
            "corr_with_fold_mean_screen": float(np.corrcoef(member, reference)[0, 1]),
        }
        members.append(
            pl.DataFrame(
                {"compound_id": blind_ids, "endpoint": endpoint, "prediction": member}
            )
        )
    pl.concat(members).write_parquet(out / "member-blind.parquet")

    base = pl.read_csv(inp / "base-predictions.csv")
    mapping = pl.read_parquet(inp / "test-mapping.parquet").select(
        "compound_id", pl.col("molecule_name").alias("Molecule_Name")
    )
    payload = base
    changes = {}
    for endpoint, weight in cfg["weights"].items():
        member = (
            pl.concat(members)
            .filter(pl.col("endpoint") == endpoint)
            .join(mapping, on="compound_id", validate="1:1")
            .select("Molecule_Name", pl.col("prediction").alias("member"))
        )
        joined = payload.join(member, on="Molecule_Name", how="left", validate="1:1")
        assert joined["member"].null_count() == 0
        old = joined[endpoint].to_numpy()
        blended = (1 - weight) * old + weight * joined["member"].to_numpy()
        new = old.mean() + (blended - blended.mean()) / blended.std(ddof=1) * old.std(ddof=1)
        payload = payload.with_columns(pl.Series(endpoint, new))
        changes[endpoint] = {
            "mean_abs": float(np.abs(new - old).mean()),
            "max_abs": float(np.abs(new - old).max()),
            "spearman": float(
                pl.DataFrame({"a": old, "b": new}).select(pl.corr("a", "b", method="spearman")).item()
            ),
            "mean": [float(old.mean()), float(new.mean())],
            "sd": [float(old.std(ddof=1)), float(new.std(ddof=1))],
        }
    for column in base.columns:
        if column.endswith("_direct_inhibition") and column not in cfg["weights"]:
            assert payload[column].equals(base[column])
    payload.write_csv(out / "predictions.csv")
    summary = {"diagnostics": diagnostics, "blind_change": changes}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
