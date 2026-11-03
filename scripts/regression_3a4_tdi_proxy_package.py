"""Swap the two anchor59 TabPFN 3A4 members for TDI-proxy refits at fixed weights.

Development evidence keeps the deployed mean Caruana weights and replaces only the two
member OOF columns, comparing original files, matched baseline refits and proxy refits
on the full incumbent cohort and on the labelled RL0191 analogs. The blind payload is
the last sent regression payload with only its 3A4 column rebuilt: the weighted atomic
sum with the two members swapped, moment-matched to the sent 3A4 mean and SD.
"""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

ENDPOINT = "CYP3A4_pIC50_direct_inhibition"


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    spec = importlib.util.spec_from_file_location("official", inp / "custom_scoring_functions.py")
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)

    def st_rae(frame, prediction):
        return float(
            scorer.rae_soft_threshold_absolute_error(
                frame["y_true"].to_numpy(),
                prediction,
                y_true_upper=frame["conf_high"].to_numpy(),
                y_true_lower=frame["conf_low"].to_numpy(),
            )
        )

    weights = cfg["weights"]
    assert abs(sum(weights.values()) - 1) < 1e-9
    cohort = (
        pl.read_parquet(inp / "control-3a4.parquet")
        .filter((pl.col("method") == "caruana") & (pl.col("endpoint") == ENDPOINT))
        .sort("compound_id")
        .select("compound_id", "fold", "y_true", "conf_low", "conf_high")
    )
    assert cohort.height == 2300

    def column(file, method=None):
        frame = pl.read_parquet(inp / file).filter(pl.col("endpoint") == ENDPOINT)
        if method is not None:
            frame = frame.filter(pl.col("method") == method)
        frame = cohort.select("compound_id").join(
            frame.select("compound_id", "prediction"), on="compound_id", how="left"
        )
        values = frame["prediction"].to_numpy()
        assert np.isfinite(values).all()
        return values

    members = {name: column(file) for name, file in cfg["member_oof"].items()}
    swaps = cfg["swapped_members"]
    variants = {"original": {}}
    for arm in ("baseline", "proxy"):
        variants[arm] = {
            name: column(spec_["file"], f"{spec_['method']}--{arm}")
            for name, spec_ in swaps.items()
        }
    analogs = set(pl.read_parquet(inp / "analog-pairs.parquet")["analog_id"])
    analog_mask = cohort["compound_id"].is_in(list(analogs)).to_numpy()
    folds = cohort["fold"].to_numpy()

    def compose(replacements):
        return sum(weights[name] * replacements.get(name, members[name]) for name in weights)

    rows = []
    composed = {name: compose(values) for name, values in variants.items()}
    for name, prediction in composed.items():
        row = {
            "variant": name,
            "full": st_rae(cohort, prediction),
            "analog": st_rae(cohort.filter(pl.Series(analog_mask)), prediction[analog_mask]),
        }
        for fold in range(5):
            mask = folds == fold
            row[f"fold{fold}"] = st_rae(cohort.filter(pl.Series(mask)), prediction[mask])
        rows.append(row)
    metrics = pl.DataFrame(rows)
    metrics.write_parquet(out / "development-metrics.parquet")
    print(metrics)

    rng = np.random.default_rng(cfg["seed"])
    groups = pl.read_parquet(inp / "parent-assignments.parquet").select("compound_id", "group")
    group_ids = cohort.join(groups, on="compound_id", how="left")["group"].to_numpy()
    unique = np.unique(group_ids)
    index = {g: np.flatnonzero(group_ids == g) for g in unique}
    deltas = []
    for _ in range(cfg["bootstrap"]):
        rows_ = np.concatenate([index[g] for g in rng.choice(unique, len(unique))])
        sample = cohort[rows_]
        deltas.append(
            st_rae(sample, composed["proxy"][rows_]) - st_rae(sample, composed["baseline"][rows_])
        )
    interval = np.quantile(deltas, [0.025, 0.5, 0.975]).tolist()

    atomic = pl.read_parquet(inp / "atomic-test.parquet").filter(pl.col("endpoint") == ENDPOINT)
    assert set(atomic["member"]) == set(weights)
    blind = atomic.pivot(on="member", index="compound_id", values="prediction").sort("compound_id")
    for name, file in cfg["proxy_test"].items():
        new = pl.read_parquet(inp / file).filter(pl.col("endpoint") == ENDPOINT)
        joined = blind.select("compound_id").join(
            new.select("compound_id", "prediction"), on="compound_id", how="left"
        )
        assert joined["prediction"].null_count() == 0
        blind = blind.with_columns(pl.Series(name, joined["prediction"].to_numpy()))
    raw = sum(weights[name] * blind[name].to_numpy() for name in weights)

    sent = pl.read_csv(inp / "sent-predictions.csv")
    mapping = pl.read_parquet(inp / "test-mapping.parquet").select(
        "compound_id", pl.col("molecule_name").alias("Molecule_Name")
    )
    new_column = (
        pl.DataFrame({"compound_id": blind["compound_id"], "raw": raw})
        .join(mapping, on="compound_id", validate="1:1")
        .select("Molecule_Name", "raw")
    )
    payload = sent.join(new_column, on="Molecule_Name", how="left", validate="1:1")
    assert payload["raw"].null_count() == 0 and payload.height == 750
    target = payload[ENDPOINT].to_numpy()
    values = payload["raw"].to_numpy()
    placed = target.mean() + (values - values.mean()) / values.std(ddof=1) * target.std(ddof=1)
    payload = payload.with_columns(pl.Series(ENDPOINT, placed)).drop("raw")
    payload.write_csv(out / "predictions.csv")
    changed = sent.select("Molecule_Name", pl.col(ENDPOINT).alias("old")).join(
        payload.select("Molecule_Name", pl.col(ENDPOINT).alias("new")), on="Molecule_Name"
    )
    for other in [c for c in sent.columns if c.endswith("_direct_inhibition") and c != ENDPOINT]:
        assert payload[other].equals(sent[other])
    old, new = changed["old"].to_numpy(), changed["new"].to_numpy()
    summary = {
        "development": metrics.to_dicts(),
        "proxy_minus_baseline_group_bootstrap": interval,
        "blind_change": {
            "mean_abs": float(np.abs(new - old).mean()),
            "max_abs": float(np.abs(new - old).max()),
            "spearman": float(
                pl.DataFrame({"a": old, "b": new})
                .select(pl.corr("a", "b", method="spearman"))
                .item()
            ),
            "mean": [float(old.mean()), float(new.mean())],
            "sd": [float(old.std(ddof=1)), float(new.std(ddof=1))],
        },
        "unchanged_columns": [c for c in sent.columns if c != ENDPOINT],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
