"""Grand-bank nested selection over every deployable atomic regression OOF column, then a payload.

Scans all regression run outputs for OOF tables (compound, endpoint, fold, prediction) that cover an
endpoint's incumbent cohort exactly, matches each to a 750-compound blind column of the same run,
drops derived (blend/calibration/oracle) columns, and runs the nested anchored greedy selection of
`regression_bank_selection` over the incumbent plus the remaining columns. An endpoint is changed in
the payload only if its nested delta passes the prespecified rule in the config or the user
listed it as adopted; the changed
column is the selected blend moment-matched back to the base column's mean and SD.
"""

import hashlib
import json
import os
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl
import regression_bank_selection as sel

SKIP = {
    "std",
    "value",
    "measurement_id",
    "expected_fold",
    "is_neighbor",
    "fold",
    "split_id",
    "y_true",
    "conf_low",
    "conf_high",
    "prediction",
    "compound_id",
    "endpoint",
}
BAD = re.compile(r"observed|oracle|measured|truth|diagnostic", re.I)
DERIVED = re.compile(
    r"caruana|equal|median|ridge|simplex|positive_linear|best_single|community|previous_global|"
    r"affine|offset|calibrat|blend|ensemble|interval_strength|neighbor|selected_only|union|stack|"
    r"crossfit|full_fit|aux8|^full$",
    re.I,
)


def group_columns(frame, extra):
    return [
        c
        for c in extra
        if c not in SKIP
        and c in frame.columns
        and frame[c].dtype in (pl.Utf8, pl.Boolean, pl.Int64, pl.Float64)
        and frame[c].n_unique() < 300
    ]


def partitions(frame, columns):
    if not columns:
        return {(): frame}
    return frame.partition_by(columns, as_dict=True)


def fingerprint(vector):
    return hashlib.md5(np.round(vector, 6).tobytes()).hexdigest()


def list_tables(root, need):
    found = []
    for run in sorted(os.listdir(root)):
        out = root / run / "outputs"
        if not out.is_dir():
            continue
        for directory, _, files in os.walk(out):
            for name in files:
                path = Path(directory) / name
                if name.endswith(".parquet") and path.stat().st_size < 200e6:
                    try:
                        schema = pl.read_parquet_schema(path)
                    except Exception:
                        continue
                    if need <= set(schema):
                        found.append((run, path, sorted(set(schema) - SKIP)))
    return found


def collect_blind(root, test_ids):
    blind = defaultdict(list)
    for run, path, extra in list_tables(root, {"compound_id", "endpoint", "prediction"}):
        frame = pl.read_parquet(path, columns=None)
        if frame.select(pl.col("compound_id").is_in(test_ids).mean().fill_null(0)).item() < 0.99:
            continue
        columns = group_columns(frame, extra)
        for endpoint in frame["endpoint"].unique().to_list():
            sub = frame.filter(pl.col("endpoint") == endpoint)
            for key, part in partitions(sub, columns).items():
                key = key if isinstance(key, tuple) else (key,)
                mean = (
                    part.group_by("compound_id")
                    .agg(pl.col("prediction").mean())
                    .sort("compound_id")
                )
                if mean["compound_id"].to_list() != test_ids:
                    continue
                vector = mean["prediction"].to_numpy()
                if np.isfinite(vector).all():
                    blind[(run, endpoint)].append(
                        (path.name, dict(zip(columns, map(str, key), strict=True)), vector)
                    )
    return blind


def collect_oof(root, cohorts, blind):
    bank = {endpoint: [] for endpoint in cohorts}
    seen = defaultdict(set)
    need = {"compound_id", "endpoint", "prediction", "y_true", "fold"}
    for run, path, extra in list_tables(root, need):
        frame = pl.read_parquet(path)
        columns = group_columns(frame, extra)
        for endpoint, cohort in cohorts.items():
            sub = frame.filter(pl.col("endpoint") == endpoint)
            if sub.is_empty():
                continue
            for key, part in partitions(sub, columns).items():
                if part.height < cohort.height:
                    continue
                key = key if isinstance(key, tuple) else (key,)
                part = part.unique(subset=["compound_id", "fold"], keep="first")
                joined = cohort.join(
                    part.select("compound_id", "fold", "prediction"),
                    on=["compound_id", "fold"],
                    how="inner",
                    validate="1:1",
                ).sort("compound_id")
                if joined.height != cohort.height:
                    continue
                vector = joined["prediction"].to_numpy()
                if not np.isfinite(vector).all() or fingerprint(vector) in seen[endpoint]:
                    continue
                seen[endpoint].add(fingerprint(vector))
                meta = dict(zip(columns, map(str, key), strict=True))
                matches = {
                    fingerprint(v): (name, v)
                    for name, bkey, v in blind.get((run, endpoint), [])
                    if (
                        set(bkey) & set(meta)
                        and all(bkey[c] == meta[c] for c in set(bkey) & set(meta))
                    )
                    or (not bkey and not meta)
                }
                if len(matches) != 1:
                    continue
                blind_name, blind_vector = next(iter(matches.values()))
                label = f"{run}|{path.name}|{json.dumps(meta, sort_keys=True)}"
                if BAD.search(label) or DERIVED.search(label.split("|", 1)[1]):
                    continue
                bank[endpoint].append(
                    {"label": label, "blind_file": blind_name, "oof": vector, "blind": blind_vector}
                )
    return bank


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    root = Path(cfg["runs_root"])
    rng = np.random.default_rng(cfg["seed"])
    test_map = pl.read_parquet(inp / "test-mapping.parquet").sort("compound_id")
    test_ids = test_map["compound_id"].to_list()

    cohorts, truth = {}, {}
    for endpoint, control in cfg["controls"].items():
        frame = pl.read_parquet(inp / control["file"])
        for column, value in control["filter"].items():
            frame = frame.filter(pl.col(column) == value)
        frame = frame.filter(pl.col("endpoint") == endpoint).sort("compound_id")
        assert frame.height == control["rows"]
        cohorts[endpoint] = frame.select("compound_id", "fold")
        truth[endpoint] = frame
    bank = collect_oof(root, cohorts, collect_blind(root, test_ids))

    base = pl.read_csv(inp / "base-predictions.csv")
    mapping = test_map.select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    payload = base
    report, manifest = [], {}
    for endpoint, frame in truth.items():
        y, low, high, folds = (
            frame[c].to_numpy() for c in ("y_true", "conf_low", "conf_high", "fold")
        )
        incumbent = frame["prediction"].to_numpy()
        pool = bank[endpoint]
        x = np.column_stack([incumbent] + [m["oof"] for m in pool])
        baseline = sel.st_rae(y, low, high, incumbent)
        matched = np.empty(len(y))
        for outer in sorted(set(folds)):
            fit, held = folds != outer, folds == outer
            weights, _ = sel.procedure_weights(
                x[fit], y[fit], low[fit], high[fit], folds[fit], "greedy_anchored", cfg
            )
            matched[held] = sel.moment_match(x[held] @ weights, incumbent[held])
        deltas = []
        for _ in range(cfg["bootstrap"]):
            pick = rng.integers(0, len(y), len(y))
            deltas.append(
                sel.st_rae(y[pick], low[pick], high[pick], matched[pick])
                - sel.st_rae(y[pick], low[pick], high[pick], incumbent[pick])
            )
        wins = sum(
            sel.st_rae(y[folds == f], low[folds == f], high[folds == f], matched[folds == f])
            < sel.st_rae(y[folds == f], low[folds == f], high[folds == f], incumbent[folds == f])
            for f in sorted(set(folds))
        )
        delta = sel.st_rae(y, low, high, matched) - baseline
        rule = cfg["rule"]
        by_rule = delta <= rule["max_delta"] and wins >= rule["min_wins"]
        by_user = endpoint in cfg.get("user_adopted", [])
        adopt = by_rule or by_user
        weights, steps = sel.procedure_weights(x, y, low, high, folds, "greedy_anchored", cfg)
        members = {
            ("incumbent" if j == 0 else pool[j - 1]["label"]): round(float(w), 4)
            for j, w in enumerate(weights)
            if w > 0
        }
        record = {
            "endpoint": endpoint,
            "pool": len(pool),
            "incumbent_st_rae": baseline,
            "nested_delta": delta,
            "ci": [float(np.quantile(deltas, 0.025)), float(np.quantile(deltas, 0.975))],
            "folds_won": int(wins),
            "adopted": bool(adopt),
            "adopted_by": "rule" if by_rule else ("user" if by_user else None),
            "steps": steps,
            "members": members,
        }
        if adopt:
            blind_columns = np.column_stack(
                [pool[j - 1]["blind"] for j in range(1, len(weights)) if weights[j] > 0]
            )
            used = [j for j in range(1, len(weights)) if weights[j] > 0]
            member_blind = blind_columns @ weights[used]
            current = payload.join(
                mapping.join(pl.DataFrame({"compound_id": test_ids}), on="compound_id"),
                on="Molecule_Name",
                how="left",
            ).sort("compound_id")
            old = current[endpoint].to_numpy()
            new_sorted = old * weights[0] + member_blind
            new_sorted = old.mean() + (new_sorted - new_sorted.mean()) / new_sorted.std(
                ddof=1
            ) * old.std(ddof=1)
            replaced = current.select("Molecule_Name").with_columns(pl.Series("new", new_sorted))
            payload = payload.join(replaced, on="Molecule_Name", how="left", validate="1:1")
            payload = payload.with_columns(pl.col("new").alias(endpoint)).drop("new")
            record["blind_change"] = {
                "mean_abs": float(np.abs(new_sorted - old).mean()),
                "max_abs": float(np.abs(new_sorted - old).max()),
            }
        report.append(record)
        manifest[endpoint] = [m["label"] + "|" + m["blind_file"] for m in pool]
        print(json.dumps({k: v for k, v in record.items() if k != "members"}), flush=True)
    for column in base.columns:
        if column.endswith("_direct_inhibition") and column not in {
            r["endpoint"] for r in report if r["adopted"]
        }:
            assert payload[column].equals(base[column])
    payload.write_csv(out / "predictions.csv")
    (out / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    (out / "bank-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
