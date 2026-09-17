"""Diagnose Jeremy-medium transfer gains without fitting or changing models."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    exact = inp / name
    paths = (
        [exact]
        if exact.exists()
        else [path for path in inp.iterdir() if path.name.partition("-")[2] == name]
    )
    if len(paths) != 1:
        raise ValueError(f"Expected one input for {name}: {paths}")
    return pl.read_parquet(paths[0])


spec = importlib.util.spec_from_file_location("official", inp / "custom_scoring_functions.py")
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)

compounds = read("parent-assignments.parquet").select("compound_id", "connectivity")
overlap = read("challenge-overlap-observations.parquet")
overlap_sets = {
    "public_overlap_any": set(overlap["connectivity"]),
    "tdc_overlap": set(overlap.filter(pl.col("source") == "tdc_adme")["connectivity"]),
    "polaris_overlap": set(
        overlap.filter(pl.col("source") == "polaris_novartis_cyp")["connectivity"]
    ),
}
same_isoform_sets = {
    isoform: set(overlap.filter(pl.col("isoform") == isoform)["connectivity"])
    for isoform in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
}

jeremy = read("jeremy-oof.parquet").join(compounds, on="compound_id", validate="m:1")
rng = np.random.default_rng(cfg["bootstrap_seed"])
metric_rows = []
bootstrap_rows = []
keys = ["compound_id", "endpoint", "fold"]
label_columns = ["y_true", "conf_low", "conf_high"]


def metric(frame: pl.DataFrame, prediction: str) -> float:
    return float(
        scorer.rae_soft_threshold_absolute_error(
            frame["y_true"].to_numpy(),
            frame[prediction].to_numpy(),
            y_true_lower=frame["conf_low"].to_numpy(),
            y_true_upper=frame["conf_high"].to_numpy(),
        )
    )


for endpoint in cfg["endpoints"]:
    isoform = endpoint.split("_", 1)[0]
    owned = read(cfg["owned_by_endpoint"][endpoint]).filter(pl.col("endpoint") == endpoint)
    table = (
        jeremy.filter(pl.col("endpoint") == endpoint)
        .rename({"prediction": "jeremy_prediction"})
        .join(
            owned.select(*keys, pl.col("prediction").alias("owned_prediction")),
            on=keys,
            validate="1:1",
        )
        .sort(keys)
    )
    subset_sets = {
        "all": None,
        "no_public_overlap": overlap_sets["public_overlap_any"],
        **overlap_sets,
        "same_isoform_overlap": same_isoform_sets[isoform],
    }
    for subset, connectivities in subset_sets.items():
        if subset == "all":
            selected = table
        elif subset == "no_public_overlap":
            selected = table.filter(~pl.col("connectivity").is_in(connectivities))
        else:
            selected = table.filter(pl.col("connectivity").is_in(connectivities))
        if selected.height < 10:
            continue
        jeremy_value = metric(selected, "jeremy_prediction")
        owned_value = metric(selected, "owned_prediction")
        row = {
            "endpoint": endpoint,
            "subset": subset,
            "n": selected.height,
            "unique_compounds": selected["compound_id"].n_unique(),
            "jeremy_ST_RAE": jeremy_value,
            "owned_ST_RAE": owned_value,
            "owned_minus_jeremy": owned_value - jeremy_value,
            "prediction_pearson": float(
                selected.select(pl.corr("jeremy_prediction", "owned_prediction")).item()
            ),
        }
        metric_rows.append(row)
        if subset in {"all", "public_overlap_any", "no_public_overlap"}:
            values = []
            for _ in range(cfg["bootstrap_replicates"]):
                indices = rng.integers(0, selected.height, selected.height)
                sample = selected[indices]
                values.append(
                    metric(sample, "owned_prediction") - metric(sample, "jeremy_prediction")
                )
            bootstrap_rows.append(
                {
                    "endpoint": endpoint,
                    "subset": subset,
                    "n": selected.height,
                    "delta_mean": float(np.mean(values)),
                    "delta_p025": float(np.quantile(values, 0.025)),
                    "delta_p975": float(np.quantile(values, 0.975)),
                }
            )


def geometry(name: str, frame: pl.DataFrame) -> tuple[dict, np.ndarray]:
    columns = sorted(column for column in frame.columns if column.startswith("emb_"))
    aligned = compounds.select("compound_id").join(
        frame.select("compound_id", *columns), on="compound_id", validate="1:1"
    )
    x = aligned.select(columns).to_numpy().astype(np.float64)
    x -= x.mean(axis=0, keepdims=True)
    std = x.std(axis=0)
    covariance = x.T @ x / max(1, x.shape[0] - 1)
    eigenvalues = np.linalg.eigvalsh(covariance)[::-1]
    eigenvalues = np.maximum(eigenvalues, 0)
    proportions = eigenvalues / eigenvalues.sum()
    entropy = -np.sum(proportions[proportions > 0] * np.log(proportions[proportions > 0]))
    cumulative = np.cumsum(proportions)
    result = {
        "name": name,
        "rows": x.shape[0],
        "dimensions": x.shape[1],
        "dead_dimensions": int(np.sum(std < 1e-8)),
        "median_feature_std": float(np.median(std)),
        "median_centered_norm": float(np.median(np.linalg.norm(x, axis=1))),
        "effective_rank": float(np.exp(entropy)),
        "participation_ratio": float(1 / np.sum(proportions**2)),
    }
    for threshold in [0.5, 0.8, 0.9, 0.95, 0.99]:
        result[f"pcs_{int(threshold * 100)}"] = int(np.searchsorted(cumulative, threshold) + 1)
    return result, x


geometry_rows = []
centered = {}
for item in cfg["geometry_embeddings"]:
    row, x = geometry(item["name"], read(item["file"]))
    geometry_rows.append(row)
    centered[item["name"]] = x

cka_rows = []
names = list(centered)
for left_index, left in enumerate(names):
    x = centered[left]
    x_norm = np.linalg.norm(x.T @ x)
    for right in names[left_index + 1 :]:
        y = centered[right]
        cka_rows.append(
            {
                "left": left,
                "right": right,
                "linear_CKA": float(
                    np.linalg.norm(x.T @ y) ** 2 / (x_norm * np.linalg.norm(y.T @ y))
                ),
            }
        )

metrics = pl.DataFrame(metric_rows).sort("endpoint", "subset")
bootstraps = pl.DataFrame(bootstrap_rows).sort("endpoint", "subset")
geometry_frame = pl.DataFrame(geometry_rows)
cka_frame = pl.DataFrame(cka_rows)
metrics.write_csv(out / "overlap-metrics.csv")
bootstraps.write_csv(out / "overlap-bootstrap.csv")
geometry_frame.write_csv(out / "embedding-geometry.csv")
cka_frame.write_csv(out / "embedding-cka.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "overlap_connectivities": {name: len(values) for name, values in overlap_sets.items()},
            "same_isoform_connectivities": {
                name: len(values) for name, values in same_isoform_sets.items()
            },
            "metrics": metrics.to_dicts(),
            "bootstraps": bootstraps.to_dicts(),
            "geometry": geometry_frame.to_dicts(),
            "cka": cka_frame.to_dicts(),
        },
        indent=2,
    )
    + "\n"
)
print(metrics)
print(bootstraps)
print(geometry_frame)
print(cka_frame)
