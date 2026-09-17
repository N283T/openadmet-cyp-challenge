"""Package an all-target candidate from target-specific flat atomic banks."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr

from cyp_challenge.tracks import validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

aliases = (
    pl.read_parquet(inp / cfg["aliases"])
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
reference = pl.read_csv(inp / cfg["reference"])
test = reference.select("SMILES", "Molecule_Name").join(
    aliases, on="Molecule_Name", validate="1:1", maintain_order="left"
)
assert test.height == test["compound_id"].n_unique() == 750
ids = test.select("compound_id").sort("compound_id")


def long_prediction(filename: str, endpoint: str) -> np.ndarray:
    frame = (
        pl.read_parquet(inp / filename)
        .filter(pl.col("endpoint") == endpoint)
        .select("compound_id", "prediction")
    )
    aligned = ids.join(frame, on="compound_id", validate="1:1", maintain_order="left")
    assert aligned["prediction"].null_count() == 0
    return aligned["prediction"].to_numpy()


sources = {**cfg["old_predictions"], **cfg["new_predictions"]}
endpoints = list(cfg["recipes"])
component_predictions = {
    endpoint: {name: long_prediction(filename, endpoint) for name, filename in sources.items()}
    for endpoint in endpoints
}
component_predictions["CYP2C9_pIC50_direct_inhibition"]["owned_broad5_2c9_specialist_lgbm"] = (
    long_prediction(cfg["owned_2c9_prediction"], "CYP2C9_pIC50_direct_inhibition")
)

d6_atomic = pl.read_parquet(inp / cfg["d6_atomic_prediction"])
d6_aligned = ids.join(d6_atomic.drop("Molecule_Name"), on="compound_id", validate="1:1")
for name in d6_atomic.columns[2:]:
    component_predictions["CYP2D6_pIC50_direct_inhibition"][name] = d6_aligned[name].to_numpy()

# Recover the two legacy 3A4 atomic model predictions separately, so the final
# ensemble remains flat rather than treating their prior aggregate as a member.
compounds = pl.read_parquet(inp / cfg["caruana_compounds"]).sort("compound_id")
test_mapping = pl.read_parquet(inp / cfg["caruana_test_mapping"])
test_ids = test_mapping.select("compound_id").unique().sort("compound_id")
assert test_ids["compound_id"].equals(ids["compound_id"])
class_table = pl.read_parquet(inp / cfg["caruana_class_features"])
class_columns = sorted(c for c in class_table.columns if c != "compound_id")
class_features = (
    compounds.select("compound_id")
    .join(class_table, on="compound_id", validate="1:1", maintain_order="left")
    .select(class_columns)
    .to_numpy()
    .astype(np.float32)
)
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan_features = np.stack(
    [fpgen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in compounds["standardized_smiles"]]
)
index = {compound_id: i for i, compound_id in enumerate(compounds["compound_id"])}
itest = np.array([index[c] for c in ids["compound_id"]])
a4 = component_predictions["CYP3A4_pIC50_direct_inhibition"]
a4["class_quantile_065_lgbm"] = lgb.Booster(model_file=str(inp / cfg["a4_class_model"])).predict(
    class_features[itest]
)
a4["morgan_quantile_065_lgbm"] = lgb.Booster(model_file=str(inp / cfg["a4_morgan_model"])).predict(
    morgan_features[itest]
)


def cv_mean_weights(
    fits_file: str, bank_config_file: str, method: str = "simplex_interval"
) -> dict[str, float]:
    fits = pl.read_csv(inp / fits_file).filter(pl.col("method") == method)
    assert fits.height == 5 and sorted(fits["outer_fold"].to_list()) == list(range(5))
    names = [
        item["name"] for item in json.loads((inp / bank_config_file).read_text())["components"]
    ]
    weights = np.mean(np.array([json.loads(value) for value in fits["weights"]]), axis=0)
    weights = weights / weights.sum()
    return {
        name: float(weight) for name, weight in zip(names, weights, strict=True) if weight > 1e-12
    }


weights = {
    "CYP1A2_pIC50_direct_inhibition": {
        name: 0.1
        for name in [
            "jeremy_frozen_tabpfn26",
            "jeremy_scores_log2fc_tabpfn26",
            "jeremy_finetuned_log2fc_lgbm",
            "chemeleon_pca256_log2fc_tabicl",
            "chemeleon_raw2048_log2fc_tabicl",
            "tabpfn35_fast_pub030_raw",
            "tabpfn35_full_pub030_raw",
            "tabpfn35_fast_jeremy_raw",
            "tabpfn35_full_jeremy_raw",
            "tabpfn35_fast_anchor59_raw_scores_logfc",
        ]
    },
    "CYP2C9_pIC50_direct_inhibition": {
        name: 0.2
        for name in [
            "jeremy_scores_log2fc_tabpfn26",
            "jeremy_finetuned_log2fc_lgbm",
            "chemeleon_pca256_log2fc_tabicl",
            "owned_broad5_2c9_specialist_lgbm",
            "tabpfn35_full_pub030_raw",
        ]
    },
    "CYP2D6_pIC50_direct_inhibition": cv_mean_weights(
        cfg["d6_fits"],
        cfg["d6_bank_config"],
        cfg["recipes"]["CYP2D6_pIC50_direct_inhibition"].get("weight_method", "simplex_interval"),
    ),
    "CYP3A4_pIC50_direct_inhibition": cv_mean_weights(
        cfg["a4_fits"],
        cfg["a4_bank_config"],
        cfg["recipes"]["CYP3A4_pIC50_direct_inhibition"].get("weight_method", "simplex_interval"),
    ),
}

predictions = {}
atomic_rows = []
for endpoint, endpoint_weights in weights.items():
    available = component_predictions[endpoint]
    missing = set(endpoint_weights) - set(available)
    assert not missing, (endpoint, missing)
    prediction = sum(endpoint_weights[name] * available[name] for name in endpoint_weights)
    assert np.isfinite(prediction).all()
    predictions[endpoint] = prediction
    for name, weight in endpoint_weights.items():
        atomic_rows.append(
            ids.with_columns(
                pl.lit(endpoint).alias("endpoint"),
                pl.lit(name).alias("member"),
                pl.lit(weight).alias("weight"),
                pl.Series("prediction", available[name]),
            )
        )

candidate = test.sort("compound_id")
for endpoint in endpoints:
    candidate = candidate.with_columns(pl.Series(endpoint, predictions[endpoint]))
candidate = candidate.sort("Molecule_Name").select("SMILES", "Molecule_Name", *endpoints)
candidate.write_parquet(out / "candidate.parquet")
candidate.write_csv(out / "candidate.csv")
pl.concat(atomic_rows).write_parquet(out / "atomic-test-predictions.parquet")

placed = candidate
placement = {}
base = pl.read_csv(inp / cfg["base_candidate"]).sort("Molecule_Name")
for endpoint in endpoints:
    raw = candidate[endpoint].to_numpy()
    old = base[endpoint].to_numpy()
    adjusted = (raw - raw.mean()) / raw.std() * old.std() + old.mean()
    placed = placed.with_columns(pl.Series(endpoint, adjusted))
    placement[endpoint] = {
        "source_mean": float(raw.mean()),
        "source_std": float(raw.std()),
        "target_mean": float(old.mean()),
        "target_std": float(old.std()),
    }
placed.write_parquet(out / "candidate-placed.parquet")
placed.write_csv(out / "candidate-placed.csv")

for path in [
    out / "candidate.parquet",
    out / "candidate.csv",
    out / "candidate-placed.parquet",
    out / "candidate-placed.csv",
]:
    report = validate(path, "regression", inp / cfg["reference"])
    if not report["valid"]:
        raise ValueError(report)

comparison = []
for endpoint in endpoints:
    old = base[endpoint].to_numpy()
    new = candidate[endpoint].to_numpy()
    assert not np.array_equal(old, new)
    comparison.append(
        {
            "endpoint": endpoint,
            "member_count": len(weights[endpoint]),
            "mean_absolute_change": float(np.abs(new - old).mean()),
            "maximum_absolute_change": float(np.abs(new - old).max()),
            "base_spearman": float(spearmanr(old, new).statistic),
        }
    )
pl.DataFrame(comparison).write_csv(out / "candidate-comparison.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "recipes": cfg["recipes"],
            "weights": weights,
            "comparison_to_base": comparison,
            "placement": placement,
            "all_four_endpoints_changed": True,
            "test_labels_used": False,
            "ensemble_topology": "target-specific flat atomic members",
        },
        indent=2,
    )
    + "\n"
)
print(json.dumps(comparison), flush=True)
