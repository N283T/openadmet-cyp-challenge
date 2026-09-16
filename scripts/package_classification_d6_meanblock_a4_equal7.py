"""Package the selected D6 PCA/raw mean-block and A4 equal-seven candidates."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata

from cyp_challenge.tracks import ENDPOINTS, validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
D6, A4 = "CYP2D6_is_TDI", "CYP3A4_is_TDI"


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def top_calls(values: np.ndarray, count: int) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    result[np.argsort(-values, kind="stable")[:count]] = True
    return result


def model_score(model: object, features: np.ndarray) -> np.ndarray:
    if hasattr(model, "predict_proba"):
        return model.predict_proba(features)[:, 1]
    return model.decision_function(features)


test_mapping = read("test_mapping").select("compound_id").unique().sort("compound_id")
test_ids = test_mapping.select("compound_id")
if test_ids.height != 750:
    raise ValueError(f"Expected 750 test compounds, found {test_ids.height}")


def aligned_fold_mean(source: str, method: str, column: str) -> np.ndarray:
    table = (
        read(source)
        .filter((pl.col("endpoint") == D6) & (pl.col("method") == method))
        .group_by("compound_id")
        .agg(pl.col(column).mean().alias("prediction"))
    )
    joined = test_ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    if not joined["prediction"].is_finite().all():
        raise ValueError(f"Incomplete D6 predictions for {source}/{method}")
    return joined["prediction"].to_numpy()


# Exact full-data seven-block selection from the controlled OOF run. The external-mixed
# block is the equal mean of its PCA256 and raw600 predictions, so it represents two
# fitted models while retaining the selected block's original weight.
core_members = [
    aligned_fold_mean("lgbm_test", "morgan_plain", "prediction"),
    aligned_fold_mean("lgbm_test", "mordred_plain", "prediction"),
    aligned_fold_mean("lgbm_test", "fp_mordred_plain", "prediction"),
    aligned_fold_mean("lgbm_test", "external_class_plain", "prediction"),
    aligned_fold_mean("lgbm_test", "external_mixed_plain", "prediction"),
    (
        aligned_fold_mean("pfn_test", "external_mixed_pca256", "probability")
        + aligned_fold_mean("raw600_test", "external_mixed_raw600", "probability")
    )
    / 2,
    aligned_fold_mean("graph_test", "tdi_only", "probability"),
]
d6_core = unit_rank(np.mean(core_members, axis=0))

compounds = read("compounds").sort("compound_id")
id_to_index = {compound_id: index for index, compound_id in enumerate(compounds["compound_id"])}
test_indices = np.asarray([id_to_index[compound_id] for compound_id in test_ids["compound_id"]])
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
).astype(np.float32)


def fixed_model_test(model_files: list[str]) -> np.ndarray:
    fold_scores = []
    for model_file in model_files:
        model = joblib.load(inp / model_file)
        fold_scores.append(unit_rank(model_score(model, morgan[test_indices])))
    return unit_rank(np.mean(fold_scores, axis=0))


logistic = fixed_model_test(cfg["logistic_models"])
rbf = fixed_model_test(cfg["rbf_models"])
complement = fixed_model_test(cfg["complement_models"])
d6_live = (2 * d6_core + logistic) / 3
d6_score = ((2 * d6_live + rbf) / 3 + (2 * d6_live + complement) / 3) / 2
d6_calls = top_calls(d6_score, cfg["d6_test_calls"])

a4_table = read("a4_test").filter(pl.col("model") == cfg["a4_model"]).sort("compound_id")
if not a4_table["compound_id"].equals(test_ids["compound_id"]):
    raise ValueError("A4 and D6 test identities differ")
a4_calls = a4_table["call"].to_numpy().astype(bool)
if int(a4_calls.sum()) != cfg["a4_test_calls"]:
    raise ValueError("Unexpected A4 call count")

aliases = (
    read("aliases")
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
reference = pl.read_csv(inp / cfg["reference"])
identities = reference.select("SMILES", "Molecule_Name").join(
    aliases, on="Molecule_Name", validate="1:1"
)
calls = test_ids.with_columns(
    pl.Series(D6, d6_calls),
    pl.Series(A4, a4_calls),
)
candidate = (
    identities.join(calls, on="compound_id", validate="1:1")
    .select("SMILES", "Molecule_Name", *ENDPOINTS["classification"])
    .sort("Molecule_Name")
)
candidate.write_parquet(out / "predictions.parquet")
candidate.write_csv(out / "predictions.csv")

score_table = test_ids.with_columns(
    pl.Series("d6_score", d6_score),
    pl.Series("d6_call", d6_calls),
    pl.Series("a4_prediction", a4_table["prediction"]),
    pl.Series("a4_call", a4_calls),
)
score_table.write_parquet(out / "test-scores.parquet")
report = validate(out / "predictions.parquet", "classification", inp / cfg["reference"])
summary = {
    "D6_model": (
        "full::mean_block from controlled D6 raw600 replacement; seven core blocks "
        "plus logistic/RBF/Complement tails"
    ),
    "D6_OOF": {"AUROC": 0.58076, "AP": 0.29686, "MCC_target_prior": 0.11011},
    "D6_positive_calls": int(d6_calls.sum()),
    "A4_model": cfg["a4_model"],
    "A4_positive_calls": int(a4_calls.sum()),
    "validation": report,
    "uploaded": False,
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
