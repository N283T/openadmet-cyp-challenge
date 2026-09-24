"""Apply five independent AP-mass90 outer ensembles to the blind set."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata, spearmanr

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
pool_cfg = json.loads((inp / cfg["pool_definition"]).read_text())


def read(key: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[key])


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def top_calls(values: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    calls[np.argsort(-values, kind="stable")[: cfg["test_positive_calls"]]] = True
    return calls


def fold_test(
    table: pl.DataFrame,
    ids: pl.DataFrame,
    predicate: pl.Expr,
    value_column: str = "prediction",
) -> np.ndarray:
    selected = table.filter(predicate)
    folds = sorted(selected["fold"].unique().to_list())
    if folds == [-1]:
        joined = ids.join(
            selected.select("compound_id", pl.col(value_column).alias("value")),
            on="compound_id", validate="1:1", maintain_order="left",
        )
        return unit_rank(joined["value"].to_numpy())
    values = []
    for fold in folds:
        one = selected.filter(pl.col("fold") == fold).select(
            "compound_id", pl.col(value_column).alias("value")
        )
        joined = ids.join(one, on="compound_id", validate="1:1", maintain_order="left")
        values.append(unit_rank(joined["value"].to_numpy()))
    if len(values) != 5:
        raise ValueError(f"Expected five folds, found {len(values)}")
    return unit_rank(np.mean(values, axis=0))


details = read("outer_details").filter(pl.col("method") == cfg["outer_method"]).sort("held_fold")
weights_by_fold = [json.loads(row) for row in details["weights"]]
member_names = sorted({name for weights in weights_by_fold for name in weights})
decision = read("decision_test")
control = decision.filter(pl.col("model") == cfg["control_model"]).sort("compound_id")
candidate = decision.filter(pl.col("model") == cfg["candidate_model"]).sort("compound_id")
support = read("support_test").filter(pl.col("model") == cfg["support_model"]).sort("compound_id")
ids = control.select("compound_id")
specs = {item["name"]: item for item in pool_cfg["jeremy_chemeleon_members"]}
sources = {
    "lgbm": read("lgbm_test"),
    "readout": read("readout_test"),
    "tabpfn": read("tabpfn_test"),
}
values: dict[str, np.ndarray] = {}
for name in member_names:
    if name in specs:
        spec = specs[name]
        values[name] = fold_test(
            sources[spec["source"]], ids,
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("split_id") == cfg["split_id"])
            & (pl.col("method") == spec["method"]),
        )
    elif name == "new::mace":
        values[name] = fold_test(
            read("mace_test"), ids,
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("strategy") == cfg["split_id"])
            & (pl.col("method") == "morgan_mace_off_emb_lr"),
        )
    elif name == "new::product_activity":
        values[name] = fold_test(
            read("activity_test"), ids,
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("strategy") == cfg["split_id"])
            & (pl.col("method") == "morgan_activity_endpoint_lgbm"),
        )
    elif name == "new::smartcyp":
        values[name] = fold_test(
            read("mechanism_test"), ids,
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("strategy") == cfg["split_id"])
            & (pl.col("method") == "morgan_smartcyp_lgbm"),
        )
    elif name == "pfn::external_mixed_pca256":
        values[name] = fold_test(
            read("pfn_pca_test"), ids,
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("method") == "external_mixed_pca256"), "probability",
        )
    elif name == "pfn::external_mixed_raw600":
        table = read("pfn_raw_test").filter(
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("method") == "external_mixed_raw600")
        ).select("compound_id", pl.col("probability").alias("value"))
        values[name] = unit_rank(
            ids.join(table, on="compound_id", validate="1:1", maintain_order="left")["value"]
            .to_numpy()
        )
    elif name == "complement_nb":
        compounds = read("compounds")
        mapping = read("test_mapping").select("compound_id").unique()
        ordered = ids.join(mapping, on="compound_id", validate="1:1").join(
            compounds.select("compound_id", "standardized_smiles"),
            on="compound_id", validate="1:1", maintain_order="left",
        )
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        features = np.stack([
            generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
            for smiles in ordered["standardized_smiles"]
        ]).astype(np.float32)
        ranked = [
            unit_rank(joblib.load(inp / model).predict_proba(features)[:, 1])
            for model in cfg["complement_models"]
        ]
        values[name] = unit_rank(np.mean(ranked, axis=0))
    else:
        raise ValueError(f"Missing blind source for {name}")

fold_scores = []
fold_calls = []
frames = []
for held, weights in zip(details["held_fold"], weights_by_fold, strict=True):
    score = sum(weight * values[name] for name, weight in weights.items())
    call = top_calls(score)
    fold_scores.append(score)
    fold_calls.append(call)
    frames.append(
        ids.with_columns(
            pl.lit(int(held)).alias("held_fold"),
            pl.Series("score", score), pl.Series("call", call),
        )
    )
score_matrix = np.column_stack(fold_scores)
call_matrix = np.column_stack(fold_calls)
frequency = call_matrix.sum(axis=1)
control_call = control["call"].to_numpy()
candidate_call = candidate["call"].to_numpy()
support_call = support["call"].to_numpy()
cohort = np.where(
    candidate_call & ~control_call, "candidate_only",
    np.where(control_call & ~candidate_call, "control_only", "unchanged"),
)
stability = ids.with_columns(
    pl.Series("outer_call_frequency", frequency),
    pl.Series("control_call", control_call),
    pl.Series("candidate_call", candidate_call),
    pl.Series("support6_call", support_call),
    pl.Series("cohort", cohort),
)
stability.write_parquet(out / "compound-stability.parquet")
pl.concat(frames).write_parquet(out / "fold-test.parquet")

pair_spearman = []
pair_jaccard = []
for left in range(5):
    for right in range(left + 1, 5):
        pair_spearman.append(
            float(spearmanr(score_matrix[:, left], score_matrix[:, right]).statistic)
        )
        union = fold_calls[left] | fold_calls[right]
        pair_jaccard.append(float((fold_calls[left] & fold_calls[right]).sum() / union.sum()))

def frequency_counts(mask: np.ndarray) -> dict[str, int]:
    return {str(value): int((frequency[mask] == value).sum()) for value in range(6)}

added = candidate_call & ~control_call
removed = control_call & ~candidate_call
summary = {
    "protocol": cfg["protocol"],
    "outer_member_union": member_names,
    "mean_pairwise_score_spearman": float(np.mean(pair_spearman)),
    "min_pairwise_score_spearman": float(np.min(pair_spearman)),
    "mean_pairwise_call_jaccard": float(np.mean(pair_jaccard)),
    "min_pairwise_call_jaccard": float(np.min(pair_jaccard)),
    "full_candidate_score_vs_outer_mean_spearman": float(
        spearmanr(candidate["score"], score_matrix.mean(axis=1)).statistic
    ),
    "candidate_positive_frequency": frequency_counts(candidate_call),
    "candidate_only_frequency": frequency_counts(added),
    "control_only_omission_frequency": {
        str(value): int(((5 - frequency[removed]) == value).sum()) for value in range(6)
    },
    "candidate_only_supported_4_of_5": int((frequency[added] >= 4).sum()),
    "control_only_omitted_4_of_5": int(((5 - frequency[removed]) >= 4).sum()),
    "candidate_only_overlap_support6": int((added & support_call).sum()),
    "control_only_overlap_support6_removal": int((removed & ~support_call).sum()),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
