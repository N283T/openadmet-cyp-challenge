"""Evaluate D6 decisions on a pseudo analog-expansion validation subset."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from cyp_challenge.classification_evaluation import (
    fixed_rate_calls,
    prepare_predictions,
    projected_mcc,
)

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
ids = compounds["compound_id"].to_list()
lookup = {compound_id: index for index, compound_id in enumerate(ids)}
if len(ids) != 6896:
    raise ValueError("Unexpected compound universe")
assignments = pl.read_parquet(inp / cfg["assignments"])
train_ids = set(assignments.filter(pl.col("strategy") == "umap")["compound_id"].to_list())
if len(train_ids) != 6146:
    raise ValueError("Unexpected training pool")
labels = pl.read_parquet(inp / cfg["regression_labels"])
anchors = []
for head in cfg["anchor_heads"]:
    ranked = (
        labels.filter((pl.col("head") == head) & pl.col("compound_id").is_in(train_ids))
        .sort("value", "compound_id", descending=[True, False])
        .select("compound_id")
    )
    chosen = []
    for compound_id in ranked["compound_id"]:
        if compound_id not in anchors:
            chosen.append(compound_id)
        if len(chosen) == cfg["anchors_per_head"]:
            break
    if len(chosen) != cfg["anchors_per_head"]:
        raise ValueError("Not enough distinct potent anchors")
    anchors.extend(chosen)
anchor_set = set(anchors)
generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=cfg["morgan_radius"], fpSize=cfg["morgan_bits"]
)
fps = []
for smiles in compounds["standardized_smiles"]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid standardized SMILES")
    fps.append(generator.GetFingerprint(molecule))
eligible = [compound_id for compound_id in ids if compound_id in train_ids - anchor_set]
eligible_fps = [fps[lookup[compound_id]] for compound_id in eligible]
pairs = []
for anchor_id in anchors:
    similarities = np.asarray(
        DataStructs.BulkTanimotoSimilarity(fps[lookup[anchor_id]], eligible_fps)
    )
    nearest = np.argsort(-similarities, kind="stable")[: cfg["analogs_per_anchor"]]
    pairs.extend(
        {
            "anchor_id": anchor_id,
            "analog_id": eligible[index],
            "similarity": float(similarities[index]),
        }
        for index in nearest
    )
pairs_table = pl.DataFrame(pairs)
pairs_table.write_parquet(out / "anchor-analog-pairs.parquet")
analog_ids = set(pairs_table["analog_id"].to_list())

results = []
for strategy, files in cfg["predictions"].items():
    bank = pl.concat([pl.read_parquet(inp / name) for name in files]).filter(
        pl.col("split_id") == strategy
    )
    models = cfg["models"][strategy]
    aligned = prepare_predictions(bank, assignments, models, strategy)
    parts = [aligned.filter(pl.col("model") == name).sort("compound_id") for name in models]
    if (
        parts[0].height != 6146
        or parts[0]["compound_id"].to_list() != parts[1]["compound_id"].to_list()
    ):
        raise ValueError("Model pool misalignment")
    folds = parts[0]["fold"].to_numpy()
    groups = parts[0]["group"].to_numpy()
    y = parts[0]["y_true"].to_numpy()
    selected = np.asarray([compound_id in analog_ids for compound_id in parts[0]["compound_id"]])
    scored = selected & np.isfinite(y)
    if scored.sum() < 30 or (y[scored] == 1).sum() < 5:
        raise ValueError("Pseudo analog validation has too few D6 labels")
    calls = [
        fixed_rate_calls(part["score"].to_numpy(), folds, cfg["historical_rate"]) for part in parts
    ]
    prior = cfg["target_prior"]
    truth = y[scored].astype(int)
    control = projected_mcc(truth, calls[0][scored], prior)
    ap = projected_mcc(truth, calls[1][scored], prior)
    indices = np.flatnonzero(scored)
    unique_groups = np.unique(groups[indices])
    group_rows = {g: indices[groups[indices] == g] for g in unique_groups}
    rng = np.random.default_rng(cfg["bootstrap_seed"])
    deltas = []
    for _ in range(cfg["bootstrap_replicates"]):
        draw = rng.choice(unique_groups, len(unique_groups), replace=True)
        rows = np.concatenate([group_rows[g] for g in draw])
        delta = projected_mcc(y[rows].astype(int), calls[1][rows], prior) - projected_mcc(
            y[rows].astype(int), calls[0][rows], prior
        )
        if np.isfinite(delta):
            deltas.append(delta)
    if not deltas:
        raise ValueError("No valid bootstrap samples")
    results.append(
        {
            "strategy": strategy,
            "anchor_rows": len(anchors),
            "analog_pairs": len(pairs),
            "unique_analogs": len(analog_ids),
            "scored_analogs": int(scored.sum()),
            "positive_analogs": int(truth.sum()),
            "control_calls_scored": int(calls[0][scored].sum()),
            "ap_calls_scored": int(calls[1][scored].sum()),
            "control_mcc": control,
            "ap_mcc": ap,
            "ap_minus_control": ap - control,
            "ci025": float(np.quantile(deltas, 0.025)),
            "ci975": float(np.quantile(deltas, 0.975)),
        }
    )
result = pl.DataFrame(results)
result.write_csv(out / "analog-validation.csv")
print(result, flush=True)
