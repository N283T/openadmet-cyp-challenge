"""Check whether a fixed train-to-blind density ratio changes D6 OOF transfer."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdFingerprintGenerator
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from cyp_challenge.classification_evaluation import fixed_rate_calls

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
test_ids = pl.read_parquet(inp / cfg["test_mapping"]).select("compound_id").unique()
if compounds.height != 6896 or test_ids.height != 750:
    raise ValueError("Unexpected molecular population")
domain = (
    compounds.select("compound_id")
    .join(test_ids.with_columns(pl.lit(1).alias("is_test")), on="compound_id", how="left")
    .with_columns(pl.col("is_test").fill_null(0))
    .sort("compound_id")
)
target = domain["is_test"].to_numpy().astype(int)
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
molecules = [Chem.MolFromSmiles(smiles) for smiles in compounds["standardized_smiles"]]
if any(molecule is None for molecule in molecules):
    raise ValueError("Invalid standardized SMILES")
fingerprints = np.stack([generator.GetFingerprintAsNumPy(mol) for mol in molecules]).astype(
    np.float32
)
descriptors = np.asarray(
    [
        [
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.TPSA(mol),
            Lipinski.NumHDonors(mol),
            Lipinski.NumHAcceptors(mol),
            Lipinski.RingCount(mol),
            Lipinski.NumRotatableBonds(mol),
            Lipinski.FractionCSP3(mol),
        ]
        for mol in molecules
    ],
    dtype=np.float32,
)
scaler = StandardScaler().fit(descriptors)
features = np.column_stack([fingerprints, scaler.transform(descriptors)])
probability = np.empty(len(target))
splitter = StratifiedKFold(
    n_splits=cfg["domain_folds"], shuffle=True, random_state=cfg["domain_seed"]
)
for fold, (train, held) in enumerate(splitter.split(features, target)):
    model = LogisticRegression(
        C=cfg["domain_c"], solver="liblinear", max_iter=3000, random_state=cfg["domain_seed"]
    )
    model.fit(features[train], target[train])
    probability[held] = model.predict_proba(features[held])[:, 1]
    print(f"domain fold {fold} complete", flush=True)
if not np.isfinite(probability).all():
    raise ValueError("Nonfinite domain probability")
ratio = probability / (1 - probability) * (len(target) - target.sum()) / target.sum()
ratio = np.clip(ratio, cfg["weight_min"], cfg["weight_max"])
ratio[target == 0] /= ratio[target == 0].mean()
weights = domain.select("compound_id").with_columns(
    pl.Series("is_test", target),
    pl.Series("test_probability", probability),
    pl.Series("density_weight", ratio),
)
weights.write_parquet(out / "domain-weights.parquet")


def expected_mcc(y: np.ndarray, calls: np.ndarray, weight: np.ndarray) -> float:
    positive = y == 1
    tpr = np.average(calls[positive], weights=weight[positive])
    fpr = np.average(calls[~positive], weights=weight[~positive])
    prior = cfg["target_prior"]
    tp, fn = prior * tpr, prior * (1 - tpr)
    fp, tn = (1 - prior) * fpr, (1 - prior) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


rows = []
for strategy in ("umap", "scaffold"):
    data = {}
    for role in ("control", "ap"):
        source = cfg[f"{role}_{strategy}"]
        model = cfg["control_model"] if role == "control" else cfg[f"ap_{strategy}_model"]
        frame = (
            pl.read_parquet(inp / source)
            .filter(
                (pl.col("endpoint") == cfg["endpoint"])
                & (pl.col("split_id") == strategy)
                & (pl.col("model") == model)
            )
            .sort("compound_id")
        )
        if frame.height != 6146 or frame["compound_id"].n_unique() != 6146:
            raise ValueError(f"Incomplete {strategy} {role} scores")
        data[role] = frame
    control, ap = data["control"], data["ap"]
    if not (
        control["compound_id"].equals(ap["compound_id"])
        and control["fold"].equals(ap["fold"])
        and control["y_true"].equals(ap["y_true"])
    ):
        raise ValueError(f"{strategy} prediction alignment failed")
    aligned = control.select("compound_id").join(weights, on="compound_id", validate="1:1")
    if aligned["is_test"].sum() > 0:
        raise ValueError("Blind molecule entered OOF")
    scored = control["y_true"].is_not_null().to_numpy()
    truth = control["y_true"].to_numpy()[scored].astype(int)
    w = aligned["density_weight"].to_numpy()[scored]
    folds = control["fold"].to_numpy()
    calls = {
        role: fixed_rate_calls(frame["score"].to_numpy(), folds, cfg["call_rate"])
        for role, frame in data.items()
    }
    for kind, use_weight in [("unweighted", np.ones(len(w))), ("test_density", w)]:
        control_mcc = expected_mcc(truth, calls["control"][scored], use_weight)
        ap_mcc = expected_mcc(truth, calls["ap"][scored], use_weight)
        rows.append(
            {
                "strategy": strategy,
                "weighting": kind,
                "control_mcc": control_mcc,
                "ap_mcc": ap_mcc,
                "ap_minus_control": ap_mcc - control_mcc,
                "scored_rows": int(scored.sum()),
                "scored_effective_n": float(use_weight.sum() ** 2 / (use_weight**2).sum()),
                "domain_auroc": float(roc_auc_score(target, probability)),
            }
        )
pl.DataFrame(rows).write_csv(out / "transfer.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "domain_auroc": float(roc_auc_score(target, probability)),
            "test_rows": int(target.sum()),
            "train_rows": int((target == 0).sum()),
            "ratio_clip": [cfg["weight_min"], cfg["weight_max"]],
            "interpretation": (
                "diagnostic covariate weighting; not a candidate-selection "
                "estimate or hidden-label result"
            ),
        },
        indent=2,
    )
    + "\n"
)
