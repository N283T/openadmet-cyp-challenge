"""Interpret Morgan-bit signal used by the retained TDI classifiers."""

import json
import os
from collections import Counter
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import fisher_exact
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def read(name):
    return pl.read_parquet(inp / cfg[name])


def bit_fragment(mol, atom_index, radius):
    if radius == 0:
        return Chem.MolFragmentToSmiles(
            mol, atomsToUse=[atom_index], rootedAtAtom=atom_index, canonical=True
        )
    bonds = list(Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_index))
    atoms = {atom_index}
    for bond_index in bonds:
        bond = mol.GetBondWithIdx(bond_index)
        atoms.add(bond.GetBeginAtomIdx())
        atoms.add(bond.GetEndAtomIdx())
    return Chem.MolFragmentToSmiles(
        mol,
        atomsToUse=sorted(atoms),
        bondsToUse=bonds,
        rootedAtAtom=atom_index,
        canonical=True,
    )


def bh_adjust(p_values):
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    adjusted = np.empty(len(values), dtype=float)
    running = 1.0
    for reverse_rank in range(len(values) - 1, -1, -1):
        index = order[reverse_rank]
        rank = reverse_rank + 1
        running = min(running, values[index] * len(values) / rank)
        adjusted[index] = running
    return adjusted


def fixed_rate_calls(values, folds, rate):
    calls = np.zeros(len(values), dtype=bool)
    for fold in range(5):
        positions = np.flatnonzero(folds == fold)
        count = int(round(rate * len(positions)))
        order = np.argsort(-values[positions], kind="stable")
        calls[positions[order[:count]]] = True
    return calls


def expected_mcc(y, binary, prevalence):
    positive = y == 1
    tpr = float(binary[positive].mean())
    fpr = float(binary[~positive].mean())
    tp = prevalence * tpr
    fn = prevalence * (1 - tpr)
    fp = (1 - prevalence) * fpr
    tn = (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator > 0 else 0.0


compounds = read("compounds").sort("compound_id")
ids = compounds.select("compound_id")
molecules = [Chem.MolFromSmiles(smiles) for smiles in compounds["standardized_smiles"]]
assert all(molecule is not None for molecule in molecules)

x = np.zeros((len(molecules), 2048), dtype=np.uint8)
x_radius = [np.zeros_like(x) for _ in range(3)]
fragment_counts = [Counter() for _ in range(2048)]
radius_counts = [Counter() for _ in range(2048)]
for row, molecule in enumerate(molecules):
    additional = rdFingerprintGenerator.AdditionalOutput()
    additional.AllocateBitInfoMap()
    fingerprint = generator.GetFingerprint(molecule, additionalOutput=additional)
    on_bits = list(fingerprint.GetOnBits())
    x[row, on_bits] = 1
    for bit, environments in additional.GetBitInfoMap().items():
        seen_radius = set()
        seen_fragment = set()
        for atom_index, radius in environments:
            seen_radius.add(radius)
            seen_fragment.add(bit_fragment(molecule, atom_index, radius))
        for radius in seen_radius:
            x_radius[radius][row, bit] = 1
            radius_counts[bit][radius] += 1
        for fragment in seen_fragment:
            fragment_counts[bit][fragment] += 1

id_to_index = {compound_id: row for row, compound_id in enumerate(ids["compound_id"])}
labels = read("labels")
assignments = read("assignments").filter(pl.col("strategy") == "umap")

attribution_rows = []
ablation_rows = []
endpoint_summaries = {}
for endpoint in cfg["endpoints"]:
    data = (
        labels.filter(pl.col("endpoint") == endpoint)
        .join(assignments.select("compound_id", "fold", "group"), on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    indices = np.asarray([id_to_index[value] for value in data["compound_id"]])
    features = x[indices].astype(np.float32)
    radius_features = [values[indices].astype(np.float32) for values in x_radius]
    y = data["value"].to_numpy().astype(int)
    folds = data["fold"].to_numpy().astype(int)

    coefficients = []
    logistic_contributions = np.zeros_like(features, dtype=np.float32)
    lgbm_contributions = np.zeros_like(features, dtype=np.float32)
    lgbm_importances = []
    ablation_predictions = {
        name: np.empty(len(y), dtype=float)
        for name in [
            "all_bits",
            "radius0_only",
            "radius1_only",
            "radius2_only",
            "radius0_1",
            "top50_train_selected",
            "top200_train_selected",
            "drop_top50_train_selected",
        ]
    }
    for fold in range(5):
        train = folds != fold
        validation = ~train
        base = LogisticRegression(
            C=cfg["logistic_c"],
            max_iter=3000,
            solver="liblinear",
            random_state=cfg["seed"] + fold,
        ).fit(features[train], y[train])
        coefficient = base.coef_[0]
        coefficients.append(coefficient)
        logistic_contributions[validation] = features[validation] * coefficient

        order = np.argsort(-np.abs(coefficient), kind="stable")
        feature_sets = {
            "all_bits": features,
            "radius0_only": radius_features[0],
            "radius1_only": radius_features[1],
            "radius2_only": radius_features[2],
            "radius0_1": np.maximum(radius_features[0], radius_features[1]),
            "top50_train_selected": features[:, order[:50]],
            "top200_train_selected": features[:, order[:200]],
            "drop_top50_train_selected": features[:, order[50:]],
        }
        for name, selected in feature_sets.items():
            model = LogisticRegression(
                C=cfg["logistic_c"],
                max_iter=3000,
                solver="liblinear",
                random_state=cfg["seed"] + fold,
            ).fit(selected[train], y[train])
            ablation_predictions[name][validation] = model.predict_proba(selected[validation])[:, 1]

        model_path = inp / cfg["lgbm_models"][endpoint][fold]
        booster = lgb.Booster(model_file=str(model_path))
        contributions = booster.predict(features[validation], pred_contrib=True)[:, :-1]
        lgbm_contributions[validation] = contributions
        gain = booster.feature_importance(importance_type="gain").astype(float)
        gain /= gain.sum() if gain.sum() else 1.0
        lgbm_importances.append(gain)

    coefficients = np.stack(coefficients)
    lgbm_importances = np.stack(lgbm_importances)
    positive = y == 1
    negative = ~positive
    p_values = []
    raw_rows = []
    for bit in range(2048):
        pos_on = int(features[positive, bit].sum())
        neg_on = int(features[negative, bit].sum())
        pos_off = int(positive.sum() - pos_on)
        neg_off = int(negative.sum() - neg_on)
        odds_ratio, p_value = fisher_exact([[pos_on, pos_off], [neg_on, neg_off]])
        delta = pos_on / positive.sum() - neg_on / negative.sum()
        fold_deltas = []
        for fold in range(5):
            heldout = folds == fold
            heldout_positive = heldout & positive
            heldout_negative = heldout & negative
            fold_deltas.append(
                float(
                    features[heldout_positive, bit].mean()
                    - features[heldout_negative, bit].mean()
                )
            )
        fragments = fragment_counts[bit]
        radii = radius_counts[bit]
        leading_fragments = [item[0] for item in fragments.most_common(5)]
        raw_rows.append(
            {
                "endpoint": endpoint,
                "bit": bit,
                "molecules_with_bit": int(features[:, bit].sum()),
                "positive_with_bit": pos_on,
                "negative_with_bit": neg_on,
                "positive_rate_with_bit": pos_on / max(pos_on + neg_on, 1),
                "positive_prevalence": float(positive.mean()),
                "prevalence_delta": delta,
                "odds_ratio": float(odds_ratio),
                "fisher_p": float(p_value),
                "enrichment_fold_sign_consistency": float(
                    np.mean(np.sign(fold_deltas) == np.sign(delta))
                ),
                "enrichment_fold_deltas": json.dumps(fold_deltas),
                "logistic_coefficient_mean": float(coefficients[:, bit].mean()),
                "logistic_coefficient_sd": float(coefficients[:, bit].std()),
                "logistic_sign_consistency": float(
                    max(np.mean(coefficients[:, bit] > 0), np.mean(coefficients[:, bit] < 0))
                ),
                "logistic_oof_mean_abs_contribution": float(
                    np.abs(logistic_contributions[:, bit]).mean()
                ),
                "lgbm_gain_mean": float(lgbm_importances[:, bit].mean()),
                "lgbm_gain_fold_count": int(np.sum(lgbm_importances[:, bit] > 0)),
                "lgbm_oof_mean_abs_contribution": float(
                    np.abs(lgbm_contributions[:, bit]).mean()
                ),
                "lgbm_oof_mean_contribution_when_present": float(
                    lgbm_contributions[features[:, bit] > 0, bit].mean()
                )
                if features[:, bit].sum()
                else 0.0,
                "radius0_occurrences": int(radii[0]),
                "radius1_occurrences": int(radii[1]),
                "radius2_occurrences": int(radii[2]),
                "distinct_fragment_count": len(fragments),
                "top_fragments": " | ".join(leading_fragments),
            }
        )
        p_values.append(p_value)
    adjusted = bh_adjust(p_values)
    for row, q_value in zip(raw_rows, adjusted, strict=True):
        row["fisher_q"] = float(q_value)
        attribution_rows.append(row)

    for name, prediction in ablation_predictions.items():
        calls = fixed_rate_calls(prediction, folds, cfg["call_rates"][endpoint])
        ablation_rows.append(
            {
                "endpoint": endpoint,
                "feature_set": name,
                "AUROC": float(roc_auc_score(y, prediction)),
                "AP": float(average_precision_score(y, prediction)),
                "source_MCC_at_live_rate": float(matthews_corrcoef(y, calls)),
                "expected_target_MCC_at_live_rate": expected_mcc(
                    y, calls, cfg["target_priors"][endpoint]
                ),
                "calls": int(calls.sum()),
            }
        )

    endpoint_rows = pl.DataFrame(raw_rows).with_columns(pl.Series("fisher_q", adjusted))
    eligible = endpoint_rows.filter(pl.col("molecules_with_bit") >= cfg["minimum_bit_count"])
    endpoint_summaries[endpoint] = {
        "n": len(y),
        "positives": int(y.sum()),
        "bits_observed": int((features.sum(axis=0) > 0).sum()),
        "bits_fdr_0_05": int((eligible["fisher_q"] <= 0.05).sum()),
        "top_positive_logistic_bits": eligible.sort(
            "logistic_coefficient_mean", descending=True
        )["bit"][:10].to_list(),
        "top_negative_logistic_bits": eligible.sort("logistic_coefficient_mean")["bit"][
            :10
        ].to_list(),
        "top_lgbm_gain_bits": eligible.sort("lgbm_gain_mean", descending=True)["bit"][
            :10
        ].to_list(),
    }

attributions = pl.DataFrame(attribution_rows).sort(
    ["endpoint", "logistic_oof_mean_abs_contribution"], descending=[False, True]
)
attributions.write_csv(out / "bit-attributions.csv")
ablations = pl.DataFrame(ablation_rows).sort(["endpoint", "AUROC"], descending=[False, True])
ablations.write_csv(out / "ablation-metrics.csv")

# Cross-endpoint comparison is computed only on bits observed often enough in both
# endpoints; this distinguishes shared generic inhibitor chemistry from endpoint-specific signal.
comparison = (
    attributions.filter(pl.col("molecules_with_bit") >= cfg["minimum_bit_count"])
    .select(
        "endpoint",
        "bit",
        "prevalence_delta",
        "logistic_coefficient_mean",
        "lgbm_gain_mean",
        "top_fragments",
    )
    .pivot(
        on="endpoint",
        index=["bit", "top_fragments"],
        values=["prevalence_delta", "logistic_coefficient_mean", "lgbm_gain_mean"],
    )
)
comparison.write_csv(out / "cross-endpoint-bits.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "fingerprint": "Morgan radius2 2048-bit",
            "logistic_c": cfg["logistic_c"],
            "split": cfg["split_run"],
            "interpretation_scope": (
                "Exact OOF attribution for the Morgan logistic and Morgan LightGBM "
                "components; not a causal explanation of the full heterogeneous ensemble."
            ),
            "endpoints": endpoint_summaries,
        },
        indent=2,
    )
    + "\n"
)
print(ablations.write_csv(), flush=True)
print(json.dumps(endpoint_summaries, indent=2), flush=True)
