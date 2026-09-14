"""Compare hashed and unhashed Morgan coverage between training and test sets."""

import json
import os
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.spatial.distance import jensenshannon
from scipy.stats import ks_2samp, spearmanr

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def read(name):
    return pl.read_parquet(inp / cfg[name])


def quantiles(values):
    return {
        "mean": float(np.mean(values)),
        "sd": float(np.std(values)),
        "p05": float(np.quantile(values, 0.05)),
        "median": float(np.median(values)),
        "p95": float(np.quantile(values, 0.95)),
    }


def flatten(prefix, values):
    return {f"{prefix}_{key}": value for key, value in quantiles(values).items()}


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


compounds = read("compounds").sort("compound_id")
labels = read("labels")
test_ids = set(read("test_mapping")["compound_id"].to_list())
molecules = [Chem.MolFromSmiles(smiles) for smiles in compounds["standardized_smiles"]]
assert all(molecule is not None for molecule in molecules)

n = len(molecules)
x = np.zeros((n, 2048), dtype=np.uint8)
x_radius = [np.zeros_like(x) for _ in range(3)]
fingerprints = []
sparse_ids = []
sparse_counts = []
fragment_counts = [Counter() for _ in range(2048)]
for row, molecule in enumerate(molecules):
    additional = rdFingerprintGenerator.AdditionalOutput()
    additional.AllocateBitInfoMap()
    fingerprint = generator.GetFingerprint(molecule, additionalOutput=additional)
    fingerprints.append(fingerprint)
    on_bits = list(fingerprint.GetOnBits())
    x[row, on_bits] = 1
    for bit, environments in additional.GetBitInfoMap().items():
        for atom_index, radius in environments:
            x_radius[radius][row, bit] = 1
            fragment_counts[bit][bit_fragment(molecule, atom_index, radius)] += 1
    elements = generator.GetSparseCountFingerprint(molecule).GetNonzeroElements()
    sparse_ids.append(set(elements))
    sparse_counts.append(elements)

compound_ids = compounds["compound_id"].to_list()
id_to_index = {compound_id: row for row, compound_id in enumerate(compound_ids)}
test_index = np.asarray([row for row, value in enumerate(compound_ids) if value in test_ids])
all_train_index = np.asarray(
    [row for row, value in enumerate(compound_ids) if value not in test_ids]
)
assert len(test_index) == 750

groups = {"all_unblinded_train": all_train_index}
for endpoint in cfg["endpoints"]:
    endpoint_ids = labels.filter(pl.col("endpoint") == endpoint)["compound_id"].to_list()
    groups[endpoint] = np.asarray([id_to_index[value] for value in endpoint_ids])

summary_rows = []
shift_rows = []
bin_rows = []
nearest_rows = []
for group_name, train_index in groups.items():
    train_x = x[train_index]
    test_x = x[test_index]
    train_density = train_x.sum(axis=1)
    test_density = test_x.sum(axis=1)
    density_ks = ks_2samp(train_density, test_density)
    train_prevalence = train_x.mean(axis=0)
    test_prevalence = test_x.mean(axis=0)
    active_train = train_prevalence > 0
    active_test = test_prevalence > 0
    novel_hashed = active_test & ~active_train
    train_distribution = train_x.sum(axis=0).astype(float) + 1e-12
    test_distribution = test_x.sum(axis=0).astype(float) + 1e-12
    train_distribution /= train_distribution.sum()
    test_distribution /= test_distribution.sum()

    train_sparse_counter = Counter()
    test_sparse_counter = Counter()
    for index in train_index:
        train_sparse_counter.update(sparse_counts[index])
    for index in test_index:
        test_sparse_counter.update(sparse_counts[index])
    train_sparse_set = set(train_sparse_counter)
    test_sparse_set = set(test_sparse_counter)
    novel_sparse = test_sparse_set - train_sparse_set

    similarities = []
    train_fingerprints = [fingerprints[index] for index in train_index]
    for test_row, index in enumerate(test_index):
        values = DataStructs.BulkTanimotoSimilarity(fingerprints[index], train_fingerprints)
        maximum = max(values)
        similarities.append(maximum)
        nearest_rows.append(
            {
                "comparison": group_name,
                "compound_id": compound_ids[index],
                "test_on_bits": int(test_density[test_row]),
                "nearest_train_tanimoto": float(maximum),
            }
        )
    similarities = np.asarray(similarities)

    train_sparse_unique = np.asarray([len(sparse_ids[index]) for index in train_index])
    test_sparse_unique = np.asarray([len(sparse_ids[index]) for index in test_index])
    train_sparse_total = np.asarray(
        [sum(sparse_counts[index].values()) for index in train_index]
    )
    test_sparse_total = np.asarray([sum(sparse_counts[index].values()) for index in test_index])
    row = {
        "comparison": group_name,
        "n_train": len(train_index),
        "n_test": len(test_index),
        **flatten("train_hashed_on_bits", train_density),
        **flatten("test_hashed_on_bits", test_density),
        "hashed_density_ks": float(density_ks.statistic),
        "hashed_density_ks_p": float(density_ks.pvalue),
        "active_train_bit_positions": int(active_train.sum()),
        "active_test_bit_positions": int(active_test.sum()),
        "shared_active_bit_positions": int((active_train & active_test).sum()),
        "test_only_bit_positions": int(novel_hashed.sum()),
        "test_occurrences_in_test_only_positions_fraction": float(
            test_x[:, novel_hashed].sum() / max(test_x.sum(), 1)
        ),
        "bit_prevalence_spearman": float(spearmanr(train_prevalence, test_prevalence).statistic),
        "bit_occurrence_js_divergence_bits": float(
            jensenshannon(train_distribution, test_distribution, base=2) ** 2
        ),
        "bits_abs_prevalence_shift_gt_0_02": int(
            (np.abs(test_prevalence - train_prevalence) > 0.02).sum()
        ),
        "bits_abs_prevalence_shift_gt_0_05": int(
            (np.abs(test_prevalence - train_prevalence) > 0.05).sum()
        ),
        "train_sparse_environment_ids": len(train_sparse_set),
        "test_sparse_environment_ids": len(test_sparse_set),
        "test_only_sparse_environment_ids": len(novel_sparse),
        "test_only_sparse_environment_ids_fraction": len(novel_sparse)
        / max(len(test_sparse_set), 1),
        "test_occurrences_from_test_only_sparse_ids_fraction": sum(
            test_sparse_counter[value] for value in novel_sparse
        )
        / max(sum(test_sparse_counter.values()), 1),
        **flatten("train_sparse_unique_ids", train_sparse_unique),
        **flatten("test_sparse_unique_ids", test_sparse_unique),
        **flatten("train_sparse_total_occurrences", train_sparse_total),
        **flatten("test_sparse_total_occurrences", test_sparse_total),
        **flatten("nearest_train_tanimoto", similarities),
        "nearest_train_tanimoto_lt_0_3_fraction": float((similarities < 0.3).mean()),
        "nearest_train_tanimoto_lt_0_4_fraction": float((similarities < 0.4).mean()),
    }
    for radius in range(3):
        row.update(
            flatten(
                f"train_radius{radius}_on_bits", x_radius[radius][train_index].sum(axis=1)
            )
        )
        row.update(
            flatten(f"test_radius{radius}_on_bits", x_radius[radius][test_index].sum(axis=1))
        )
    summary_rows.append(row)

    shifts = np.argsort(-np.abs(test_prevalence - train_prevalence), kind="stable")[:50]
    for rank, bit in enumerate(shifts, start=1):
        shift_rows.append(
            {
                "comparison": group_name,
                "rank": rank,
                "bit": int(bit),
                "train_prevalence": float(train_prevalence[bit]),
                "test_prevalence": float(test_prevalence[bit]),
                "test_minus_train": float(test_prevalence[bit] - train_prevalence[bit]),
                "train_count": int(train_x[:, bit].sum()),
                "test_count": int(test_x[:, bit].sum()),
                "top_fragments": " | ".join(
                    value for value, _ in fragment_counts[bit].most_common(5)
                ),
            }
        )

    for bin_index in range(16):
        start = bin_index * 128
        stop = start + 128
        bin_rows.append(
            {
                "comparison": group_name,
                "bit_start": start,
                "bit_stop_exclusive": stop,
                "train_active_positions": int(active_train[start:stop].sum()),
                "test_active_positions": int(active_test[start:stop].sum()),
                "train_occurrence_fraction": float(train_distribution[start:stop].sum()),
                "test_occurrence_fraction": float(test_distribution[start:stop].sum()),
            }
        )

summary = pl.DataFrame(summary_rows)
summary.write_csv(out / "summary.csv")
pl.DataFrame(shift_rows).write_csv(out / "top-bit-prevalence-shifts.csv")
pl.DataFrame(bin_rows).write_csv(out / "bit-index-bins.csv")
pl.DataFrame(nearest_rows).write_parquet(out / "test-nearest-train.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "fingerprint": "Morgan radius2 2048-bit plus unhashed sparse count IDs",
            "comparisons": list(groups),
            "diagnostic_only": True,
            "test_labels_used": False,
        },
        indent=2,
    )
    + "\n"
)
print(summary.write_csv(), flush=True)
