"""Shared molecule-level split candidates; no test structures fitted into UMAP."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.cluster import KMeans
from sklearn.model_selection import KFold
from umap import UMAP

from .storage import digest, write_json


def grouped_folds(groups, n_folds):
    unique = sorted(set(groups), key=lambda group: (-sum(g == group for g in groups), str(group)))
    if len(unique) < n_folds:
        raise ValueError("Not enough independent groups")
    counts = np.zeros(n_folds, dtype=int)
    assignment = np.full(len(groups), -1, dtype=int)
    groups = np.asarray(groups)
    for group in unique:
        mask = groups == group
        fold = int(counts.argmin())
        assignment[mask] = fold
        counts[fold] += int(mask.sum())
    return assignment


def main():
    inputs, output = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    compounds = pl.read_parquet(next(inputs.glob("*compounds.parquet"))).sort("compound_id")
    aliases = pl.read_parquet(next(inputs.glob("*aliases.parquet")))
    test_ids = set(aliases.filter(pl.col("role") == "test")["compound_id"])
    train_ids = set(aliases.filter(pl.col("role") != "test")["compound_id"].drop_nulls())
    if test_ids & train_ids:
        raise ValueError("Resolve train/test standardized overlap before splitting")
    fpconfig = config["fingerprint"]
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=fpconfig["radius"],
        fpSize=fpconfig["size"],
        includeChirality=fpconfig["include_chirality"],
    )
    mols = [Chem.MolFromSmiles(s) for s in compounds["standardized_smiles"]]
    fingerprints = np.array([generator.GetFingerprintAsNumPy(m) for m in mols], dtype=np.uint8)
    ids = compounds["compound_id"].to_list()
    train_index = [i for i, cid in enumerate(ids) if cid in train_ids]
    training_ids = [ids[i] for i in train_index]
    scaffold = [
        MurckoScaffold.MurckoScaffoldSmiles(mol=mols[i], includeChirality=False)
        for i in train_index
    ]
    np.savez_compressed(output / "fingerprints.npz", ids=np.array(ids), fingerprints=fingerprints)
    n = len(training_ids)
    assignments = {}
    random = np.full(n, -1, dtype=int)
    for fold, (_, val) in enumerate(
        KFold(config["n_folds"], shuffle=True, random_state=config["seed"]).split(training_ids)
    ):
        random[val] = fold
    assignments["random"] = (random, [str(i) for i in range(n)])
    assignments["scaffold"] = (grouped_folds(scaffold, config["n_folds"]), scaffold)
    uc = config["umap"]
    embedding = UMAP(
        n_components=uc["n_components"],
        n_neighbors=uc["n_neighbors"],
        metric=uc["metric"],
        random_state=config["seed"],
        n_jobs=1,
    ).fit_transform(fingerprints[train_index])
    clusters = KMeans(
        n_clusters=uc["n_clusters"], n_init=uc["kmeans_n_init"], random_state=config["seed"]
    ).fit_predict(embedding)
    assignments["umap"] = (
        grouped_folds(clusters.tolist(), config["n_folds"]),
        [str(c) for c in clusters],
    )
    pl.DataFrame(
        {
            "compound_id": training_ids,
            "cluster": clusters,
            **{f"umap_{i}": embedding[:, i] for i in range(embedding.shape[1])},
        }
    ).write_parquet(output / "umap_embedding.parquet")
    rows = []
    for strategy, (folds, groups) in assignments.items():
        assert (folds >= 0).all() and len(set(folds)) == config["n_folds"]
        for cid, fold, group in zip(training_ids, folds, groups, strict=True):
            rows.append(
                {"strategy": strategy, "compound_id": cid, "fold": int(fold), "group": group}
            )
    table = pl.DataFrame(rows).sort("strategy", "compound_id")
    table.write_parquet(output / "assignments.parquet")
    write_json(
        output / "summary.json",
        {
            "config": config,
            "train_aux_compounds": n,
            "test_compounds": len(test_ids),
            "fit_test_structures": False,
            "assignments_sha256": digest(output / "assignments.parquet"),
            "fold_counts": table.group_by("strategy", "fold")
            .len()
            .sort("strategy", "fold")
            .to_dicts(),
            "aux_policy": (
                "Future fold models must exclude all held-out fold compound IDs "
                "from auxiliary supervised pretraining."
            ),
        },
    )
    print(table.group_by("strategy", "fold").len().sort("strategy", "fold"), flush=True)


if __name__ == "__main__":
    main()
