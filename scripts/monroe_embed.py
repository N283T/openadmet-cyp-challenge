"""Frozen MONROE graph embeddings, averaged over seeded single-conformer draws.

Upstream `predict_structure` leaves the ETKDG seed unset, so repeated calls give
different conformers. This export pins one seed per draw and also keeps the
first-seed vector, so the conformer-noise contribution can be audited later.
"""

import json
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl
import torch
from monroe.eval import embed as monroe_embed
from monroe.model import featurizer
from monroe.model.ckpt import load_ckpt
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from torch_geometric.data import Batch

SEED = None


def seeded_predict_structure(mol, n_confs=1):
    """Upstream `predict_structure` with an explicit ETKDG random seed."""
    params = rdDistGeom.ETKDGv3()
    params.enforceChirality = True
    params.useRandomCoords = True
    params.numThreads = 1
    params.maxIterations = 500
    params.randomSeed = SEED
    mol.RemoveAllConformers()
    confs = rdDistGeom.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    if len(confs) == 0:
        raise ValueError("Failed to predict structure with ETKDGv3.")
    props = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol, "MMFF94s")
    if props is not None:
        field = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, props)
        results = rdForceFieldHelpers.OptimizeMoleculeConfs(mol, field, maxIters=200, numThreads=1)
    elif rdForceFieldHelpers.UFFHasAllMoleculeParams(mol):
        results = rdForceFieldHelpers.UFFOptimizeMoleculeConfs(mol, maxIters=200, numThreads=1)
    else:
        raise ValueError("MMFF and UFF parameters not available.")
    best = int(np.argmin([energy for _, energy in results]))
    return mol.GetConformer(int(confs[best]))


def featurize(item):
    compound_id, smiles = item
    built = monroe_embed._featurize_one(smiles)
    if built is None:
        return compound_id, None, False
    graph = built[1]
    flat = float(np.asarray(graph["pos_rdkit"])[:, 2].std()) < 1e-6
    return compound_id, graph, flat


def main():
    global SEED
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    scratch = Path(os.environ["TMPDIR"]) / "checkpoint"
    scratch.mkdir(parents=True, exist_ok=True)
    for name in ("config.json", "weights.pt"):
        link = scratch / name
        if not link.exists():
            link.symlink_to(inp / cfg["checkpoint"][name])
    featurizer.predict_structure = seeded_predict_structure
    monroe_embed._graphs.tqdm.write = lambda *args, **kwargs: None

    compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
    assert compounds["compound_id"].n_unique() == compounds.height == cfg["expected_compounds"]
    items = list(compounds.select("compound_id", "standardized_smiles").iter_rows())
    encoder = load_ckpt(str(scratch)).to(cfg["device"]).eval()

    names = None
    draws, audit = {}, []
    for seed in cfg["seeds"]:
        SEED = seed
        with ProcessPoolExecutor(cfg["workers"], mp_context=mp.get_context("fork")) as pool:
            built = list(pool.map(featurize, items, chunksize=8))
        ids = [cid for cid, graph, _ in built if graph is not None]
        graphs = [monroe_embed.to_pyg(graph) for _, graph, _ in built if graph is not None]
        vectors = []
        with torch.no_grad():
            for start in range(0, len(graphs), cfg["batch_size"]):
                batch = Batch.from_data_list(graphs[start : start + cfg["batch_size"]])
                vectors.append(encoder(batch.to(cfg["device"]))[0].float().cpu().numpy())
        matrix = np.concatenate(vectors).astype(np.float32)
        names = [f"emb_{i:03d}" for i in range(matrix.shape[1])]
        draws[seed] = dict(zip(ids, matrix, strict=True))
        # Persist each draw so a later failure cannot discard finished featurization.
        pl.DataFrame({"compound_id": ids}).hstack(pl.DataFrame(matrix, schema=names)).write_parquet(
            out / f"draw-seed{seed}.parquet"
        )
        for cid, graph, flat in built:
            audit.append(
                {"seed": seed, "compound_id": cid, "built": graph is not None, "flat": flat}
            )
        print(f"seed {seed}: {len(ids)}/{len(items)} built", flush=True)

    audit = pl.DataFrame(audit)
    audit.write_parquet(out / "featurization-audit.parquet")
    flat = {(row["seed"], row["compound_id"]) for row in audit.filter("flat").iter_rows(named=True)}

    # Average only real 3D draws; a compound with no 3D draw keeps its flat draws.
    mean_ids, mean_rows, draw_counts = [], [], []
    for cid, _ in items:
        available = [seed for seed in cfg["seeds"] if cid in draws[seed]]
        real = [seed for seed in available if (seed, cid) not in flat]
        chosen = real or available
        if not chosen:
            continue
        mean_ids.append(cid)
        mean_rows.append(np.mean([draws[seed][cid] for seed in chosen], axis=0))
        draw_counts.append(len(chosen))

    def table(ids, rows):
        matrix = np.stack(rows).astype(np.float32)
        return pl.DataFrame({"compound_id": ids}).hstack(pl.DataFrame(matrix, schema=names))

    table(mean_ids, mean_rows).write_parquet(out / "external-embeddings.parquet")
    first = draws[cfg["seeds"][0]]
    first_ids = [cid for cid, _ in items if cid in first]
    table(first_ids, [first[cid] for cid in first_ids]).write_parquet(
        out / "external-embeddings-first-seed.parquet"
    )
    dim = len(names)

    missing = sorted(set(compounds["compound_id"]) - set(mean_ids))
    summary = {
        "compounds": compounds.height,
        "embedded": len(mean_ids),
        "missing": missing,
        "dimension": dim,
        "seeds": cfg["seeds"],
        "draws_per_compound": {str(k): draw_counts.count(k) for k in sorted(set(draw_counts))},
        "flat_draws": int(audit["flat"].sum()),
        "flat_only_compounds": int(
            audit.group_by("compound_id").agg(pl.col("flat").all()).filter(pl.col("flat")).height
        ),
        "torch": torch.__version__,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
