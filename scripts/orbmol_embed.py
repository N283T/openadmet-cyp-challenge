"""Frozen OrbMol embeddings: pooled atomic features from seeded ETKDGv3+MMFF94s conformers.

Each draw embeds one explicit-hydrogen conformer per compound with a fixed ETKDG seed. The
molecule vector is the concatenation of mean- and max-pooled final-layer node features plus
summary statistics of the latent atomic charges; the exported vector is the mean over draws.
"""

import json
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl
import torch
from ase import Atoms
from orb_models.forcefield import pretrained
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

SEED = None


def conformer(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = rdDistGeom.ETKDGv3()
    params.enforceChirality = True
    params.useRandomCoords = True
    params.numThreads = 1
    params.maxIterations = 500
    params.randomSeed = SEED
    if rdDistGeom.EmbedMolecule(mol, params) < 0:
        params.ignoreSmoothingFailures = True
        params.enforceChirality = False
        params.maxIterations = 0
        if rdDistGeom.EmbedMolecule(mol, params) < 0:
            return None
    props = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol, "MMFF94s")
    if props is not None:
        rdForceFieldHelpers.MMFFOptimizeMolecule(mol, maxIters=200, mmffVariant="MMFF94s")
    elif rdForceFieldHelpers.UFFHasAllMoleculeParams(mol):
        rdForceFieldHelpers.UFFOptimizeMolecule(mol, maxIters=200)
    charge = int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))
    electrons = sum(atom.GetAtomicNum() for atom in mol.GetAtoms()) - charge
    return (
        [a.GetAtomicNum() for a in mol.GetAtoms()],
        mol.GetConformer().GetPositions(),
        charge,
        1 + electrons % 2,
    )


def build(item):
    compound_id, smiles = item
    built = conformer(smiles)
    if built is None:
        return compound_id, None
    numbers, positions, charge, spin = built
    atoms = Atoms(numbers=numbers, positions=positions)
    atoms.info["charge"] = charge
    atoms.info["spin"] = spin
    return compound_id, atoms


def main():
    global SEED
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
    assert compounds["compound_id"].n_unique() == compounds.height
    if cfg.get("limit"):
        compounds = compounds.head(cfg["limit"])
    else:
        assert compounds.height == cfg["expected_compounds"]
    items = list(compounds.select("compound_id", "standardized_smiles").iter_rows())
    model, adapter = pretrained.orbmol_v2(device=cfg["device"], precision="float32-high")
    model.eval()
    ids = [cid for cid, _ in items]
    draws = []
    for seed in cfg["seeds"]:
        SEED = seed
        with ProcessPoolExecutor(cfg["workers"], mp_context=mp.get_context("fork")) as pool:
            built = dict(pool.map(build, items, chunksize=32))
        missing = [cid for cid in ids if built[cid] is None]
        vectors = {}
        order = [cid for cid in ids if built[cid] is not None]
        for start in range(0, len(order), cfg["batch_size"]):
            chunk = order[start : start + cfg["batch_size"]]
            graph = adapter.from_ase_atoms_list([built[c] for c in chunk], device=cfg["device"])
            with torch.no_grad():
                nodes = model.model(graph)["node_features"].float()
                charges = model.heads["latent_charges"](nodes, graph).float().reshape(-1)
            offsets = torch.cumsum(torch.cat([torch.zeros(1, dtype=torch.long), graph.n_node.cpu()]), 0)
            for i, cid in enumerate(chunk):
                a, b = int(offsets[i]), int(offsets[i + 1])
                q = charges[a:b]
                stats = torch.stack([q.min(), q.max(), q.mean(), q.std(unbiased=False)])
                vectors[cid] = torch.cat([nodes[a:b].mean(0), nodes[a:b].max(0).values, stats])
                vectors[cid] = vectors[cid].cpu().numpy()
        table = pl.DataFrame(
            {"compound_id": order, **{f"emb_{j:04d}": [vectors[c][j] for c in order] for j in range(len(vectors[order[0]]))}}
        )
        table.write_parquet(out / f"draw-seed{seed}.parquet")
        draws.append(vectors)
        print(f"seed {seed}: embedded {len(order)}, missing {len(missing)}", flush=True)
    result, unembedded = {}, []
    for cid in ids:
        have = [d[cid] for d in draws if cid in d]
        if not have:
            unembedded.append(cid)
            continue
        result[cid] = np.mean(have, axis=0)
    assert not unembedded, f"no conformer for {unembedded}"
    width = len(next(iter(result.values())))
    pl.DataFrame(
        {"compound_id": ids, **{f"emb_{j:04d}": [float(result[c][j]) for c in ids] for j in range(width)}}
    ).write_parquet(out / "external-embeddings.parquet")
    (out / "summary.json").write_text(json.dumps({"compounds": len(ids), "width": width, "seeds": cfg["seeds"]}))


if __name__ == "__main__":
    main()
