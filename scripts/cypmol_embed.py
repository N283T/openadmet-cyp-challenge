"""Frozen CYPMol (inhibitor ensemble member) molecule-by-isoform embeddings.

CYPMol fuses a Uni-Mol molecule encoder with ESM-2 residue embeddings of one CYP isoform by
cross-attention. The isoform side is identical for every molecule, so its ESM-2 embedding is
computed once. For each compound, isoform and seeded conformer this exports the 512-d vector
that feeds the final classifier and the inhibitor probability; vectors are averaged over draws.
The released checkpoint holds every weight, so no ESM-2 or Uni-Mol download is needed.
"""

import importlib.util
import io
import json
import multiprocessing as mp
import os
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import esm
import numpy as np
import polars as pl
import torch
import torch.nn as nn
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers

SEED = None


def heavy_atom_conformer(smiles):
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
    if rdForceFieldHelpers.MMFFHasAllMoleculeParams(mol):
        rdForceFieldHelpers.MMFFOptimizeMolecule(mol, maxIters=200, mmffVariant="MMFF94s")
    elif rdForceFieldHelpers.UFFHasAllMoleculeParams(mol):
        rdForceFieldHelpers.UFFOptimizeMolecule(mol, maxIters=200)
    positions = mol.GetConformer().GetPositions()
    keep = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    return [mol.GetAtomWithIdx(i).GetSymbol() for i in keep], positions[keep]


def build(item):
    compound_id, smiles = item
    return compound_id, heavy_atom_conformer(smiles)


def load_model(module, ckpt_path, dictionary):
    """EsmUnimolClassifier assembled without its hard-coded weight paths."""
    model = module.EsmUnimolClassifier.__new__(module.EsmUnimolClassifier)
    nn.Module.__init__(model)
    model.molecule_encoder = module.UniMolModel()
    alphabet = esm.data.Alphabet.from_architecture("ESM-1b")
    model.protein_encoder = esm.model.esm2.ESM2(
        num_layers=33, embed_dim=1280, attention_heads=20, alphabet=alphabet, token_dropout=True
    )
    model.alphabet = alphabet
    model.batch_converter = alphabet.get_batch_converter(truncation_seq_length=2048)
    model.transformer_layer_cross_attention = module.TransformerDecoderLayer()
    model.transformer_layer_self_attention = module.TransformerEncoderLayer()
    model.mlp = nn.Sequential(nn.Linear(512, 128), nn.ReLU(), nn.Linear(128, 2))
    state = torch.load(ckpt_path, map_location="cpu")
    model.load_state_dict(state, strict=True)
    return model.cuda().eval()


def main():
    global SEED
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    spec = importlib.util.spec_from_file_location("cypmol_model", inp / "cypmol_model.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    dict_path = str(inp / "dict.txt")
    original_load = module.Dictionary.load
    module.Dictionary.load = staticmethod(lambda _path: original_load(dict_path))
    dictionary = module.Dictionary.load("")
    dictionary.add_symbol("[MASK]", is_special=True)
    n_dict = len(dictionary)

    archive = zipfile.ZipFile(inp / "inhibitor.zip")
    training = pl.read_csv(io.BytesIO(archive.read("inhibitor.csv")), infer_schema_length=0)
    isoforms = {}
    for name in cfg["isoforms"]:
        row = training.filter(pl.col("CYP") == name).row(0, named=True)
        isoforms[name] = (row["protein_seq"], row["pocket_site"])

    compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
    if cfg.get("limit"):
        compounds = compounds.head(cfg["limit"])
    else:
        assert compounds.height == cfg["expected_compounds"]
    items = list(compounds.select("compound_id", "standardized_smiles").iter_rows())
    ids = [cid for cid, _ in items]
    conformers = []
    for seed in cfg["seeds"]:
        SEED = seed
        with ProcessPoolExecutor(cfg["workers"], mp_context=mp.get_context("fork")) as pool:
            conformers.append(dict(pool.map(build, items, chunksize=32)))
    for cid in ids:
        assert any(table[cid] is not None for table in conformers), f"no conformer for {cid}"

    sums = {n: np.zeros((len(ids), 513)) for n in isoforms}
    counts = np.zeros(len(ids))
    for member in cfg["checkpoints"]:
        model = load_model(module, inp / member, dictionary)
        protein = {}
        with torch.no_grad():
            for name, (sequence, pocket) in isoforms.items():
                _, _, tokens = model.batch_converter([("", sequence)])
                rep = model.protein_encoder(tokens.cuda(), repr_layers=[33], return_contacts=False)
                rep = rep["representations"][33]
                sites = [int(p) for p in pocket.replace("(", "").replace(")", "").split(", ")]
                protein[name] = torch.cat([rep, rep[:, [s + 1 for s in sites], :]], dim=1)
        for table in conformers:
            for index, cid in enumerate(ids):
                if table[cid] is None:
                    continue
                symbols, coords = table[cid]
                coords = coords[:256] - coords[:256].mean(0, keepdims=True)
                symbols = symbols[:256]
                tokens = (
                    [dictionary.bos()] + [dictionary.index(s) for s in symbols] + [dictionary.eos()]
                )
                tokens = torch.tensor(tokens).cuda()[None]
                xyz = torch.cat(
                    [
                        torch.zeros(1, 3),
                        torch.tensor(coords, dtype=torch.float32),
                        torch.zeros(1, 3),
                    ]
                )
                xyz = xyz.cuda()[None]
                sample = {
                    "input": {
                        "src_tokens": tokens,
                        "src_distance": torch.cdist(xyz, xyz),
                        "src_coord": xyz,
                        "src_edge_type": tokens[:, :, None] * n_dict + tokens[:, None, :],
                        "smiles": [""],
                    }
                }
                with torch.no_grad():
                    molecule = model.molecule_encoder(sample)["molecule_embedding"]
                    for name, memory in protein.items():
                        x = model.transformer_layer_cross_attention(molecule, memory, None)
                        x = model.transformer_layer_self_attention(x)[:, 0, :]
                        prob = torch.softmax(model.mlp(x), dim=1)[:, 1:2]
                        vector = torch.cat([x, prob], dim=1).float().cpu().numpy()[0]
                        sums[name][index] += vector
                counts[index] += 1
        print(f"{member}: done", flush=True)
        del model
        torch.cuda.empty_cache()
    columns = {"compound_id": ids}
    for name in isoforms:
        mean = sums[name] / counts[:, None]
        for j in range(512):
            columns[f"cypmol_{name}_{j:03d}"] = mean[:, j]
        columns[f"cypmol_{name}_prob"] = mean[:, 512]
    pl.DataFrame(columns).write_parquet(out / "external-embeddings.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {"compounds": len(ids), "isoforms": list(isoforms), "members": cfg["checkpoints"]}
        )
    )


if __name__ == "__main__":
    main()
