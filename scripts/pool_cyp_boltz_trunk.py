"""Pool frozen CYP/HEM trunk tensors into reusable keyed feature blocks."""

from __future__ import annotations

import json
import os
import pickle
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl


def summarize(prefix, values, stats):
    values = values.reshape(-1, values.shape[-1]).astype(np.float32)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError(f"Invalid block: {prefix}")
    result = {}
    for stat in stats:
        if stat == "mean":
            vector = values.mean(0)
        elif stat == "max":
            vector = values.max(0)
        elif stat == "std":
            vector = values.std(0)
        else:
            vector = np.quantile(values, float(stat[1:]) / 100, axis=0)
        result.update({f"{prefix}_{stat}_{i:03d}": float(v) for i, v in enumerate(vector)})
    return result


def pool_one(item):
    path, compound_id, role, pocket, fe_index = item
    with np.load(path, allow_pickle=False) as f:
        s, z = f["s"], f["z"]
        asym, residue = f["asym_id"], f["residue_index"]
        if s.ndim != 2 or s.shape[1] != 384 or z.shape != (len(s), len(s), 128):
            raise ValueError(f"Unexpected tensor shape: {path}")
        if int(f["recycling_steps"]) != 1 or str(f["embedding_stage"]) != "post":
            raise ValueError(f"Wrong representation: {path}")
        protein, hem, ligand = [np.flatnonzero(asym == i) for i in range(3)]
        if len(hem) != 43 or not len(ligand) or len(protein) != pocket["length"]:
            raise ValueError(f"Unexpected chains: {path}")
        if not np.array_equal(residue[protein], np.arange(len(protein))):
            raise ValueError(f"Unexpected canonical numbering: {path}")
        regions = {
            "protein": protein,
            "hem": hem,
            "fe": hem[[fe_index]],
            "hem_other": np.delete(hem, fe_index),
        }
        for name in ["core", "extended"]:
            positions = pocket[f"{name}_residues"]
            regions[name] = protein[np.array(positions) - 1]
            regions[f"{name}_hem"] = np.concatenate([regions[name], hem])
        result = {"compound_id": compound_id, "role": role, "ligand_tokens": len(ligand)}
        for name in ["protein", "core", "extended", "hem"]:
            result.update(summarize(f"s_{name}", s[regions[name]], ["mean", "std", "max"]))
        result.update(summarize("s_ligand", s[ligand], ["mean", "std", "max"]))
        for name, indices in regions.items():
            for direction, a, b in [("pl", indices, ligand), ("lp", ligand, indices)]:
                result.update(
                    summarize(
                        f"z_{name}_{direction}",
                        z[a[:, None], b[None, :]],
                        ["mean", "max", "std", "q10", "q90"],
                    )
                )
        result.update(
            summarize(
                "z_ll", z[ligand[:, None], ligand[None, :]], ["mean", "max", "std", "q10", "q90"]
            )
        )
    return result


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    pocket = json.loads((inp / "pocket.json").read_text())
    pocket["length"] = cfg["protein_length"]
    # The archived CCD and parser/tokenizer source establish non-H heavy-atom order.
    with (inp / "HEM.pkl").open("rb") as stream:
        hem = pickle.load(stream)
    atoms = [a for a in hem.GetAtoms() if a.GetAtomicNum() != 1]
    fe = [i for i, a in enumerate(atoms) if a.GetSymbol() == "Fe"]
    if len(atoms) != 43 or len(fe) != 1:
        raise ValueError("Unexpected HEM reference")
    manifest = pl.read_parquet(inp / "manifest.parquet").sort("record_index")
    if manifest.select("compound_id", "role").unique().height != manifest.height:
        raise ValueError("Duplicate compound/role")
    recorded = json.loads((Path(os.environ["CYP_RUN_DIR"]) / "manifest.json").read_text())
    hashes = {Path(x["snapshot"]).name: x["sha256"] for x in recorded["inputs"]}
    jobs = []
    for r in manifest.iter_rows(named=True):
        alias = f"record_{r['record_index']:06d}.npz"
        if hashes[alias] != r["sha256"]:
            raise ValueError(f"Input hash differs from canonical manifest: {alias}")
        jobs.append((inp / alias, r["compound_id"], r["role"], pocket, fe[0]))
    (out / "batches").mkdir()
    frames, rows = [], []
    with ProcessPoolExecutor(max_workers=cfg["workers"]) as executor:
        for i, row in enumerate(executor.map(pool_one, jobs, chunksize=2), 1):
            rows.append(row)
            if len(rows) == 128 or i == len(jobs):
                frame = pl.DataFrame(rows).with_columns(pl.col(pl.Float64).cast(pl.Float32))
                p = out / "batches" / f"batch-{i:06d}.parquet"
                frame.write_parquet(p, compression="zstd")
                frames.append(p)
                rows = []
                print(f"Pooled {i}/{len(jobs)}", flush=True)
    table = pl.concat([pl.read_parquet(p) for p in frames])
    table.write_parquet(out / "features.parquet", compression="zstd")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "target": cfg["target"],
                "rows": table.height,
                "features": table.width - 2,
                "fe_heavy_atom_index": fe[0],
                "manifest_hash_checked": True,
                "coordinate_features": False,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"Completed pooling: {table.shape}", flush=True)


if __name__ == "__main__":
    main()
