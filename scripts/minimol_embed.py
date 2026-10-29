"""Frozen MiniMol 512-d fingerprints for every compound, on CPU.

Upstream `Minimol.__call__` returns fewer rows than inputs when featurization fails,
so each batch is checked for a one-to-one result and retried molecule by molecule
otherwise; unembeddable compounds are reported instead of shifting later rows.
"""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from minimol import Minimol


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    torch.set_num_threads(cfg["threads"])
    compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
    assert compounds["compound_id"].n_unique() == compounds.height == cfg["expected_compounds"]
    ids, smiles = compounds["compound_id"].to_list(), compounds["standardized_smiles"].to_list()
    model = Minimol(batch_size=cfg["batch_size"])

    def embed(chunk):
        try:
            result = model(chunk)
        except Exception:
            return None
        return result if len(result) == len(chunk) else None

    vectors, failed = {}, []
    size = cfg["batch_size"]
    for start in range(0, len(ids), size):
        chunk_ids, chunk = ids[start : start + size], smiles[start : start + size]
        result = embed(chunk)
        if result is None:
            result = []
            for cid, smi in zip(chunk_ids, chunk, strict=True):
                single = embed([smi])
                if single is None:
                    failed.append(cid)
                result.append(None if single is None else single[0])
        for cid, vector in zip(chunk_ids, result, strict=True):
            if vector is not None:
                vectors[cid] = vector.detach().float().numpy()
        print(
            f"{min(start + size, len(ids))}/{len(ids)} embedded, {len(failed)} failed", flush=True
        )

    kept = sorted(vectors)
    matrix = np.stack([vectors[cid] for cid in kept]).astype(np.float32)
    assert np.isfinite(matrix).all()
    names = [f"emb_{i:03d}" for i in range(matrix.shape[1])]
    pl.DataFrame({"compound_id": kept}).hstack(pl.DataFrame(matrix, schema=names)).write_parquet(
        out / "external-embeddings.parquet"
    )
    summary = {
        "compounds": len(ids),
        "embedded": len(kept),
        "failed": failed,
        "dimension": matrix.shape[1],
        "torch": torch.__version__,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
