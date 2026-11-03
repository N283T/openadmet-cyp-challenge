"""Cache MONROE input graphs for every compound under several seeded conformers.

Uses the same seeded featurization as `monroe_embed.py`, so cached graphs match the
frozen-embedding export; fine-tuning samples one conformer per molecule per step.
"""

import json
import multiprocessing as mp
import os
import pickle
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import monroe_embed
import polars as pl
from monroe.eval import embed as upstream
from monroe.model import featurizer


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    featurizer.predict_structure = monroe_embed.seeded_predict_structure
    upstream._graphs.tqdm.write = lambda *args, **kwargs: None
    compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
    assert compounds["compound_id"].n_unique() == compounds.height == cfg["expected_compounds"]
    items = list(compounds.select("compound_id", "standardized_smiles").iter_rows())
    summary = {}
    for seed in cfg["seeds"]:
        monroe_embed.SEED = seed
        with ProcessPoolExecutor(cfg["workers"], mp_context=mp.get_context("fork")) as pool:
            built = list(pool.map(monroe_embed.featurize, items, chunksize=8))
        graphs = {cid: graph for cid, graph, _ in built if graph is not None}
        with open(out / f"graphs-seed{seed}.pkl", "wb") as handle:
            pickle.dump(graphs, handle, protocol=pickle.HIGHEST_PROTOCOL)
        summary[seed] = {"built": len(graphs), "flat": sum(flat for *_, flat in built)}
        print(seed, summary[seed], flush=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
