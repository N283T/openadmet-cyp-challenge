"""Export intermediate states of the five fold-safe Jeremy FT600 encoders."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers, nn
from chemprop.data import BatchMolGraph
from extract_chemeleon_depths import states
from rdkit import Chem


def main():
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    torch.set_num_threads(cfg["threads"])
    compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
    assert compounds.height == compounds["compound_id"].n_unique() == 6896
    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
    columns = [f"emb_{i:03d}" for i in range(cfg["width"])]
    summary = []
    for fold in range(5):
        checkpoint = torch.load(
            inp / cfg["checkpoints"][str(fold)], map_location="cpu", weights_only=False
        )
        assert checkpoint["fold"] == fold and checkpoint["split_id"] == cfg["split_id"]
        architecture = checkpoint["encoder_architecture"]
        assert architecture["depth"] == 4 and architecture["aggregation"] == "norm"
        bias = "mp.W_i.bias" in checkpoint["state_dict"]
        mp = nn.BondMessagePassing(d_h=cfg["width"], depth=architecture["depth"], bias=bias)
        mp.load_state_dict(
            {
                key.removeprefix("mp."): value
                for key, value in checkpoint["state_dict"].items()
                if key.startswith("mp.")
            },
            strict=True,
        )
        mp = mp.to(cfg["device"]).eval().requires_grad_(False)
        aggregation = nn.NormAggregation(norm=architecture["aggregation_norm"]).to(cfg["device"])
        frames = {depth: [] for depth in cfg["depths"]}
        max_forward_difference = 0.0
        with torch.inference_mode():
            for start in range(0, compounds.height, cfg["batch_size"]):
                batch = compounds.slice(start, cfg["batch_size"])
                graphs = [
                    featurizer(Chem.MolFromSmiles(smiles))
                    for smiles in batch["standardized_smiles"]
                ]
                graph = BatchMolGraph(graphs)
                graph.to(cfg["device"])
                depth_states = states(mp, graph, cfg["depths"])
                reference = aggregation(mp(graph), graph.batch).cpu().numpy()
                for depth in cfg["depths"]:
                    values = aggregation(depth_states[depth], graph.batch).cpu().numpy()
                    assert values.shape == (batch.height, cfg["width"])
                    assert np.isfinite(values).all()
                    frames[depth].append(
                        pl.DataFrame(values, schema=columns).insert_column(0, batch["compound_id"])
                    )
                    if depth == architecture["depth"]:
                        max_forward_difference = max(
                            max_forward_difference, float(np.max(np.abs(values - reference)))
                        )
        final = pl.concat(frames[4]).sort("compound_id")
        saved = pl.read_parquet(inp / cfg["reference_embeddings"][str(fold)]).sort("compound_id")
        assert final["compound_id"].equals(saved["compound_id"])
        max_saved_difference = float(
            np.max(np.abs(final.select(columns).to_numpy() - saved.select(columns).to_numpy()))
        )
        assert max_forward_difference <= cfg["validation_atol"]
        assert max_saved_difference <= cfg["validation_atol"]
        for depth in cfg["depths"]:
            pl.concat(frames[depth]).with_columns(
                pl.lit(fold).alias("fold"), pl.lit(cfg["split_id"]).alias("split_id")
            ).write_parquet(out / f"fold{fold}-depth{depth}.parquet")
        summary.append(
            {
                "fold": fold,
                "rows": final.height,
                "width": cfg["width"],
                "max_forward_difference": max_forward_difference,
                "max_saved_final_difference": max_saved_difference,
            }
        )
        print(json.dumps(summary[-1]), flush=True)
    (out / "summary.json").write_text(
        json.dumps({"config": cfg, "folds": summary}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
