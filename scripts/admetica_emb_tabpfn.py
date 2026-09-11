"""Frozen matched-CYP Admetica embeddings alone into TabPFN 2.6."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers, nn
from chemprop.data import BatchMolGraph
from rdkit import Chem
from tabpfn import TabPFNRegressor
from tabpfn.constants import ModelVersion

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name):
    return pl.read_parquet(next(inp.glob("*-" + name + ".parquet")))


compounds = read("compounds").sort("compound_id")
labels = read("measurements")
a = read("parent-assignments")
torch.set_num_threads(4)
torch.backends.cuda.matmul.allow_tf32 = False
f = featurizers.SimpleMoleculeMolGraphFeaturizer()
graphs = [f(Chem.MolFromSmiles(s)) for s in compounds["standardized_smiles"]]
rows, audits = [], []
for path in sorted(inp.glob("*-encoder.pt")):
    enzyme = next(e for e in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"] if e.lower() in path.name)
    ep = enzyme + "_pIC50_direct_inhibition"
    with torch.serialization.safe_globals([nn.Activation]):
        ckpt = torch.load(path, map_location="cpu", weights_only=True)
    mp = nn.BondMessagePassing(**ckpt["hyper_parameters"]).cuda().eval().requires_grad_(False)
    mp.load_state_dict(ckpt["state_dict"], strict=True)
    blocks = []
    with torch.inference_mode():
        for start in range(0, len(graphs), 128):
            graph = BatchMolGraph(graphs[start : start + 128])
            graph.to("cuda")
            blocks.append(nn.MeanAggregation()(mp(graph), graph.batch).cpu().numpy())
    values = np.concatenate(blocks)
    assert values.shape == (6896, 300) and np.isfinite(values).all()
    cols = [f"emb_{i:03d}" for i in range(300)]
    emb = pl.DataFrame(values, schema=cols).insert_column(0, compounds["compound_id"])
    emb.write_parquet(out / f"{enzyme}-embeddings.parquet")
    h = (
        labels.filter(pl.col("endpoint") == ep)
        .join(a, on="compound_id", validate="1:1")
        .join(emb, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    x = h.select(cols).to_numpy()
    y = h["value"].to_numpy()
    for fold in range(5):
        va = (h["fold"].to_numpy() == fold) & ~h["is_anchor"].to_numpy()
        tr = ~va
        model = TabPFNRegressor.create_default_for_version(
            ModelVersion.V2_6,
            n_estimators=8,
            model_path=str(inp / cfg["weights"]),
            device="cuda",
            random_state=42 + fold,
            auto_scale_n_estimators=False,
            n_preprocessing_jobs=1,
            memory_saving_mode="auto",
            show_progress_bar=False,
            ignore_pretraining_limits=True,
            fit_mode="fit_preprocessors",
        )
        model.fit(x[tr], y[tr])
        pred = model.predict(x[va]).astype(float)
        assert np.isfinite(pred).all()
        model.save_fit_state(out / f"{ep}-fold{fold}.tabpfn_fit")
        rows.append(
            h.filter(pl.Series(va))
            .select(
                "compound_id",
                "endpoint",
                "fold",
                "is_neighbor",
                "conf_low",
                "conf_high",
                pl.col("value").alias("y_true"),
            )
            .with_columns(
                pl.lit("admetica_emb_tabpfn").alias("method"),
                pl.lit(cfg["split_run"]).alias("split_id"),
                pl.Series("prediction", pred),
            )
        )
        audits.append(
            dict(
                endpoint=ep,
                fold=fold,
                training=int(tr.sum()),
                validation=int(va.sum()),
                features=300,
            )
        )
        print(json.dumps(audits[-1]), flush=True)
pl.concat(rows).write_parquet(out / "admetica-emb-oof.parquet")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            fits=audits,
            features="mean pooled frozen message-passing output before BatchNorm, no other columns",
            encoder_finetuning=False,
        ),
        indent=2,
    )
)
