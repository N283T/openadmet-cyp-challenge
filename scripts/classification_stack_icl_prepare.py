"""Exact historical outer and independently reconstructed inner feature matrices, CPU only."""

import copy
import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
torch.set_num_threads(4)


def read(name):
    return pl.read_parquet(inp / name)


labels = read(cfg["labels"])
assign = read(cfg["assignments"]).filter(pl.col("strategy") == "umap")
emb = pl.concat([read(n) for n in cfg["embedding"]])
ecols = sorted(c for c in emb.columns if c.startswith("chemeleon_"))
desc = pl.concat([read(n) for n in cfg["mordred"]])
test = read(cfg["test"]).select("compound_id").sort("compound_id")
rawtest = (
    test.join(emb, on="compound_id", validate="1:1", maintain_order="left")
    .select(ecols)
    .to_numpy()
    .astype(np.float32)
)
aux = {}
for item in cfg["auxiliary"]:
    summary = json.loads((inp / item["summary"]).read_text())
    heads = summary["heads"]
    aux[tuple(item["excluded"])] = pl.concat([read(n) for n in item["predictions"]]).select(
        "compound_id", *heads
    )
records = []
for outer in range(5):
    inner = (outer + 1) % 5
    for ep in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
        frame = (
            labels.filter(pl.col("endpoint") == ep)
            .join(assign.select("compound_id", "fold"), on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        ff = frame["fold"].to_numpy()
        for phase in ["outer", "inner"]:
            name = f"{phase}-{ep}-fold{outer}"
            item = next(
                i
                for i in cfg["readouts"]
                if i["outer"] == outer and i["endpoint"] == ep and i["phase"] == phase
            )
            if phase == "outer":
                train = ff != outer
                val = ff == outer
                saved = read(item["features"]).sort("compound_id")
                assert saved["compound_id"].equals(frame["compound_id"])
                base = np.array(saved["features"].to_list(), dtype=np.float32)
                assert base.shape == (frame.height, 260)
                pca = joblib.load(inp / item["pca"])
                imputer = joblib.load(inp / item["imputer"])
                columns = json.loads((inp / item["columns"]).read_text())
                av = np.zeros((750, 4), dtype=np.float32)
                for other in range(5):
                    if other == outer:
                        continue
                    ap = aux[tuple(sorted([outer, other]))]
                    av += (
                        test.join(ap, on="compound_id", validate="1:1", maintain_order="left")
                        .select(heads)
                        .to_numpy()
                        / 4
                    )
                testbase = np.column_stack([pca.transform(rawtest), av])
            else:
                train = ~np.isin(ff, [outer, inner])
                val = ff == inner
                prep = joblib.load(inp / item["preprocessors"])
                pca = prep["pca"]
                imputer = prep["imputer"]
                columns = prep["descriptor_columns"]
                raw = (
                    frame.select("compound_id")
                    .join(emb, on="compound_id", validate="1:1", maintain_order="left")
                    .select(ecols)
                    .to_numpy()
                    .astype(np.float32)
                )
                # Historical fit_transform uses the fitted SVD U*S on training rows.
                # Recreate that operation rather than substitute transform(X_train).
                original_pca = copy.deepcopy(pca)
                reduced = pca.transform(raw)
                reduced[train] = pca.fit_transform(raw[train])
                reduced[val] = pca.transform(raw[val])
                np.testing.assert_allclose(
                    pca.components_, original_pca.components_, rtol=0, atol=1e-7
                )
                av = np.zeros((len(frame), 4), dtype=np.float32)
                for other in range(5):
                    if other in [outer, inner]:
                        continue
                    ap = (
                        frame.select("compound_id")
                        .join(
                            aux[tuple(sorted([outer, inner, other]))],
                            on="compound_id",
                            validate="1:1",
                            maintain_order="left",
                        )
                        .select(heads)
                        .to_numpy()
                    )
                    av[ff == other] = ap[ff == other]
                    av[val] += ap[val] / 3
                base = np.column_stack([reduced, av])
                testbase = None

            def descriptors(ids, columns=columns, imputer=imputer):
                d = (
                    ids.select("compound_id")
                    .join(desc, on="compound_id", validate="1:1", maintain_order="left")
                    .select(columns)
                    .to_numpy()
                    .astype(np.float64)
                )
                d[~np.isfinite(d)] = np.nan
                return imputer.transform(d)

            x = np.column_stack([base, descriptors(frame)])
            matrices = dict(
                train=x[train],
                validation=x[val],
                y_train=frame["value"].to_numpy()[train],
                y_validation=frame["value"].to_numpy()[val],
            )
            if phase == "outer":
                matrices["test"] = np.column_stack([testbase, descriptors(test)])
            assert all(np.isfinite(v).all() for v in matrices.values())
            np.savez_compressed(out / f"{name}.npz", **matrices)
            frame.select(
                "compound_id", "endpoint", "fold", pl.col("value").alias("y_true")
            ).with_columns(
                pl.Series("training", train),
                pl.Series("validation", val),
                pl.lit(outer).alias("outer_fold"),
            ).write_parquet(out / f"{name}-membership.parquet")
            records.append(
                dict(
                    name=name,
                    endpoint=ep,
                    outer=outer,
                    inner=inner,
                    phase=phase,
                    features=f"{name}.npz",
                    membership=f"{name}-membership.parquet",
                    dimensions=x.shape[1],
                )
            )
            print(name, x.shape[1], flush=True)
test.write_parquet(out / "test-identities.parquet")
(out / "matrices.json").write_text(json.dumps(records, indent=2))
