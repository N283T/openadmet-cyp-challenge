"""Matched CheMeleon regression/log2fc models with or without external auxiliary losses."""

import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop.data import BatchMolGraph
from rdkit import Chem
from split_baselines import metrics

from cyp_challenge.external_training import grouped_membership, snapshot_input
from cyp_challenge.log2fc_pretrain import Network, target_scaler
from cyp_challenge.neural_losses import masked_loss
from cyp_challenge.storage import digest, write_json
from cyp_challenge.tracks import ENDPOINTS


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    torch.set_num_threads(4)
    torch.manual_seed(cfg["seed"] + cfg["fold"])
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    device = torch.device("cuda")

    def read(name):
        return pl.read_parquet(snapshot_input(inp, f"{name}.parquet"))

    checkpoint_path = next(inp.glob("*-chemeleon_mp.pt"))
    assert digest(checkpoint_path) == cfg["checkpoint_sha256"]
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    compounds, labels = read("compounds"), read("measurements")
    assignments = read("assignments").filter(pl.col("strategy") == "umap")
    connections = read("challenge-connectivity")
    compound_map = compounds.join(connections, on="compound_id", validate="1:1")
    assert compound_map.height == compounds.height
    direct = ENDPOINTS["regression"]
    auxiliary = read("targets").select("compound_id", pl.col("head").alias("endpoint"), "value")
    aux_heads = sorted(auxiliary["endpoint"].unique())
    ext = read("external-targets")
    chembl_heads = sorted(ext.filter(pl.col("source") == "chembl_ic50")["endpoint"].unique())
    pubchem_heads = sorted(
        ext.filter(pl.col("source") == "pubchem_max_response")["endpoint"].unique()
    )
    heads = direct + aux_heads + chembl_heads + pubchem_heads
    assert len(heads) == 18
    target = pl.concat(
        [
            labels.filter(pl.col("endpoint").is_in(direct)).select(
                "compound_id", "endpoint", "value"
            ),
            auxiliary,
        ]
    )
    assert target.select("compound_id", "endpoint").n_unique() == target.height
    challenge = (
        target.pivot(on="endpoint", index="compound_id", values="value")
        .join(
            assignments.select("compound_id", "fold"), on="compound_id", how="left", validate="1:1"
        )
        .join(compound_map, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    assert challenge["fold"].null_count() == 0
    ids = challenge["compound_id"].to_list()
    outer_keys = assignments.filter(pl.col("fold") == cfg["fold"]).join(
        connections, on="compound_id"
    )["connectivity"]
    test_ids = compounds.join(assignments.select("compound_id"), on="compound_id", how="anti")[
        "compound_id"
    ]
    test_keys = connections.filter(pl.col("compound_id").is_in(test_ids.implode()))["connectivity"]
    train, val, held = grouped_membership(
        ids,
        challenge["connectivity"],
        challenge["fold"],
        cfg["fold"],
        cfg["seed"],
        cfg["inner_validation_fraction"],
        outer_keys,
        test_keys,
    )
    ec = read("external-compounds")
    assert set(ec["connectivity"]).isdisjoint(connections["connectivity"])
    external = (
        ext.pivot(on="endpoint", index="compound_id", values="value")
        .join(ec, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    y = challenge.select(direct + aux_heads).to_numpy().astype(np.float32)
    ey = external.select(chembl_heads + pubchem_heads).to_numpy().astype(np.float32)
    means, stds = target_scaler(y, train)
    emean, estd = target_scaler(ey, np.arange(len(ey)))
    z, ez = (y - means) / stds, (ey - emean) / estd
    role = np.full(len(ids), "excluded_connection", dtype=object)
    role[train], role[val], role[held] = "train", "inner_validation", "outer_holdout"
    challenge.select("compound_id", "connectivity", "fold").with_columns(
        pl.Series("role", role.tolist()), pl.lit(cfg["split_run"] + ":umap").alias("split_id")
    ).write_parquet(out / "membership.parquet")
    feature = featurizers.SimpleMoleculeMolGraphFeaturizer()
    graph_map = {
        cid: feature(Chem.MolFromSmiles(smi))
        for cid, smi in compounds.select("compound_id", "standardized_smiles").iter_rows()
    }
    eids = external["compound_id"].to_list()
    if cfg["external_weight"]:
        for i, (cid, smi) in enumerate(ec.select("compound_id", "standardized_smiles").iter_rows()):
            graph_map[cid] = feature(Chem.MolFromSmiles(smi))
            if i % 5000 == 0:
                print(f"External graphs {i}/{ec.height}", flush=True)

    def graph(batch_ids):
        result = BatchMolGraph([graph_map[cid] for cid in batch_ids])
        result.to(device)
        return result

    model = Network(checkpoint, cfg, len(heads)).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg["learning_rate"])
    rng, erng = np.random.default_rng(cfg["seed"]), np.random.default_rng(cfg["seed"] + 1000)
    pools = [
        np.flatnonzero(np.isfinite(ey[:, :5]).any(axis=1)),
        np.flatnonzero(np.isfinite(ey[:, 5:]).any(axis=1)),
    ]
    metadata = dict(
        config=cfg,
        heads=heads,
        means=means.tolist() + emean.tolist(),
        stds=stds.tolist() + estd.tolist(),
        checkpoint_sha256=digest(checkpoint_path),
        train=len(train),
        validation=len(val),
        heldout=len(held),
        external_compounds=ec.height,
        external_labels=ext.height,
        external_policy="All challenge connectivity excluded; equal-source random batches",
        stopping="Inner validation direct-regression standardized per-head MSE only",
        torch=torch.__version__,
        device=torch.cuda.get_device_name(),
    )
    write_json(out / "metadata.json", metadata)

    def predict(batch_ids):
        model.eval()
        predictions = []
        with torch.inference_mode():
            for start in range(0, len(batch_ids), cfg["batch_size"]):
                predictions.append(
                    model(graph(batch_ids[start : start + cfg["batch_size"]])).cpu().numpy()
                )
        return np.concatenate(predictions)

    best, best_epoch, stale = float("inf"), -1, 0
    for epoch in range(cfg["max_epochs"]):
        started = time.monotonic()
        model.train()
        losses = []
        for batch in np.array_split(
            rng.permutation(train), int(np.ceil(len(train) / cfg["batch_size"]))
        ):
            batch_ids = [ids[i] for i in batch]
            if cfg["external_weight"]:
                eb = np.concatenate(
                    [
                        erng.choice(pool, cfg["external_batch_per_source"], replace=True)
                        for pool in pools
                    ]
                )
                batch_ids += [eids[i] for i in eb]
            optimizer.zero_grad(set_to_none=True)
            pred = model(graph(batch_ids))
            truth = torch.tensor(z[batch], device=device)
            loss = masked_loss(pred[: len(batch), :4], truth[:, :4], False)
            loss = loss + cfg["log2fc_weight"] * masked_loss(
                pred[: len(batch), 4:8], truth[:, 4:8], False
            )
            if cfg["external_weight"]:
                et = torch.tensor(ez[eb], device=device)
                ep = pred[len(batch) :, 8:]
                eloss = (
                    masked_loss(ep[:, :5], et[:, :5], False)
                    + masked_loss(ep[:, 5:], et[:, 5:], False)
                ) / 2
                loss = loss + cfg["external_weight"] * eloss
            assert torch.isfinite(loss)
            loss.backward()
            optimizer.step()
            losses.append(float(loss.detach()))
        vp = predict([ids[i] for i in val])[:, :4]
        vl = float(masked_loss(torch.tensor(vp), torch.tensor(z[val, :4]), False))
        assert np.isfinite(vl)
        if vl < best:
            best, best_epoch, stale = vl, epoch, 0
            torch.save(
                dict(state_dict=model.state_dict(), metadata=metadata, epoch=epoch, val_loss=vl),
                out / "best.pt",
            )
        else:
            stale += 1
        entry = dict(
            epoch=epoch,
            train_loss=float(np.mean(losses)),
            validation_mse=vl,
            best_epoch=best_epoch,
            seconds=time.monotonic() - started,
        )
        with (out / "history.jsonl").open("a") as stream:
            stream.write(json.dumps(entry) + "\n")
        print(json.dumps(entry), flush=True)
        if stale >= cfg["patience"]:
            break
    torch.save(
        dict(
            state_dict=model.state_dict(),
            optimizer=optimizer.state_dict(),
            epoch=epoch,
            metadata=metadata,
        ),
        out / "last.pt",
    )
    saved = torch.load(out / "best.pt", map_location="cpu", weights_only=False)
    model.load_state_dict(saved["state_dict"])
    hp = predict([ids[i] for i in held])[:, :4] * stds[:4] + means[:4]
    rows, scores = [], []
    for j, ep in enumerate(direct):
        selected = np.isfinite(y[held, j])
        f = labels.filter(pl.col("endpoint") == ep).join(
            pl.DataFrame(
                {"compound_id": [ids[i] for i in held[selected]], "prediction": hp[selected, j]}
            ),
            on="compound_id",
            validate="1:1",
        )
        scores.append(
            dict(
                endpoint=ep,
                fold=cfg["fold"],
                n=f.height,
                **metrics(
                    f["value"].to_numpy(),
                    f["prediction"].to_numpy(),
                    f["conf_low"].to_numpy(),
                    f["conf_high"].to_numpy(),
                    False,
                ),
            )
        )
        rows.append(
            f.select("compound_id", pl.col("value").alias("y_true"), "prediction").with_columns(
                pl.lit(ep).alias("endpoint"),
                pl.lit(cfg["fold"]).alias("fold"),
                pl.lit(cfg["split_run"] + ":umap").alias("split_id"),
            )
        )
    pl.concat(rows).write_parquet(out / "oof.parquet")
    test_ids = sorted(test_ids)
    tp = predict(test_ids)[:, :4] * stds[:4] + means[:4]
    pl.concat(
        [
            pl.DataFrame({"compound_id": test_ids, "prediction": tp[:, j]}).with_columns(
                pl.lit(ep).alias("endpoint"),
                pl.lit(cfg["fold"]).alias("fold"),
                pl.lit(cfg["split_run"] + ":umap").alias("split_id"),
            )
            for j, ep in enumerate(direct)
        ]
    ).write_parquet(out / "test-predictions.parquet")
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            metrics=scores,
            best_epoch=best_epoch,
            best_validation_mse=best,
            epochs=epoch + 1,
            metadata=metadata,
        ),
    )
    print("COMPLETE", json.dumps(scores), flush=True)


if __name__ == "__main__":
    main()
