"""Verify persisted model inference and create an exact two-target reviewable payload."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from torch import nn

from cyp_challenge.tracks import validate


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    torch.set_num_threads(4)
    policy = pl.read_parquet(inp / "blind-policy.parquet")
    folds = pl.read_parquet(inp / "blind-folds.parquet")
    features = pl.read_parquet(inp / "features.parquet")
    columns = [c for c in features.columns if c.startswith("emb_")]
    embedded = (
        policy.select("compound_id")
        .join(features, on="compound_id", how="left", validate="1:1", maintain_order="left")
        .select(columns)
        .to_numpy()
        .astype(np.float32)
    )
    assert np.isfinite(embedded).all()
    scores, margins, manifest = [], [], []
    for filename in cfg["model_files"]:
        path = inp / filename
        saved = torch.load(path, map_location="cpu", weights_only=False)
        network = nn.Sequential(
            nn.Linear(2048, 256),
            nn.ReLU(),
            nn.Dropout(0.15),
            nn.Linear(256, 64),
            nn.ReLU(),
            nn.Linear(64, saved["architecture"][-1]),
        )
        network.load_state_dict(saved["network_state"])
        network.eval()
        normalized = ((embedded - saved["feature_mean"]) / saved["feature_scale"]).astype(
            np.float32
        )
        with torch.no_grad():
            score = torch.sigmoid(network(torch.as_tensor(normalized))[:, 0]).numpy()
        recorded = folds.filter(
            (pl.col("seed") == saved["seed"]) & (pl.col("fold") == saved["fold"])
        )
        recorded = policy.select("compound_id").join(
            recorded, on="compound_id", how="left", validate="1:1", maintain_order="left"
        )
        np.testing.assert_allclose(score, recorded["score"].to_numpy(), atol=1e-7, rtol=0)
        assert np.isfinite(saved["cut"])
        scores.append(score)
        margins.append(score - saved["cut"])
        manifest.append(
            dict(
                file=filename,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                seed=saved["seed"],
                fold=saved["fold"],
                cut=saved["cut"],
            )
        )
    assert len(scores) == 10
    averaged = np.mean(margins, axis=0)
    np.testing.assert_allclose(averaged, policy["margin"].to_numpy(), atol=1e-7, rtol=0)
    np.testing.assert_array_equal(averaged > 0, policy["call"].to_numpy())
    last = pl.read_parquet(inp / "last-sent.parquet")
    best = pl.read_parquet(inp / "best-control.parquet")
    keyed = last.join(
        policy.rename({"molecule_name": "Molecule_Name"}),
        on="Molecule_Name",
        how="left",
        validate="1:1",
        maintain_order="left",
    )
    assert keyed["margin"].is_finite().all() and keyed.height == 750
    payload = keyed.select(last.columns).with_columns(pl.Series("CYP2D6_is_TDI", keyed["call"]))
    assert payload["CYP3A4_is_TDI"].equals(last["CYP3A4_is_TDI"])
    payload.write_parquet(out / "predictions.parquet")
    payload.write_csv(out / "predictions.csv")
    keyed.write_parquet(out / "blind-review.parquet")
    validation = validate(out / "predictions.parquet", "classification", inp / "reference.csv")
    changes = []
    for name, reference in [("last_sent", last), ("historical_best", best)]:
        reference = payload.select("Molecule_Name").join(
            reference, on="Molecule_Name", validate="1:1", maintain_order="left"
        )
        for endpoint in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
            a, b = reference[endpoint].to_numpy(), payload[endpoint].to_numpy()
            changes.append(
                dict(
                    reference=name,
                    endpoint=endpoint,
                    positives=int(b.sum()),
                    added=int((b & ~a).sum()),
                    removed=int((a & ~b).sum()),
                    retained=int((a & b).sum()),
                )
            )
    pl.DataFrame(changes).write_csv(out / "target-changes.csv")
    (out / "model-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                validation=validation,
                target_changes=changes,
                model_roundtrip_verified=True,
                models=manifest,
                no_upload=True,
                decision="Fresh CheMeleon plus four-primary-assay 20-epoch readout, "
                "ten equal-margin bag models. "
                "A4 retains the latest exact leaderboard-best payload. "
                "No D6 incumbent score or positive budget.",
                limits="Not a demonstrated Live improvement "
                "or identification of historical hidden labels. "
                "Encoder original pretraining corpus unknown; source provenance retained.",
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
