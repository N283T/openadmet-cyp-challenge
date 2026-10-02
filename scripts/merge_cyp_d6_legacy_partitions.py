"""Merge complete predetermined split/seed blocks with unchanged protocol and inputs."""

import json
import os
from pathlib import Path

import polars as pl


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    common = [
        "endpoint",
        "fixed_assignments",
        "outer_folds",
        "bag_folds",
        "inner_folds",
        "bootstrap_samples",
        "cpu_threads",
        "mordred_files",
        "graph_max_epochs",
        "graph_patience",
        "model_name",
        "members",
    ]
    base, input_hashes = None, None
    predictions, thresholds, fits, metrics, leakage = [], [], [], [], []
    summaries = []
    for i, spec in enumerate(cfg["partitions"]):
        config = json.loads((inp / f"p{i}-config.json").read_text())
        result = json.loads((inp / f"p{i}-result.json").read_text())
        manifest = json.loads((inp / f"p{i}-manifest.json").read_text())
        current = {k: config[k] for k in common}
        hashes = {Path(v["snapshot"]).name: v["sha256"] for v in manifest["inputs"]}
        if base is None:
            base, input_hashes = current, hashes
        assert current == base and hashes == input_hashes
        assert result["exit_code"] == 0 or (
            spec.get("allow_completed_interrupted_block", False)
            and result["exit_code"] in [-2, 130, -15, 143]
        )

        def select(frame, spec=spec):
            return frame.filter(
                (pl.col("split_id") == spec["split_id"]) & (pl.col("seed") == spec["seed"])
            )

        pred = select(pl.read_parquet(inp / f"p{i}-nested-oof.parquet"))
        assert pred.height == 6146 and pred["compound_id"].n_unique() == 6146
        assert pred["score"].is_finite().all() and pred["margin"].is_finite().all()
        assert pred.filter((pl.col("margin") > 0) != pl.col("call")).height == 0
        assert sorted(pred["fold"].unique()) == list(range(5))
        cuts = select(pl.read_csv(inp / f"p{i}-threshold-scale-audit.csv"))
        audit = select(pl.read_csv(inp / f"p{i}-fit-audit.csv"))
        fold_metrics = select(pl.read_csv(inp / f"p{i}-fold-metrics.csv"))
        leak = select(pl.read_csv(inp / f"p{i}-leakage-and-operating-points.csv"))
        assert cuts.height == 25 and audit.height == 100 and fold_metrics.height == leak.height == 5
        assert not audit["outer_test_used"].any()
        assert not leak["outer_test_used_in_fit_or_cut"].any()
        for outer in range(5):
            assert cuts.filter(pl.col("outer_fold") == outer).height == 5
            for bag in range(5):
                block = audit.filter((pl.col("outer_fold") == outer) & (pl.col("bag_fold") == bag))
                assert block.height == 4
                inner = block.filter(pl.col("phase") == "cut_inner")
                assert sorted(inner["inner_fold"].to_list()) == [0, 1, 2]
                assert block.filter(pl.col("phase") == "bag_refit").height == 1
        for collection, frame in zip(
            [predictions, thresholds, fits, metrics, leakage],
            [pred, cuts, audit, fold_metrics, leak],
            strict=True,
        ):
            collection.append(frame)
        summaries.append(
            {
                **spec,
                "exit_code": result["exit_code"],
                "retained_complete_folds": 5,
                "retained_fits": 100,
            }
        )
    combined = pl.concat(predictions)
    assert combined.height == 4 * 6146
    assert not combined.select("compound_id", "split_id", "seed").is_duplicated().any()
    assert sorted(set(zip(combined["split_id"], combined["seed"], strict=True))) == [
        ("scaffold", 20261002),
        ("scaffold", 20261017),
        ("umap", 20261002),
        ("umap", 20261017),
    ]
    combined.write_parquet(out / "nested-oof.parquet")
    for parts, filename in [
        (thresholds, "threshold-scale-audit.csv"),
        (fits, "fit-audit.csv"),
        (metrics, "fold-metrics.csv"),
        (leakage, "leakage-and-operating-points.csv"),
    ]:
        pl.concat(parts).write_csv(out / filename)
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                common_protocol=base,
                partitions=summaries,
                fits=400,
                complete_outer_blocks=20,
                unchanged_input_hashes=input_hashes,
                no_upload=True,
                no_pretrained_assets=False,
                reused_external_assets=[
                    "class-embeddings.parquet",
                    "mixed-embeddings.parquet",
                    "external-best.pt",
                    "tabpfn-v2.6.ckpt",
                ],
                metadata_correction=(
                    "The generic source representation.json flag inspected only the "
                    "fresh embedding_file option. The legacy manifests explicitly "
                    "declare external pretrained assets. Corrected here without "
                    "altering source outputs or fitted models."
                ),
                reason="Predetermined compute blocks; no metric-dependent selection.",
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
