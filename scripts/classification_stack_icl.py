"""Change only the classifier on frozen historical nested-stack feature matrices."""

import gc
import json
import os
import time
from pathlib import Path

import numpy as np
import polars as pl
import torch
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from tabicl import TabICLClassifier

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
torch.set_num_threads(4)
test = pl.read_parquet(inp / cfg["test"])
method = "historical_stack_icl"
oofs = []
inners = []
tests = []
fits = []
cuts = []
for fold in range(5):
    for ep in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
        pair = {}
        for phase in ["inner", "outer"]:
            item = next(
                i
                for i in cfg["matrices"]
                if i["outer"] == fold and i["endpoint"] == ep and i["phase"] == phase
            )
            d = np.load(inp / item["features"])
            members = pl.read_parquet(inp / item["membership"])
            train = members.filter(pl.col("training"))
            val = members.filter(pl.col("validation"))
            np.testing.assert_array_equal(train["y_true"], d["y_train"])
            np.testing.assert_array_equal(val["y_true"], d["y_validation"])
            assert not set(train["compound_id"]) & set(val["compound_id"])
            assert (train["fold"] != fold).all()
            if phase == "inner":
                assert (train["fold"] != (fold + 1) % 5).all()
            start = time.monotonic()
            model = TabICLClassifier(
                n_estimators=8,
                batch_size=1,
                device="cuda",
                random_state=42 + fold,
                n_jobs=2,
                model_path=str(inp / cfg["weights"]),
                allow_auto_download=False,
                use_amp="auto",
                use_fa3="auto",
                offload_mode="auto",
                disk_offload_dir=os.environ["TMPDIR"],
            )
            model.fit(d["train"], d["y_train"])
            positive = list(model.classes_).index(1)
            p = model.predict_proba(d["validation"])[:, positive].astype(float)
            assert np.isfinite(p).all() and np.all((p >= 0) & (p <= 1))
            t = val.select("compound_id", "endpoint", "fold", "y_true").with_columns(
                pl.Series("probability", p),
                pl.lit(method).alias("method"),
                pl.lit("20260909T135130268493Z-edba6fc0").alias("split_id"),
            )
            if phase == "inner":
                t = t.with_columns(
                    pl.lit(fold).alias("outer_fold"), pl.lit((fold + 1) % 5).alias("inner_fold")
                )
                inners.append(t)
                threshold = float(
                    max(
                        np.linspace(0.01, 0.99, 99),
                        key=lambda cut: (
                            matthews_corrcoef(t["y_true"], p >= cut),
                            -abs(cut - 0.5),
                            -cut,
                        ),
                    )
                )
                cuts.append(dict(endpoint=ep, fold=fold, method=method, threshold=threshold))
            else:
                t = t.with_columns(pl.lit(threshold).alias("threshold"))
                oofs.append(t)
                tp = model.predict_proba(d["test"])[:, positive].astype(float)
                tests.append(
                    test.with_columns(
                        pl.lit(ep).alias("endpoint"),
                        pl.lit(fold).alias("fold"),
                        pl.lit(method).alias("method"),
                        pl.lit(threshold).alias("threshold"),
                        pl.Series("probability", tp),
                    )
                )
            model.save(
                out / f"{item['name']}.tabicl_fit",
                save_model_weights=False,
                save_training_data=True,
                save_kv_cache=False,
            )
            fits.append(
                dict(
                    name=item["name"],
                    endpoint=ep,
                    fold=fold,
                    phase=phase,
                    dimensions=d["train"].shape[1],
                    n_train=len(train),
                    n_validation=len(val),
                    classes=model.classes_.tolist(),
                    positive_index=positive,
                    seconds=time.monotonic() - start,
                )
            )
            members.write_parquet(out / f"{item['name']}-membership.parquet")
            (out / "fit-history.json").write_text(json.dumps(fits, indent=2))
            if inners:
                pl.concat(inners).write_parquet(out / "inner-probabilities.parquet")
            if oofs:
                pl.concat(oofs).write_parquet(out / "oof.parquet")
            if tests:
                pl.concat(tests).write_parquet(out / "fold-test-probabilities.parquet")
            print(json.dumps(fits[-1]), flush=True)
            del model
            gc.collect()
            torch.cuda.empty_cache()
oof = pl.concat(oofs)
inner = pl.concat(inners)
testfold = pl.concat(tests)
assert len(oof) == len(inner) == 5081 and len(testfold) == 7500
rows = []
for (ep,), t in oof.group_by("endpoint"):
    y = t["y_true"].to_numpy()
    p = t["probability"].to_numpy()
    for policy, cut in [("fixed_05", 0.5), ("inner_mcc", t["threshold"].to_numpy())]:
        rows.append(
            dict(
                method=method,
                endpoint=ep,
                policy=policy,
                n=len(t),
                AUROC=float(roc_auc_score(y, p)),
                AP=float(average_precision_score(y, p)),
                MCC=float(matthews_corrcoef(y, p >= cut)),
            )
        )
pl.DataFrame(rows).write_csv(out / "metrics.csv")
pl.DataFrame(cuts).write_csv(out / "thresholds.csv")
testfold.group_by("compound_id", "endpoint", "method").agg(
    pl.col("probability").mean(), pl.col("threshold").median()
).write_parquet(out / "test-probabilities.parquet")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            fits=20,
            outer_rows=5081,
            inner_rows=5081,
            fold_test_rows=7500,
            memory=dict(batch_size=1, n_estimators=8, offload="auto"),
            feature_policy=(
                "Exact historical outer features; replay original inner PCA fit_transform and "
                "triple-excluded logfc; no added feature scaling or dimensionality reduction"
            ),
        ),
        indent=2,
    )
)
print(pl.DataFrame(rows).write_csv())
