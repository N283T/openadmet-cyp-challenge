"""Global, endpoint-specific Caruana on a frozen catalog of saved OOF predictions."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(i):
    return pl.read_parquet(next(inp.glob(f"{i:03d}-*")))


keys = ["compound_id", "endpoint", "fold"]
cols = keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]
ref = read(0).select(cols + ["prediction"]).sort(keys)
assert ref.height == 6384 and not ref.select(keys).is_duplicated().any()
endpoints = sorted(ref["endpoint"].unique())
raw, rejected = [], []
for i, source in enumerate(cfg["sources"], 2):
    t = read(i)
    if "y_true" not in t.columns and "value" in t.columns:
        t = t.rename({"value": "y_true"})
    if not all(c in t.columns for c in keys + ["y_true", "prediction"]):
        rejected.append(dict(source=source, reason="unsupported schema", columns=t.columns))
        continue
    discriminators = [c for c in ["method", "model", "variant"] if c in t.columns]
    groups = t.partition_by(discriminators, as_dict=True) if discriminators else {(): t}
    for label, h in groups.items():
        label = "|".join(str(x) for x in label)
        raw.append(dict(source=source, label=label, table=h))
# The early component experiment wrote one file per fold. Assemble each named
# component from the latest completed source per fold, without looking at scores.
components = [r for r in raw if r["source"]["script"] == "ensemble_components.py"]
for label in sorted({r["label"] for r in components}):
    pieces, provenance = [], []
    for fold in range(5):
        choices = [r for r in components if r["label"] == label and fold in r["table"]["fold"]]
        if choices:
            selected = sorted(choices, key=lambda r: r["source"]["run"])[-1]
            pieces.append(selected["table"].filter(pl.col("fold") == fold))
            provenance.append(selected["source"])
    if len(pieces) == 5:
        raw.append(
            dict(
                source=dict(
                    run="assembled_components",
                    script="ensemble_components.py",
                    purpose=label,
                    inputs=provenance,
                ),
                label=label,
                table=pl.concat(pieces, how="diagonal_relaxed"),
            )
        )

vectors, catalog, aliases = {}, {}, []
for ep in endpoints:
    base = ref.filter(pl.col("endpoint") == ep)
    seen = {}
    for r in raw:
        h = r["table"].filter(pl.col("endpoint") == ep).sort(keys)
        reason = None
        if h.height != base.height or h.select(keys).is_duplicated().any():
            reason = "incomplete endpoint coverage or duplicate keys"
        elif any(not h[c].equals(base[c]) for c in keys + ["y_true"]):
            reason = "key/fold/label mismatch"
        elif any(
            c in h.columns and not h[c].equals(base[c])
            for c in ["conf_low", "conf_high", "is_neighbor"]
        ):
            reason = "confidence bounds or neighbor mismatch"
        if reason:
            rejected.append(dict(source=r["source"], label=r["label"], endpoint=ep, reason=reason))
            continue
        pred = h["prediction"].to_numpy().astype(float)
        if not np.isfinite(pred).all():
            rejected.append(
                dict(source=r["source"], label=r["label"], endpoint=ep, reason="nonfinite")
            )
            continue
        digest = hashlib.sha256(pred.tobytes()).hexdigest()
        identity = json.dumps([r["source"]["run"], r["label"], r["source"].get("path")])
        if digest in seen:
            aliases.append(
                dict(endpoint=ep, retained=seen[digest], source=r["source"], label=r["label"])
            )
            continue
        if identity not in catalog:
            catalog[identity] = dict(
                column=f"m{len(catalog):04d}", source=r["source"], label=r["label"], endpoints=[]
            )
        entry = catalog[identity]
        entry["endpoints"].append(ep)
        name = entry["column"]
        seen[digest] = name
        if name not in vectors:
            vectors[name] = np.full(ref.height, np.nan)
        vectors[name][ref["endpoint"].to_numpy() == ep] = pred
assert len(vectors) > 10, len(vectors)
aligned = ref.select(cols).with_columns([pl.Series(k, v) for k, v in vectors.items()])
aligned.write_parquet(out / "aligned-pool.parquet")
(out / "catalog.json").write_text(json.dumps(list(catalog.values()), indent=2))
(out / "excluded.json").write_text(json.dumps(rejected, indent=2))
(out / "aliases.json").write_text(json.dumps(aliases, indent=2))


def greedy(x, lo, hi):
    total = np.zeros(len(lo))
    counts = np.zeros(x.shape[1], dtype=int)
    best_loss, best_counts, history = np.inf, None, []
    for step in range(1, cfg["steps"] + 1):
        pred = (total[:, None] + x) / step
        losses = (np.maximum(lo[:, None] - pred, 0) + np.maximum(pred - hi[:, None], 0)).mean(
            axis=0
        )
        choice = int(np.argmin(losses))
        counts[choice] += 1
        total += x[:, choice]
        score = float(losses[choice])
        history.append(dict(step=step, column_index=choice, training_interval_mae=score))
        if score < best_loss - 1e-14:
            best_loss, best_counts = score, counts.copy()
    weights = best_counts / best_counts.sum()
    return weights, history


rows = [ref.with_columns(pl.lit("latest").alias("method"))]
weight_records, histories, metrics = [], [], []
for ep in endpoints:
    h = aligned.filter(pl.col("endpoint") == ep)
    names = [n for n in vectors if np.isfinite(h[n].to_numpy()).all()]
    assert len(names) > 3
    x = h.select(names).to_numpy()
    lo, hi, folds = (h[c].to_numpy() for c in ["conf_low", "conf_high", "fold"])
    print(ep, "eligible_candidates=", len(names), flush=True)
    for method in ["caruana_full_fit", "caruana_weight_crossfit"]:
        pred = np.full(h.height, np.nan)
        for held in [None] if method == "caruana_full_fit" else range(5):
            train = np.ones(h.height, dtype=bool) if held is None else folds != held
            apply = np.ones(h.height, dtype=bool) if held is None else folds == held
            w, history = greedy(x[train], lo[train], hi[train])
            pred[apply] = x[apply] @ w
            weight_records.append(
                dict(
                    endpoint=ep,
                    method=method,
                    held_fold=held,
                    selection_fold_ids=sorted(set(folds[train].tolist())),
                    weights={n: float(v) for n, v in zip(names, w, strict=True) if v > 0},
                )
            )
            histories.append(
                dict(endpoint=ep, method=method, held_fold=held, names=names, history=history)
            )
        assert np.isfinite(pred).all()
        rows.append(
            h.select(cols).with_columns(
                pl.Series("prediction", pred), pl.lit(method).alias("method")
            )
        )
    # Pairwise prediction correlations remain inspectable by endpoint.
    correlations = np.corrcoef(x, rowvar=False)
    pl.DataFrame(correlations, schema=names).with_columns(pl.Series("candidate", names)).write_csv(
        out / f"{ep}-correlations.csv"
    )
pl.concat(rows, how="diagonal_relaxed").sort(["method"] + keys).write_parquet(
    out / "global-oof.parquet"
)
(out / "weights.json").write_text(json.dumps(weight_records, indent=2))
(out / "selection-history.json").write_text(json.dumps(histories, indent=2))
(out / "summary.json").write_text(
    json.dumps(
        dict(
            sources=len(cfg["sources"]),
            candidates=len(catalog),
            endpoints={ep: sum(ep in c["endpoints"] for c in catalog.values()) for ep in endpoints},
            policy=cfg["policy"],
            crossfit_interpretation=cfg["crossfit_interpretation"],
        ),
        indent=2,
    )
)
