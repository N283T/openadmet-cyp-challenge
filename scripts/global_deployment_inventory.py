"""Metadata-only support ledger for the selected global126 deployment recipe."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
catalog = json.loads(next(inp.glob("000-*")).read_text())
weights = json.loads(next(inp.glob("001-*")).read_text())
records = {}
for r in cfg["records"]:
    metadata = {k: json.loads((inp / v).read_text()) for k, v in r["files"].items()}
    artifacts = [f for f in metadata["result"].get("files", []) if f["path"].startswith("outputs/")]
    models = [
        f
        for f in artifacts
        if f["path"].endswith((".txt", ".pt", ".ckpt", ".tabpfn_fit", ".joblib", ".cbm"))
    ]
    tests = [
        f
        for f in artifacts
        if "test" in Path(f["path"]).name.lower() and f["path"].endswith((".parquet", ".csv"))
    ]
    selections = [
        f
        for f in artifacts
        if Path(f["path"]).name
        in ["selected.json", "selection.json", "weights.json", "candidates.json"]
    ]
    records[r["run"]] = dict(
        run=r["run"],
        command=metadata["manifest"]["command"],
        config=metadata["config"],
        inputs=metadata["manifest"]["inputs"],
        model_artifacts=models,
        test_artifacts=tests,
        selection_artifacts=selections,
        status=metadata["result"]["status"],
        mode="Saved cross-validation assets; not a full-data refit"
        if models
        else "No saved model listed in source run; resolve upstream or refit",
    )
rows = []
for w in weights:
    if w["method"] != "caruana_full_fit" or "2D6" in w["endpoint"]:
        continue
    for mid, weight in w["weights"].items():
        c = next(c for c in catalog if c["column"] == mid)
        r = records[c["source"]["run"]]
        children = [
            x["run_id"]
            for x in r["config"].get("components", [])
            if isinstance(x, dict) and "run_id" in x
        ]
        rows.append(
            dict(
                endpoint=w["endpoint"],
                candidate=mid,
                weight=weight,
                method=c["label"],
                source_run=r["run"],
                source_script=c["source"]["script"],
                saved_model_artifacts=len(r["model_artifacts"]),
                saved_test_artifacts=len(r["test_artifacts"]),
                nested_components=";".join(children),
                status="Nested aggregate: recursively combine exact component recipes"
                if children
                else r["mode"],
            )
        )
pl.DataFrame(rows).write_csv(out / "selected-support.csv")
(out / "support-manifest.json").write_text(json.dumps(records, indent=2))
(out / "remaining-preconditions.json").write_text(
    json.dumps(
        dict(
            aggregate_test="Global126/174 have no completed aggregate test predictions",
            feature_routing="Rebuild test auxiliary features from all6896 pair predictions; "
            "reuse exact feature order and preprocessing",
            nested_aggregates="m0032/m0123 retain selected base models; m0079/m0080 need exact "
            "recursive child means, not old submission substitutes",
            baseline="m0002 upstream five-fold models are saved; see source audit chain",
            full_data="Most assets are fivefold CV fits. Averaging them is not an all-label refit; "
            "resolve deployment mode explicitly",
            calibration="Global test absent; global calibration requires new inference. "
            "Native D6 test is ready in5d73ca01",
            ready=False,
        ),
        indent=2,
    )
)
print(
    json.dumps(
        dict(
            selected_rows=len(rows),
            unique_candidates=len({r["candidate"] for r in rows}),
            source_runs=len(records),
            ready=False,
        )
    )
)
