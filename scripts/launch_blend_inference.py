"""Bind source artifacts and launch the two promising blend members per fold."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/shared"
old = runs / "20260910T132141471698Z-343eb5ec"
base = json.loads((old / "config.json").read_text())
pretrained = [
    "20260910T111808494092Z-67df22cf",
    "20260910T111905846314Z-8e3ab14b",
    "20260910T111959350302Z-5529e407",
    "20260910T112052970926Z-dbe145e9",
    "20260910T112147200339Z-2bbfe757",
]
for fold in range(5):
    paths = []

    def add(path, paths=paths):
        path = Path(path)
        if path not in paths:
            paths.append(path)
        return f"{paths.index(path):03d}-{path.name}"

    item = base["folds"][fold]
    cfg = dict(
        fold=fold,
        quantiles=[0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7],
        aliases=add(old / "inputs" / base["aliases"]),
        weights=add(old / "inputs" / base["weights"]["regressor"]),
        foundation=[add(old / "inputs" / p) for p in base["foundation"]],
        auxiliary=[],
    )
    for src in item["auxiliary"]:
        cfg["auxiliary"].append(
            dict(
                summary=add(old / "inputs" / src["summary"]),
                predictions=[add(old / "inputs" / p) for p in src["predictions"]],
            )
        )
    for name, parent in [
        ("compact", runs / item["stack_run"]),
        ("pretrained", runs / pretrained[fold]),
    ]:
        cfg[name] = dict(run=parent.name, oof=add(parent / "outputs/oof.parquet"), endpoints={})
        if name == "pretrained":
            enc = json.loads((parent / "config.json").read_text())["encoder_run"]
            cfg[name]["embeddings"] = [
                add(p) for p in sorted((runs / enc / "outputs").glob("embeddings-*.parquet"))
            ]
        for ep in item["endpoints"]:
            if ep.endswith("is_TDI"):
                continue
            e = dict(model=add(parent / "outputs" / (ep + ".tabpfn_fit")))
            if name == "compact":
                e.update(
                    features=add(parent / "outputs" / (ep + "-features.parquet")),
                    pca=add(parent / "outputs" / (ep + "-pca.joblib")),
                )
            cfg[name]["endpoints"][ep] = e
    path = root / f"configs/calibration/blend-inference-fold{fold}-v1.json"
    with path.open("x") as f:
        f.write(json.dumps(cfg, indent=2) + "\n")
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Reload compact and pretrained regression members fold {fold}; "
        "verify OOF, predict test distributions",
        "--config",
        str(path),
    ]
    for p in paths:
        cmd += ["--input", str(p)]
    subprocess.run(cmd + ["--", "python", "scripts/infer_blend_members.py"], check=True, cwd=root)
