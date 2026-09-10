"""Snapshot all saved fold models and preprocessing for baseline inference."""

import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs/shared"
READOUT = [
    "20260910T123025501770Z-492cf5f6",
    "20260910T123205649940Z-20f4f93f",
    "20260910T123336884511Z-63b24af3",
    "20260910T123510230339Z-c92321ca",
    "20260910T123640049734Z-b387c093",
]
STACK = [
    "20260910T121031603901Z-fbc2b2ea",
    "20260910T121149831493Z-a62f51ef",
    "20260910T121402665787Z-08f7afb5",
    "20260910T121548828209Z-abe37170",
    "20260910T121624900089Z-271684c2",
]


def main():
    files = []

    def add(path):
        path = Path(path)
        if path not in files:
            files.append(path)
        return f"{files.index(path):03d}-{path.name}"

    cfg = dict(
        model="TabPFN-2.6 PCA256 plus predicted log2fc4 plus Mordred",
        aggregation="mean of five existing UMAP fold models",
        threshold=0.5,
        folds=[],
    )
    cfg["reference"] = add(ROOT / "data/raw/cyp-challenge-TEST-BLINDED.csv")
    cfg["aliases"] = add(RUNS / "20260909T132734119233Z-a7063ccd/outputs/aliases.parquet")
    cfg["foundation"] = [
        add(p)
        for p in sorted(
            (RUNS / "20260910T012630927160Z-e9576045/outputs").glob("chemeleon-*.parquet")
        )
    ]
    cfg["mordred"] = [
        add(p)
        for p in sorted(
            (RUNS / "20260909T132918829974Z-bc81af58/outputs").glob("mordred-*.parquet")
        )
    ]
    cfg["weights"] = {
        k: add(Path.home() / f".cache/tabpfn/tabpfn-v2.6-{k}-v2.6_default.ckpt")
        for k in ["classifier", "regressor"]
    }
    for fold, (run, stack) in enumerate(zip(READOUT, STACK, strict=True)):
        source = RUNS / run / "outputs"
        parent = RUNS / stack / "outputs"
        summary = json.loads((parent / "summary.json").read_text())
        item = dict(
            fold=fold,
            run=run,
            stack_run=stack,
            readout_summary=add(source / "summary.json"),
            stack_summary=add(parent / "summary.json"),
            oof=add(source / "oof.parquet"),
            auxiliary=[],
            endpoints={},
        )
        for aux in summary["config"]["auxiliary"]:
            base = RUNS / aux["run"] / "outputs"
            item["auxiliary"].append(
                dict(
                    run=aux["run"],
                    summary=add(base / "summary.json"),
                    predictions=[add(p) for p in sorted(base.glob("predictions-*.parquet"))],
                )
            )
        for ep in summary["config"]["endpoints"]:
            item["endpoints"][ep] = dict(
                model=add(source / f"{ep}.tabpfn_fit"),
                columns=add(source / f"{ep}-columns.json"),
                imputer=add(source / f"{ep}-imputer.joblib"),
                pca=add(parent / f"{ep}-pca.joblib"),
                validation_features=add(parent / f"{ep}-features.parquet"),
            )
        cfg["folds"].append(item)
    path = ROOT / "configs/submission/baseline-inference.json"
    path.parent.mkdir(exist_ok=True)
    assert not path.exists()
    path.write_text(json.dumps(cfg, indent=2) + "\n")
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        "Baseline submission inference: five fold TabPFN 2.6 ensemble; saved OOF reload audit",
        "--config",
        str(path),
    ]
    for p in files:
        cmd += ["--input", str(p)]
    subprocess.run(cmd + ["--", "python", "scripts/baseline_inference.py"], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
