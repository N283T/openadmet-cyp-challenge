"""Full-fit selected TabPFN 3.5 atomic components for all regression endpoints."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs"
VIEWS = RUNS / "regression/20260917T111518065286Z-801ab060/outputs"
AUXILIARY = RUNS / "shared/20260912T074339967571Z-ff194fcf/outputs/full-aux.parquet"
FAST_WEIGHTS = Path.home() / ".cache/tabpfn/tabpfn-v3.5-fast-20260909.safetensors"
FULL_WEIGHTS = Path.home() / ".cache/tabpfn/tabpfn-v3.5-20260909.safetensors"

views = [
    ("pub030_raw", VIEWS / "jeremy_pubchem_anchor030-raw.parquet", "raw_", [], 600, None),
    (
        "jeremy_raw",
        RUNS / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
        "emb_",
        [],
        600,
        None,
    ),
    (
        "anchor59_raw_scores_logfc",
        VIEWS / "jeremy_init_task_anchor-raw-scores.parquet",
        "raw_",
        ["score_"],
        659,
        AUXILIARY,
    ),
    (
        "jeremy_scores_logfc",
        RUNS / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
        "score_",
        [],
        47,
        AUXILIARY,
    ),
]

for version_name, version, weights in [
    ("fast", "v3.5-fast", FAST_WEIGHTS),
    ("full", "v3.5", FULL_WEIGHTS),
]:
    for view_name, features, prefix, extra_prefixes, dimensions, auxiliary in views:
        method = f"tabpfn35_{version_name}_{view_name}_fullfit_v1"
        completed = False
        for result_path in (RUNS / "regression").glob("*/result.json"):
            result = __import__("json").loads(result_path.read_text())
            config_path = result_path.parent / "config.json"
            if result.get("status") == "completed" and config_path.exists():
                completed = (
                    __import__("json").loads(config_path.read_text()).get("method") == method
                )
                if completed:
                    break
        if completed:
            continue
        command = [
            "pixi",
            "run",
            "python",
            "scripts/launch_full_submission_readout.py",
            "--method",
            method,
            "--features",
            str(features),
            "--prefix",
            prefix,
            "--dimensions",
            str(dimensions),
            "--backend",
            "tabpfn",
            "--tabpfn-version",
            version,
            "--tabpfn-weights",
            str(weights),
        ]
        for item in extra_prefixes:
            command.extend(["--extra-prefix", item])
        if auxiliary is not None:
            command.extend(["--auxiliary", str(auxiliary)])
        subprocess.run(command, cwd=ROOT, check=True)
