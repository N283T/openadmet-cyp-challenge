"""Launch readout screens of OpenADMET CheMeleon checkpoints at chosen message-passing depths.

Usage: launch_openadmet_depth_screen.py name:depth [name:depth ...]. Each argument becomes one
recorded run of `regression_embedding_member.py` on that checkpoint's pooled depth-state vectors,
scored against the incumbent parent-cohort OOF with the same folds.
"""

import argparse
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared, regression = root / "artifacts/runs/shared", root / "artifacts/runs/regression"


CONTROLS = {
    "1a2": regression
    / "20260918T034159381060Z-46314b50/outputs/1A2/add_dense_and_fast_pool/oof.parquet",
    "2c9": regression / "20260917T123601175420Z-e1560701/outputs/oof.parquet",
    "2d6": regression / "20260926T085959697226Z-0dbf57d4/outputs/oof.parquet",
    "3a4": regression / "20260917T130718808881Z-dfa02730/outputs/oof.parquet",
}


def depth_runs():
    """Map checkpoint name to the shared run that exported its depths."""
    found = {}
    for run in sorted(shared.iterdir()):
        intent = run / "intent.json"
        if not intent.exists() or not (run / "outputs/depth-6.parquet").exists():
            continue
        purpose = json.loads(intent.read_text()).get("purpose", "")
        if purpose.startswith("Export all six frozen message-passing depths of the "):
            found[purpose.split(" of the ")[1].split(" CheMeleon")[0]] = run
    return found


TABPFN_WEIGHTS = {
    "tabpfn.ckpt": shared
    / "20260910T122343854890Z-3a03f908/inputs/056-tabpfn-v2.6-regressor-v2.6_default.ckpt",
    "tabpfn-v3.5-20260909.safetensors": regression
    / "20260917T122531642138Z-54f25a12/inputs/tabpfn-v3.5-20260909.safetensors",
    "tabpfn-v3.5-fast-20260909.safetensors": regression
    / "20260917T121017505228Z-433d9acb/inputs/tabpfn-v3.5-fast-20260909.safetensors",
}


def launch(name, depth, runs, config, controls):
    settings = json.loads((root / config).read_text())
    weights = settings.get("tabpfn_weights", "tabpfn.ckpt")
    inputs = {
        "embeddings.parquet": runs[name] / f"outputs/depth-{depth}.parquet",
        "test-mapping.parquet": shared
        / "20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
        weights: TABPFN_WEIGHTS[weights],
        "tabicl.ckpt": shared
        / "20260911T093102177307Z-beddd36f/outputs/tabicl-regressor-v2-20260212.ckpt",
    }
    inputs.update({f"control-{key}.parquet": CONTROLS[key] for key in controls})
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Screen {name} CheMeleon depth {depth} {settings.get('purpose_tag', '')}readouts "
        "as regression members",
        "--config",
        str(root / config),
    ]
    for alias, path in inputs.items():
        command += ["--input", f"{alias}={path}"]
    subprocess.run(
        command + ["--", "python", "scripts/regression_embedding_member.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/regression/openadmet-depth-screen-v1.json")
    parser.add_argument("--controls", default="1a2,2c9,2d6,3a4")
    parser.add_argument("tokens", nargs="+", help="name:depth")
    args = parser.parse_args()
    runs = depth_runs()
    for token in args.tokens:
        name, depth = token.split(":")
        launch(name, int(depth), runs, args.config, args.controls.split(","))
