"""Launch cached-graph MONROE embedding with a given checkpoint."""

import argparse
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("graph_run")
    parser.add_argument("--config", required=True)
    parser.add_argument("--purpose", required=True)
    parser.add_argument("--init-config", required=True)
    parser.add_argument("--init-weights", required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    shared = root / "artifacts/runs/shared"
    inputs = {
        f"graphs-seed{s}.pkl": shared / args.graph_run / f"outputs/graphs-seed{s}.pkl"
        for s in range(1, 6)
    }
    inputs["monroe-config.json"] = Path(args.init_config)
    inputs["monroe-weights.pt"] = Path(args.init_weights)
    command = ["pixi", "run", "cyp", "run", "--track", "shared", "--purpose", args.purpose]
    command += ["--config", str(root / args.config)]
    for alias, path in inputs.items():
        command += ["--input", f"{alias}={path}"]
    python = str(root / ".pixi/envs/monroe/bin/python")
    subprocess.run(command + ["--", python, "scripts/monroe_embed_cached.py"], cwd=root, check=True)


if __name__ == "__main__":
    main()
