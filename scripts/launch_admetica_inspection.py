"""Record external checkpoint and source inputs for metadata inspection."""

import subprocess
from pathlib import Path

root = Path("/home/nagaet/ghq/github.com/datagrok-ai/admetica")
files = [
    root / f"ADMET/metabolism/{ep}-inhibitor/{ep}-inhibitor.ckpt"
    for ep in ["cyp1a2", "cyp2c9", "cyp2d6", "cyp3a4"]
]
files += [
    root / "examples/cyp_challenge/models/multitask_pac50.ckpt",
    root / "README.md",
    root / "ADMET/metabolism/metabolism.md",
    root / "examples/cyp_challenge/scripts/train_cyp_multitask.py",
]
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Inspect Admetica CYP models and challenge weights at commit 64eacdb7",
]
for p in files:
    cmd += ["--input", str(p)]
subprocess.run(cmd + ["--", "python", "scripts/inspect_admetica_weights.py"], check=True)
