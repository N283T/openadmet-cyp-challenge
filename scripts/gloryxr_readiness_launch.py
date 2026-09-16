"""Launch a recorded GLORYxR model/rule readiness audit."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source = Path.home() / "ghq/github.com/molinfo-vienna/GLORYxR"
inputs = [
    source / "README.rst",
    source / "pyproject.toml",
    source / "docs/source/tutorials/python_api.ipynb",
    source / "docs/source/tutorials/example_model.joblib",
    source / "src/gloryxr/rules_data/adams_rules.csv",
    source / "src/gloryxr/rules_data/gloryx_reactionrules_connect.csv",
    source / "src/gloryxr/rules_data/gloryx_reactionrules_fixed.csv",
    source / "src/gloryxr/rules_data/gloryx_reactionrules_orig.csv",
]
config = root / "configs/classification-round/gloryxr-readiness-v1.json"
for path in [config, *inputs]:
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Audit pinned GLORYxR reaction rules and availability of scientifically valid scoring models",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(command + ["--", "python", "scripts/gloryxr_readiness.py"], cwd=root, check=True)
