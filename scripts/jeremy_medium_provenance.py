"""Record public checkpoint history and distinguish upstream versus challenge FT outputs."""

import hashlib
import json
import os
import subprocess
from pathlib import Path

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
repo = Path(cfg["repository_path"])


def git(*args):
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, check=True).stdout


history = git(
    "log", cfg["revision"], "--format=%H %aI %cI %s", "--follow", "--", cfg["checkpoint_path"]
).decode()
(out / "checkpoint-history.txt").write_text(history)
first = history.strip().splitlines()[-1].split()[0]
original_blob = git("show", first + ":" + cfg["checkpoint_path"])
checkpoint = next(p for p in inp.iterdir() if p.name.partition("-")[2] == "chemprop_medium.pt")
assert (
    hashlib.sha256(original_blob).hexdigest() == hashlib.sha256(checkpoint.read_bytes()).hexdigest()
)
original_readme = git("show", first + ":CYP_Challenge/README.md")
(out / "first-publication-README.md").write_bytes(original_readme)
(out / "readme-history.txt").write_bytes(
    git("log", cfg["revision"], "--format=%H %aI %s", "--", "CYP_Challenge/README.md")
)
summary = dict(
    config=cfg,
    first_checkpoint_commit=first,
    history=history.splitlines(),
    initial_and_current_checkpoint_sha256=hashlib.sha256(original_blob).hexdigest(),
    tracked_checkpoint_versions=len(history.strip().splitlines()),
    upstream_input="CHECKPOINTS_DIR / chemprop_medium.pt",
    challenge_finetuned_output="RESULTS_DIR / chemprop_multitask_final_model.pt",
    challenge_finetune_targets=17,
    source_statement="Sibling pretraining: public ChEMBL/TDC/Polaris plus Novartis-surrogate CYP",
    explicit_challenge_pIC50_exclusion_statement_found=False,
    upstream_training_date_available=False,
    conclusion=(
        "Author labels it upstream public-data pretraining; "
        "code keeps separate challenge FT outputs. "
        "No evidence of challenge-label use in this upstream checkpoint was found. "
        "Exact training corpus is unavailable."
    ),
)
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary), flush=True)
