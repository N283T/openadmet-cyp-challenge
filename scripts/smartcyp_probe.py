"""Recorded compatibility probe for the vendored SMARTCyp 2.4.2 JAR."""

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
scratch = Path(os.environ["TMPDIR"]) / "smartcyp-probe"
scratch.mkdir(parents=True, exist_ok=False)
config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
jar = inp / config["jar"]

digest = hashlib.sha256(jar.read_bytes()).hexdigest()
if digest != config["jar_sha256"]:
    raise ValueError(f"SMARTCyp JAR hash mismatch: {digest}")

command = [
    "java",
    "-jar",
    str(jar),
    "-printall",
    "-smiles",
    config["probe_smiles"],
]
result = subprocess.run(
    command,
    cwd=scratch,
    text=True,
    capture_output=True,
    timeout=config["timeout_seconds"],
    check=False,
)
(out / "stdout.txt").write_text(result.stdout)
(out / "stderr.txt").write_text(result.stderr)
files = []
for source in sorted(scratch.iterdir()):
    if source.is_file():
        destination = out / source.name
        shutil.copy2(source, destination)
        files.append(
            {
                "name": source.name,
                "bytes": source.stat().st_size,
                "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            }
        )

summary = {
    "command": command,
    "returncode": result.returncode,
    "jar_sha256": digest,
    "source_commit": config["source_commit"],
    "files": files,
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
if result.returncode != 0:
    raise SystemExit(result.returncode)
if not any(item["name"].endswith(".csv") for item in files):
    raise RuntimeError("SMARTCyp completed without a CSV result")
