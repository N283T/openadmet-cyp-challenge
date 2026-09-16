"""Run one pinned BioTransformer CYP prediction and inspect its products."""

import hashlib
import json
import os
import shutil
import subprocess
import tarfile
from pathlib import Path

from rdkit import Chem

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
scratch = Path(os.environ["TMPDIR"]) / "biotransformer-probe"
scratch.mkdir(parents=True, exist_ok=False)
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
archive = inp / cfg["archive"]
if hashlib.sha256(archive.read_bytes()).hexdigest() != cfg["archive_sha256"]:
    raise ValueError("BioTransformer source archive hash mismatch")
with tarfile.open(archive) as handle:
    handle.extractall(scratch, filter="data")
jar = scratch / cfg["jar"]
if hashlib.sha256(jar.read_bytes()).hexdigest() != cfg["jar_sha256"]:
    raise ValueError("BioTransformer JAR hash mismatch")
raw_output = scratch / "products.sdf"
command = [
    "java",
    "-Xmx4g",
    "-jar",
    str(jar),
    "-k",
    "pred",
    "-b",
    cfg["mode"],
    "-ismi",
    cfg["smiles"],
    "-osdf",
    str(raw_output),
    "-s",
    str(cfg["steps"]),
    "-cm",
    str(cfg["cyp450_mode"]),
]
completed = subprocess.run(command, cwd=scratch, text=True, capture_output=True, check=True)
shutil.copy2(raw_output, out / "products.sdf")
(out / "stdout.txt").write_text(completed.stdout)
(out / "stderr.txt").write_text(completed.stderr)
supplier = Chem.SDMolSupplier(str(raw_output), removeHs=False)
products = [molecule for molecule in supplier if molecule is not None]
properties = sorted({name for molecule in products for name in molecule.GetPropNames()})
summary = {
    "source_commit": cfg["source_commit"],
    "command": command,
    "product_count": len(products),
    "properties": properties,
    "product_smiles": sorted({Chem.MolToSmiles(molecule) for molecule in products}),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
