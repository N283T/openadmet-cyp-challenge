"""Set up and run resumable CYP trunk extraction on a Colab A100.

The script is sent to an already authenticated Colab session with ``colab exec``.
It keeps the system PyTorch/CUDA build, installs Boltz in an isolated venv, runs
small local chunks, and only removes a chunk after its archive is durable on Drive.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time
from collections import deque
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path

DRIVE_ROOT = Path("/content/drive/MyDrive/openadmet-cyp-boltz-colab")
WORK_ROOT = Path("/content/openadmet-cyp-boltz-colab")
TARGET = os.environ.get("COLAB_TARGET", "cyp1a2").lower()
EXPECTED_RECORDS = {
    "cyp1a2": 2163,
    "cyp3a4": 4334,
}[TARGET]
SOURCE_ARCHIVE = DRIVE_ROOT / "boltz-trunk-source.tar.gz"
INPUT_ARCHIVE = DRIVE_ROOT / f"{TARGET}-inputs.tar.gz"
CHECKPOINT_SOURCE = DRIVE_ROOT / "boltz2_conf.ckpt"
SOURCE_ROOT = WORK_ROOT / "boltz-trunk"
INPUT_ROOT = WORK_ROOT / TARGET
CHECKPOINT = WORK_ROOT / "boltz2_conf.ckpt"
VENV_ROOT = WORK_ROOT / "venv"
RESULT_ROOT = DRIVE_ROOT / f"{TARGET}-r1-a100"

EXPECTED_SHA256 = {
    "boltz-trunk-source.tar.gz": "e4dd4fa861a7cb8114d1e4edc8bb0ca25f1531dbf2928172881a487b1b7cb2a7",
    "boltz2_conf.ckpt": "090e82ac8c92f5e943fa1b39e7410a44027bea7243c0bbb3caa67a77fc1428e1",
}

DEPENDENCIES = [
    "numpy>=1.26",
    "pytorch-lightning>=2.5.0",
    "lightning-utilities>=0.10.0",
    "torchmetrics>=0.7.0",
    "rdkit>=2024.3.2",
    "requests>=2.32.3",
    "pandas>=2.2.2",
    "types-requests",
    "einops>=0.8.0",
    "mashumaro>=3.14",
    "click>=8.1.7",
    "pyyaml>=6.0.2",
    "biopython>=1.84",
    "scipy>=1.13.1",
    "numba>=0.61.0",
    "llvmlite>=0.49.0",
    "gemmi>=0.6.5",
    "scikit-learn>=1.6.1",
    "chembl_structure_pipeline>=1.2.2",
    "omegaconf>=2.3.0",
    "cuequivariance>=0.5.0",
    "cuequivariance_ops_cu12>=0.5.0",
    "cuequivariance_ops_torch_cu12>=0.5.0",
    "cuequivariance_torch>=0.5.0",
]


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command: list[str], *, cwd: Path | None = None) -> None:
    log("RUN " + " ".join(command))
    process = subprocess.Popen(
        command,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert process.stdout is not None
    for line in process.stdout:
        print(line, end="", flush=True)
    return_code = process.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, command)


def verify_sources() -> None:
    for path in (SOURCE_ARCHIVE, INPUT_ARCHIVE, CHECKPOINT_SOURCE):
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = sha256(path)
        if path == INPUT_ARCHIVE:
            expected = os.environ.get("COLAB_INPUT_SHA256")
            if not expected:
                raise ValueError("COLAB_INPUT_SHA256 is required")
        else:
            expected = EXPECTED_SHA256[path.name]
        if actual != expected:
            raise ValueError(f"SHA-256 mismatch for {path}: {actual} != {expected}")
        log(f"verified {path.name}: {actual}")


def safe_extract(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as bundle:
        destination_resolved = destination.resolve()
        for member in bundle.getmembers():
            target = (destination / member.name).resolve()
            if not target.is_relative_to(destination_resolved):
                raise ValueError(f"Unsafe archive member: {member.name}")
        bundle.extractall(destination, filter="data")


def setup() -> None:
    verify_sources()
    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    if not (SOURCE_ROOT / "pyproject.toml").is_file():
        SOURCE_ROOT.mkdir(parents=True, exist_ok=True)
        safe_extract(SOURCE_ARCHIVE, SOURCE_ROOT)
    if not (INPUT_ROOT / "target.a3m").is_file():
        safe_extract(INPUT_ARCHIVE, WORK_ROOT)
    if not CHECKPOINT.is_file() or sha256(CHECKPOINT) != EXPECTED_SHA256[CHECKPOINT_SOURCE.name]:
        temp_checkpoint = CHECKPOINT.with_suffix(".ckpt.partial")
        shutil.copyfile(CHECKPOINT_SOURCE, temp_checkpoint)
        os.replace(temp_checkpoint, CHECKPOINT)

    ready_marker = VENV_ROOT / ".boltz-ready"
    if not ready_marker.is_file():
        if VENV_ROOT.exists():
            shutil.rmtree(VENV_ROOT)
        run(
            [
                "uv",
                "venv",
                "--system-site-packages",
                "--python",
                sys.executable,
                str(VENV_ROOT),
            ]
        )
    python = str(VENV_ROOT / "bin/python")
    # uv does not consider system-site packages during dependency resolution and
    # would otherwise replace Colab's CUDA 12 torch with a CUDA 13 build.
    run(["uv", "pip", "install", "--python", python, "--no-deps", *DEPENDENCIES])
    run(
        [
            "uv",
            "pip",
            "install",
            "--python",
            python,
            "--no-deps",
            "-e",
            str(SOURCE_ROOT),
        ]
    )
    run(
        [
            python,
            "-c",
            (
                "import torch, cuequivariance_torch; "
                "assert torch.cuda.is_available(); "
                "assert torch.version.cuda.startswith('12.'); "
                "print({'torch': torch.__version__, 'cuda': torch.version.cuda, "
                "'gpu': torch.cuda.get_device_name(0)})"
            ),
        ]
    )
    ready_marker.touch()
    metadata = {
        "schema_version": 1,
        "target": TARGET,
        "source_sha256": EXPECTED_SHA256[SOURCE_ARCHIVE.name],
        "inputs_sha256": os.environ["COLAB_INPUT_SHA256"],
        "checkpoint_sha256": EXPECTED_SHA256[CHECKPOINT_SOURCE.name],
        "python": sys.version,
    }
    RESULT_ROOT.mkdir(parents=True, exist_ok=True)
    (RESULT_ROOT / "setup.json").write_text(json.dumps(metadata, indent=2) + "\n")
    log("setup complete")


def source_records() -> list[Path]:
    records = list((INPUT_ROOT / "inputs").glob("*.yaml"))
    records.sort(key=lambda path: int(path.name.split("_", 1)[0]))
    if len(records) != EXPECTED_RECORDS:
        raise ValueError(
            f"Expected {EXPECTED_RECORDS} {TARGET.upper()} inputs, found {len(records)}"
        )
    return records


def ensure_structure_free_conformer_fallback() -> None:
    """Hot-patch an older Drive source snapshot with the validated 2D fallback."""
    schema = SOURCE_ROOT / "src/boltz/data/parse/schema.py"
    source = schema.read_text()
    marker = 'conformer.SetProp("coord_generation", "RDKit2DFallback")'
    if marker in source:
        return
    needle = """    if conf_id != -1:\n        conformer = mol.GetConformer(conf_id)\n"""
    replacement = (
        "    if conf_id == -1:\n"
        "        print(\n"
        '            "WARNING: RDKit 3D embedding failed after the random-coordinate "\n'
        '            "retry; using deterministic 2D coordinates as a structure-free "\n'
        '            "fallback."\n'
        "        )\n"
        "        conf_id = AllChem.Compute2DCoords(mol, clearConfs=False)\n"
        "        conformer = mol.GetConformer(conf_id)\n"
        "        conformer.Set3D(False)\n"
        '        conformer.SetProp("name", "Computed2DFallback")\n'
        '        conformer.SetProp("coord_generation", "RDKit2DFallback")\n'
        "        return True\n"
        "\n"
        "    if conf_id != -1:\n"
        "        conformer = mol.GetConformer(conf_id)\n"
    )
    if source.count(needle) != 1:
        raise ValueError("Could not apply the validated conformer fallback to Colab source")
    schema.write_text(source.replace(needle, replacement, 1))
    log("applied validated deterministic RDKit 2D conformer fallback")


def prepare_records(paths: list[Path], destination: Path) -> list[str]:
    import yaml

    destination.mkdir(parents=True, exist_ok=False)
    record_ids = []
    for source in paths:
        record_number = int(source.name.split("_", 1)[0])
        record_id = f"record_{record_number:06d}"
        payload = yaml.safe_load(source.read_text())
        proteins = [item["protein"] for item in payload["sequences"] if "protein" in item]
        if len(proteins) != 1:
            raise ValueError(f"Expected one protein in {source}")
        proteins[0]["msa"] = str((INPUT_ROOT / "target.a3m").resolve())
        (destination / f"{record_id}.yaml").write_text(
            yaml.safe_dump(payload, sort_keys=False, default_flow_style=False)
        )
        record_ids.append(record_id)
    return record_ids


def validate_outputs(output_dir: Path, record_ids: list[str]) -> None:
    import numpy as np

    actual = sorted(
        path.stem.removeprefix("embeddings_") for path in output_dir.glob("embeddings_*.npz")
    )
    if actual != sorted(record_ids):
        raise ValueError(f"Output record mismatch: expected {len(record_ids)}, got {len(actual)}")
    for record_id in record_ids:
        with np.load(output_dir / f"embeddings_{record_id}.npz", allow_pickle=False) as data:
            if data["s"].dtype != np.float16 or data["z"].dtype != np.float16:
                raise ValueError(f"Unexpected tensor dtype for {record_id}")
            if int(data["recycling_steps"]) != 1:
                raise ValueError(f"Unexpected recycling steps for {record_id}")
            if str(data["embedding_stage"]) != "post":
                raise ValueError(f"Unexpected embedding stage for {record_id}")
            if str(data["checkpoint_sha256"]) != EXPECTED_SHA256[CHECKPOINT_SOURCE.name]:
                raise ValueError(f"Unexpected checkpoint for {record_id}")
    config = json.loads((output_dir / "trunk_run_config.json").read_text())
    if config["recycling_steps"] != 1 or not config["use_kernels"]:
        raise ValueError(f"Unexpected runtime config: {config}")


def build_archive(
    output_dir: Path,
    chunk_name: str,
    record_ids: list[str],
) -> tuple[Path, str, Path]:
    archive = WORK_ROOT / f"{chunk_name}.tar"
    with tarfile.open(archive, "w") as bundle:
        for path in sorted(output_dir.glob("embeddings_*.npz")):
            bundle.add(path, arcname=path.name)
        bundle.add(output_dir / "trunk_manifest.json", arcname="trunk_manifest.json")
        bundle.add(output_dir / "trunk_run_config.json", arcname="trunk_run_config.json")
        manifest = WORK_ROOT / f"{chunk_name}.json"
        manifest.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "chunk": chunk_name,
                    "record_ids": record_ids,
                    "record_count": len(record_ids),
                    "archive": archive.name,
                },
                indent=2,
            )
            + "\n"
        )
        bundle.add(manifest, arcname="colab_chunk_manifest.json")
    return archive, sha256(archive), manifest


def upload_archive(
    archive: Path,
    archive_hash: str,
    manifest: Path,
    chunk_root: Path,
    record_count: int,
) -> None:
    partial = RESULT_ROOT / f".{archive.name}.partial"
    shutil.copyfile(archive, partial)
    if sha256(partial) != archive_hash:
        raise ValueError(f"Drive copy verification failed for {archive.name}")
    final = RESULT_ROOT / archive.name
    os.replace(partial, final)
    (RESULT_ROOT / f"{archive.name}.sha256").write_text(f"{archive_hash}  {archive.name}\n")
    log(f"saved {final} ({record_count} records, sha256={archive_hash})")
    archive.unlink()
    manifest.unlink()
    shutil.rmtree(chunk_root)


def extract_chunk(
    paths: list[Path],
    chunk_name: str,
    uploader: ThreadPoolExecutor,
) -> Future[None] | None:
    final = RESULT_ROOT / f"{chunk_name}.tar"
    checksum = RESULT_ROOT / f"{chunk_name}.tar.sha256"
    if final.is_file() and checksum.is_file():
        expected = checksum.read_text().split()[0]
        if sha256(final) == expected:
            log(f"skip completed {chunk_name}")
            return None
        raise ValueError(f"Existing Drive archive failed verification: {final}")

    chunk_root = WORK_ROOT / "chunks" / chunk_name
    if chunk_root.exists():
        shutil.rmtree(chunk_root)
    input_dir = chunk_root / "inputs"
    output_dir = chunk_root / "outputs"
    record_ids = prepare_records(paths, input_dir)
    command = [
        str(VENV_ROOT / "bin/boltz-stfree"),
        "extract-trunk",
        str(input_dir),
        "--output-dir",
        str(output_dir),
        "--checkpoint",
        str(CHECKPOINT),
        "--recycling-steps",
        "1",
        "--device",
        "cuda",
        "--max-msa-seqs",
        "8192",
        "--subsample-msa",
        "--num-subsampled-msa",
        "1024",
        "--seed",
        "20260928",
        "--ligand-conformer-seed",
        "20260928",
        "--embedding-stage",
        "post",
        "--output-dtype",
        "float16",
    ]
    run(command, cwd=SOURCE_ROOT)
    validate_outputs(output_dir, record_ids)
    archive, archive_hash, manifest = build_archive(output_dir, chunk_name, record_ids)
    return uploader.submit(
        upload_archive,
        archive,
        archive_hash,
        manifest,
        chunk_root,
        len(record_ids),
    )


def extract(smoke: bool) -> None:
    verify_sources()
    if not (VENV_ROOT / "bin/boltz-stfree").is_file():
        raise FileNotFoundError("Run COLAB_MODE=setup first")
    ensure_structure_free_conformer_fallback()
    records = source_records()
    start_record = int(os.environ.get("COLAB_START_RECORD", "1"))
    end_record = int(os.environ.get("COLAB_END_RECORD", str(EXPECTED_RECORDS)))
    records = [
        path for path in records if start_record <= int(path.name.split("_", 1)[0]) <= end_record
    ]
    if not records:
        raise ValueError(f"Empty record range: {start_record}..{end_record}")
    log(f"selected records {start_record}..{end_record} ({len(records)} total)")
    pending: deque[Future[None]] = deque()
    with ThreadPoolExecutor(max_workers=1, thread_name_prefix="drive-upload") as uploader:
        if smoke:
            future = extract_chunk(records[:2], "smoke-000001-000002", uploader)
            if future is not None:
                future.result()
            return
        chunk_size = int(os.environ.get("COLAB_CHUNK_SIZE", "250"))
        for start in range(0, len(records), chunk_size):
            while pending and pending[0].done():
                pending.popleft().result()
            # Bound local storage while still allowing the next GPU chunk to run
            # concurrently with the preceding Drive copy and verification.
            if len(pending) >= 2:
                pending.popleft().result()
            selected = records[start : start + chunk_size]
            first = int(selected[0].name.split("_", 1)[0])
            last = int(selected[-1].name.split("_", 1)[0])
            future = extract_chunk(
                selected,
                f"chunk-{first:06d}-{last:06d}",
                uploader,
            )
            if future is not None:
                pending.append(future)
        while pending:
            pending.popleft().result()
    log(f"all {TARGET.upper()} chunks complete")


def main() -> None:
    mode = os.environ.get("COLAB_MODE", "setup")
    if mode == "inspect":
        log(f"python={sys.executable}")
        for executable in ("python", "python3", "pip", "pip3", "uv"):
            log(f"{executable}={shutil.which(executable)}")
        for command in ([sys.executable, "-m", "pip", "--version"], ["pip", "--version"]):
            completed = subprocess.run(command, check=False, capture_output=True, text=True)
            log(
                f"diagnostic rc={completed.returncode}: "
                f"stdout={completed.stdout.strip()!r} stderr={completed.stderr.strip()!r}"
            )
    elif mode == "setup":
        setup()
    elif mode == "smoke":
        extract(smoke=True)
    elif mode == "full":
        extract(smoke=False)
    elif mode == "results":
        for archive in sorted(RESULT_ROOT.glob("chunk-*.tar")):
            checksum = archive.with_suffix(archive.suffix + ".sha256")
            log(
                f"result {archive.name} bytes={archive.stat().st_size} "
                f"checksum={checksum.is_file()}"
            )
    else:
        raise ValueError(f"Unknown COLAB_MODE: {mode}")


if __name__ == "__main__":
    main()
