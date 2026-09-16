"""Launch generic frozen KERMT or MoLFormer feature extraction."""

import argparse
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("encoder", choices=["kermt_v2", "molformer_c3", "molformer_c3_cls_mean"])
parser.add_argument("mode", choices=["smoke", "full"])
args = parser.parse_args()

root = Path(__file__).resolve().parents[1]
compounds = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet"
if args.encoder == "kermt_v2":
    name = "generic-kermt-v2"
    snapshot = (
        root
        / ".cache/huggingface/hub/models--nvidia--NV-KERMT-70M-v2/snapshots"
        / "7df5eb3179235fdea1e8124db73215da33d77dce"
    )
    inputs = [
        compounds,
        snapshot / "kermt_contrastive_v2.0.pt",
        snapshot / "pretrain_atom_vocab.json",
        snapshot / "pretrain_bond_vocab.json",
        snapshot / "pretrain_smiles_vocab.pkl",
    ]
else:
    name = (
        "generic-molformer-c3-cls-mean"
        if args.encoder == "molformer_c3_cls_mean"
        else "generic-molformer-c3"
    )
    snapshot = (
        root
        / ".cache/huggingface/hub/models--DeepChem--MoLFormer-c3-1.1B/snapshots"
        / "9f1b9ea3590833bd0ea1a70e789c5d3da11ba7ed"
    )
    code_snapshot = Path(
        "/home/nagaet/.cache/huggingface/modules/transformers_modules/ibm/"
        "MoLFormer_hyphen_XL_hyphen_both_hyphen_10pct/"
        "7b12d946c181a37f6012b9dc3b002275de070314"
    )
    model_files = [
        "config.json",
        "model.safetensors",
        "special_tokens_map.json",
        "tokenizer.json",
        "tokenizer_config.json",
        "vocab.json",
    ]
    code_files = [
        "__init__.py",
        "configuration_molformer.py",
        "modeling_molformer.py",
        "tokenization_molformer.py",
        "tokenization_molformer_fast.py",
    ]
    inputs = [
        compounds,
        *(snapshot / name for name in model_files),
        *(code_snapshot / name for name in code_files),
    ]

config = root / f"configs/regression/{name}-features-{args.mode}-v1.json"
for path in [config, *inputs]:
    if not path.is_file():
        raise FileNotFoundError(path)

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    f"Extract untouched frozen {args.encoder} molecular features ({args.mode})",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/generic_encoder_features.py"], cwd=root, check=True
)
