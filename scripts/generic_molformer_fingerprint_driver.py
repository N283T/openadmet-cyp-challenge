"""Extract frozen generic MoLFormer-c3 pooled representations to NPZ."""

import argparse

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset
from transformers import AutoModel, AutoTokenizer


class SmilesDataset(Dataset):
    def __init__(self, smiles: list[str], tokenizer, max_length: int):
        self.smiles = smiles
        self.tokenizer = tokenizer
        self.max_length = max_length

    def __len__(self) -> int:
        return len(self.smiles)

    def __getitem__(self, index: int) -> dict[str, torch.Tensor]:
        encoded = self.tokenizer(
            self.smiles[index],
            padding="max_length",
            truncation=True,
            max_length=self.max_length,
            return_tensors="pt",
        )
        return {key: value.squeeze(0) for key, value in encoded.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-path", required=True)
    parser.add_argument("--model-path", required=True)
    parser.add_argument("--output-path", required=True)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--pooling", choices=["cls", "mean", "cls_mean"], default="cls")
    args = parser.parse_args()
    torch.manual_seed(42)
    frame = pd.read_csv(args.data_path)
    if args.limit:
        frame = frame.head(args.limit)
    tokenizer = AutoTokenizer.from_pretrained(
        args.model_path, trust_remote_code=True, local_files_only=True
    )
    model = AutoModel.from_pretrained(
        args.model_path, trust_remote_code=True, local_files_only=True
    )
    for layer in model.encoder.layer:
        rotary = layer.attention.self.rotary_embeddings
        rotary.inv_freq = 1.0 / (
            rotary.base ** (torch.arange(0, rotary.dim, 2).float() / rotary.dim)
        )
        rotary._set_cos_sin_cache(
            seq_len=rotary.max_position_embeddings,
            device=rotary.inv_freq.device,
            dtype=torch.float32,
        )
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device).eval()
    loader = DataLoader(
        SmilesDataset(frame["smiles"].tolist(), tokenizer, 202),
        batch_size=args.batch_size,
        shuffle=False,
    )
    chunks = []
    with torch.no_grad():
        for batch in loader:
            batch = {key: value.to(device) for key, value in batch.items()}
            hidden = model(**batch).last_hidden_state
            cls = hidden[:, 0, :]
            mask = batch["attention_mask"].unsqueeze(-1).to(hidden.dtype)
            mean = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp_min(1)
            pooled = {"cls": cls, "mean": mean, "cls_mean": torch.cat([cls, mean], dim=1)}
            chunks.append(pooled[args.pooling].cpu().numpy())
    fingerprints = np.concatenate(chunks, axis=0)
    if not np.isfinite(fingerprints).all():
        raise ValueError("Non-finite MoLFormer fingerprints")
    np.savez_compressed(args.output_path, fps=fingerprints)
    print(f"MoLFormer fingerprints: {fingerprints.shape}", flush=True)


if __name__ == "__main__":
    main()
