"""Generate KERMT fingerprints while handling the upstream tuple-return loader."""

import random

import numpy as np
import torch
from rdkit import RDLogger


def main() -> None:
    RDLogger.logger().setLevel(RDLogger.CRITICAL)
    from kermt.data import MolCollator, MoleculeDataset
    from kermt.data.torchvocab import MolVocab  # noqa: F401
    from kermt.util.parsing import parse_args
    from kermt.util.utils import create_logger, get_data, load_checkpoint
    from torch.utils.data import DataLoader

    args = parse_args()
    torch.manual_seed(args.seed)
    torch.cuda.manual_seed_all(args.seed)
    np.random.seed(args.seed)
    random.seed(args.seed)
    torch.backends.cudnn.deterministic = True
    torch.use_deterministic_algorithms(mode=True)
    logger = create_logger(name="fingerprint", save_dir=None, quiet=False)
    data = MoleculeDataset(
        get_data(
            path=args.data_path,
            args=args,
            use_compound_names=False,
            max_data_size=float("inf"),
            skip_invalid_smiles=False,
        )
    )
    model, _state = load_checkpoint(
        args.checkpoint_paths[0], cuda=args.cuda, current_args=args, logger=logger
    )
    model.eval()
    args.bond_drop_rate = 0
    loader = DataLoader(
        data,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=0,
        collate_fn=MolCollator(args=args, shared_dict={}),
    )
    chunks = []
    for _, batch, features_batch, _, _ in loader:
        with torch.no_grad():
            chunks.append(model(batch, features_batch).detach().cpu().numpy())
    fingerprints = np.concatenate(chunks, axis=0)
    if not np.isfinite(fingerprints).all():
        raise ValueError("Non-finite KERMT fingerprints")
    np.savez_compressed(args.output_path, fps=fingerprints)
    print(f"KERMT fingerprints: {fingerprints.shape}", flush=True)


if __name__ == "__main__":
    main()
