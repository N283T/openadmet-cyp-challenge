"""CPU molecular features for the external broad-panel tree control."""

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator


def featurize(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = generator.GetFingerprintAsNumPy(molecule).astype(np.float32)
    descriptor_values = Descriptors.CalcMolDescriptors(
        molecule, missingVal=float("nan"), silent=True
    )
    descriptor = np.array(
        [descriptor_values[name] for name, _ in Descriptors.descList], dtype=np.float32
    )
    descriptor[~np.isfinite(descriptor)] = np.nan
    return np.concatenate([fp, descriptor])


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def read(name):
        paths = [p for p in inp.iterdir() if p.name.partition("-")[2] == name]
        assert len(paths) == 1
        return pl.read_parquet(paths[0])

    external, challenge = read("external-compounds.parquet"), read("compounds.parquet")
    assert not set(external["compound_id"]) & set(challenge["compound_id"])
    columns = ["compound_id", "standardized_smiles"]
    all_compounds = pl.concat(
        [
            external.select(columns).with_columns(pl.lit("external").alias("role")),
            challenge.select(columns).with_columns(pl.lit("challenge").alias("role")),
        ]
    ).sort("compound_id")
    assert all_compounds["compound_id"].n_unique() == all_compounds.height
    all_compounds.write_parquet(out / "feature-ids.parquet")
    names = [f"morgan_{i:04d}" for i in range(2048)] + [name for name, _ in Descriptors.descList]
    (out / "feature-columns.json").write_text(json.dumps(names))
    with ProcessPoolExecutor(max_workers=cfg.get("descriptor_workers", 2)) as pool:
        features = []
        for index, value in enumerate(
            pool.map(featurize, all_compounds["standardized_smiles"], chunksize=32)
        ):
            features.append(value)
            if (index + 1) % 2000 == 0:
                print(f"Features {index + 1}/{all_compounds.height}", flush=True)
    matrix = np.stack(features)
    assert matrix.shape == (all_compounds.height, len(names))
    assert np.isfinite(matrix[:, :2048]).all()
    np.save(out / "molecular-features.npy", matrix)
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                rows=all_compounds.height,
                external=external.height,
                challenge=challenge.height,
                dimensions=len(names),
                descriptor_nonfinite=int(np.isnan(matrix[:, 2048:]).sum()),
                descriptor_missing_policy="NaN supplied directly to LightGBM",
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
