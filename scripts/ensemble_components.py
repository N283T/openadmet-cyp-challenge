"""Three separate regression components on matched parent-retained folds."""

import json
import os
from pathlib import Path

import joblib
import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.impute import SimpleImputer
from tabpfn import TabPFNRegressor
from tabpfn.constants import ModelVersion

from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
fold = cfg["fold"]
torch.set_num_threads(4)


def read(name):
    return pl.read_parquet(next(inp.glob("*-" + name + ".parquet")))


a, labels, compounds = read("parent-assignments"), read("measurements"), read("compounds")
emb = pl.concat(pl.read_parquet(p) for p in sorted(inp.glob("*-chemeleon-*.parquet")))
desc = pl.concat(pl.read_parquet(p) for p in sorted(inp.glob("*-mordred-*.parquet")))
ecols = sorted(c for c in emb.columns if c != "compound_id")
dcols = [c for c in desc.columns if c != "compound_id"]
auxiliary = []
for item in cfg["auxiliary"]:
    meta = json.loads((inp / item["summary"]).read_text())
    pair = meta["config"]["excluded_folds"]
    assert fold in pair
    members = pl.read_parquet(inp / item["membership"]).join(
        a.select("compound_id", "fold"), on="compound_id", validate="1:1"
    )
    assert (members["fold"].is_in(pair) == (members["role"] == "outer_holdout")).all()
    auxiliary.append(
        (
            next(i for i in pair if i != fold),
            pl.concat(pl.read_parquet(inp / p) for p in item["predictions"]),
            meta["heads"],
        )
    )
assert len(auxiliary) == 4
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
rows, audits = [], []
for ep in sorted(
    labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))["endpoint"].unique()
):
    f = (
        labels.filter(pl.col("endpoint") == ep)
        .join(a, on="compound_id", validate="1:1")
        .join(emb, on="compound_id", validate="1:1")
        .join(desc, on="compound_id", validate="1:1")
        .join(compounds, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    y = f["value"].to_numpy()
    ff = f["fold"].to_numpy()
    va = (ff == fold) & ~f["is_anchor"].to_numpy()
    tr = ~va
    aux = np.zeros((f.height, 4), dtype=np.float32)
    for other, pred, heads in auxiliary:
        v = (
            f.select("compound_id")
            .join(pred, on="compound_id", validate="1:1", maintain_order="left")
            .select(heads)
            .to_numpy()
        )
        aux[ff == other] = v[ff == other]
        aux[ff == fold] += v[ff == fold] / 4
    dx = f.select(dcols).to_numpy().astype(float)
    dx[~np.isfinite(dx)] = np.nan
    keep = np.isfinite(dx[tr]).any(axis=0) & (np.nan_to_num(np.nanstd(dx[tr], axis=0), nan=0) > 0)
    imputer = SimpleImputer(strategy="median").fit(dx[tr][:, keep])
    dd = imputer.transform(dx[:, keep])
    fp = np.stack(
        [fpgen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in f["standardized_smiles"]]
    )
    xs = {
        "fp_mordred_lgbm": np.column_stack([fp, dd]),
        "emb_tabpfn": f.select(ecols).to_numpy().astype(np.float32),
        "pred_log2fc_tabpfn": aux,
    }
    joblib.dump(
        dict(imputer=imputer, columns=[c for c, k in zip(dcols, keep, strict=True) if k]),
        out / f"{ep}-preprocessing.joblib",
    )
    for name, x in xs.items():
        if name == "fp_mordred_lgbm":
            model = lgb.LGBMRegressor(
                objective="regression_l1",
                n_estimators=300,
                num_leaves=15,
                min_child_samples=20,
                learning_rate=0.05,
                reg_lambda=1,
                n_jobs=4,
                random_state=42 + fold,
                verbosity=-1,
                deterministic=True,
                force_col_wise=True,
            )
        else:
            model = TabPFNRegressor.create_default_for_version(
                ModelVersion.V2_6,
                n_estimators=8,
                model_path=str(inp / cfg["weights"]),
                device="cuda",
                random_state=42 + fold,
                auto_scale_n_estimators=False,
                n_preprocessing_jobs=1,
                memory_saving_mode="auto",
                show_progress_bar=False,
                ignore_pretraining_limits=True,
                fit_mode="fit_preprocessors",
            )
        model.fit(x[tr], y[tr])
        pred = model.predict(x[va])
        assert np.isfinite(pred).all()
        if name == "fp_mordred_lgbm":
            model.booster_.save_model(str(out / f"{ep}-{name}.txt"))
        else:
            model.save_fit_state(out / f"{ep}-{name}.tabpfn_fit")
        rows.append(
            f.filter(pl.Series(va))
            .select(
                "compound_id",
                "endpoint",
                "fold",
                "is_neighbor",
                "conf_low",
                "conf_high",
                pl.col("value").alias("y_true"),
            )
            .with_columns(
                pl.lit(name).alias("method"),
                pl.lit(cfg["split_run"]).alias("split_id"),
                    pl.Series("prediction", pred.astype(np.float64)),
            )
        )
        audits.append(
            dict(
                endpoint=ep,
                method=name,
                fold=fold,
                train_n=int(tr.sum()),
                validation_n=int(va.sum()),
                features=x.shape[1],
            )
        )
        print(json.dumps(audits[-1]), flush=True)
pl.concat(rows).write_parquet(out / "component-oof.parquet")
write_json(out / "summary.json", dict(config=cfg, fits=audits))
