"""Challenge-disjoint substrate classifiers; frozen six-score auxiliary features."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

from cyp_challenge.preprocessing import standardize

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


RDLogger.DisableLog("rdApp.*")
challenge = pl.read_parquet(path(0)).sort("compound_id")
fp = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def conn(s):
    return Chem.MolToInchiKey(Chem.MolFromSmiles(s)).split("-")[0]


def fingerprints(smiles):
    return np.asarray([fp.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in smiles])


forbidden = {conn(s) for s in challenge["standardized_smiles"]}
xchallenge = fingerprints(challenge["standardized_smiles"])
features = challenge.select("compound_id")
records = []
for i, iso in enumerate(cfg["isoforms"], 1):
    for r in pl.read_csv(path(i), separator="\t").iter_rows(named=True):
        s = standardize(r["Drug"])
        c = conn(s["standardized_smiles"]) if s["compound_id"] else None
        records.append(
            dict(
                isoform=iso,
                source_name=r["Drug_ID"],
                raw_smiles=r["Drug"],
                y=int(r["Y"]),
                compound_id=s["compound_id"],
                standardized_smiles=s["standardized_smiles"],
                connectivity=c,
                status=s["status"],
                overlap=c in forbidden,
            )
        )
raw = pl.DataFrame(records)
raw.write_parquet(out / "curation-records.parquet")
valid = raw.filter(pl.col("connectivity").is_not_null() & ~pl.col("overlap"))
conflicts = (
    valid.group_by("isoform", "connectivity")
    .agg(pl.col("y").n_unique().alias("labels"))
    .filter(pl.col("labels") > 1)
)
conflicts.write_parquet(out / "conflicts.parquet")
clean = (
    valid.join(
        conflicts.select("isoform", "connectivity"), on=["isoform", "connectivity"], how="anti"
    )
    .unique(["isoform", "connectivity"])
    .sort("isoform", "connectivity")
)
assert set(clean["connectivity"]).isdisjoint(forbidden)
clean.write_parquet(out / "external-substrates.parquet")
metrics = []
members = []
counts = []
for iso in cfg["isoforms"]:
    h = clean.filter(pl.col("isoform") == iso)
    x = fingerprints(h["standardized_smiles"])
    y = h["y"].to_numpy()
    groups = [
        MurckoScaffold.MurckoScaffoldSmiles(smiles=s) or c
        for s, c in h.select("standardized_smiles", "connectivity").iter_rows()
    ]
    counts.append(
        dict(
            isoform=iso,
            raw=raw.filter(pl.col("isoform") == iso).height,
            retained=h.height,
            positive=int(y.sum()),
            excluded_overlap=raw.filter((pl.col("isoform") == iso) & pl.col("overlap")).height,
        )
    )
    for name in ["logistic", "forest"]:

        def model(name=name):
            return (
                LogisticRegression(C=1.0, max_iter=2000, random_state=42)
                if name == "logistic"
                else RandomForestClassifier(
                    n_estimators=500,
                    min_samples_leaf=3,
                    max_features="sqrt",
                    n_jobs=3,
                    random_state=42,
                )
            )

        pred = np.zeros(len(y))
        for fold, (tr, va) in enumerate(
            StratifiedGroupKFold(3, shuffle=True, random_state=42).split(x, y, groups)
        ):
            m = model().fit(x[tr], y[tr])
            pred[va] = m.predict_proba(x[va])[:, 1]
            joblib.dump(m, out / f"{iso}-{name}-cv{fold}.joblib")
            members.append(
                h.select("compound_id", "isoform").with_columns(
                    pl.lit(name).alias("method"),
                    pl.lit(fold).alias("fold"),
                    pl.Series("validation", np.isin(np.arange(len(y)), va)),
                )
            )
        metrics.append(
            dict(
                isoform=iso,
                method=name,
                n=len(y),
                ROC_AUC=roc_auc_score(y, pred),
                AUPRC=average_precision_score(y, pred),
                positive_fraction=float(y.mean()),
            )
        )
        m = model().fit(x, y)
        joblib.dump(m, out / f"{iso}-{name}-full.joblib")
        features = features.with_columns(
            pl.Series(f"substrate_{iso}_{name}", m.predict_proba(xchallenge)[:, 1])
        )
        print(metrics[-1], flush=True)
features.write_parquet(out / "substrate-features.parquet")
pl.DataFrame(metrics).write_csv(out / "external-cv.csv")
pl.concat(members).write_parquet(out / "external-memberships.parquet")
metadata = json.loads(path(4).read_text())["datasetVersion"]
assert {4259584, 4259578, 4259581} <= {f["dataFile"]["id"] for f in metadata["files"]}
(out / "summary.json").write_text(
    json.dumps(
        dict(
            counts=counts,
            license=metadata["license"],
            config=cfg,
            features=features.columns,
            challenge_exclusion="all 6896 compound connectivity blocks",
            label_semantics="observed substrate classification; neither pIC50 nor atom-level SoM",
        ),
        indent=2,
    )
)
