"""Faithful Polars reconstruction of public author's cohort and connectivity joins."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, RDLogger

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


def csv(i):
    return pl.read_csv(path(i), infer_schema_length=None)


inc, screen, tdi, emax, chembl, veith = [csv(i) for i in range(6)]
aliases, compounds, labels, testmap, assignments = [pl.read_parquet(path(i)) for i in range(6, 11)]
RDLogger.DisableLog("rdApp.*")


def key(s):
    mol = Chem.MolFromSmiles(s)
    return Chem.MolToInchiKey(mol).split("-")[0] if mol else None


iso = ["cyp3a4", "cyp2c9", "cyp2d6", "cyp1a2"]
pubiso = iso + ["cyp2c19"]
direct = [i + "_pic50_direct_inhibition" for i in iso]
logfc = [i + "_log2fc" for i in iso]
chemcols = [i + "_pic50_chembl" for i in pubiso]
maxcols = [i + "_max_response" for i in pubiso]
name_map = aliases.select("molecule_name", "compound_id").unique()
assert name_map["molecule_name"].n_unique() == name_map.height
allnames = dict(name_map.iter_rows())
assert all(n in allnames for n in inc["molecule_name"])
wide = (
    screen.group_by("molecule_name", "enzyme")
    .agg(pl.col("log2fc_estimate").mean())
    .pivot(on="enzyme", index="molecule_name", values="log2fc_estimate")
    .rename({i.upper(): i + "_log2fc" for i in iso})
)
base = inc.join(
    wide.select("molecule_name", *logfc),
    on="molecule_name",
    how="left",
    validate="1:1",
    maintain_order="left",
).with_columns(pl.col("smiles").map_elements(key, return_dtype=pl.String).alias("connectivity"))
base = base.with_columns(pl.lit("challenge_inhibition").alias("source"))
assert base.height == 4905
c = (
    chembl.select("chembl_id", "smiles", "inchi_key", *[i + "_pic50" for i in pubiso])
    .rename({i + "_pic50": i + "_pic50_chembl" for i in pubiso})
    .with_columns(pl.col("inchi_key").str.split("-").list.first().alias("connectivity"))
    .drop("inchi_key")
    .drop_nulls("connectivity")
    .unique("connectivity", keep="first", maintain_order=True)
)
v = (
    veith.group_by("smiles", "isoform")
    .agg(pl.col("max_response").mean())
    .pivot(on="isoform", index="smiles", values="max_response")
    .rename({i: i + "_max_response" for i in pubiso})
    .sort("smiles")
    .with_columns([pl.col(col).clip(-150, 50) for col in maxcols])
    .with_columns(pl.col("smiles").map_elements(key, return_dtype=pl.String).alias("connectivity"))
    .drop_nulls("connectivity")
    .unique("connectivity", keep="first", maintain_order=True)
)


def union(challenge):
    ck = set(challenge["connectivity"])
    new = (
        c.filter(~pl.col("connectivity").is_in(list(ck)))
        .with_columns(
            pl.col("chembl_id").alias("molecule_name"), pl.lit("chembl37").alias("source")
        )
        .drop("chembl_id")
    )
    merged = challenge.join(
        c.select("connectivity", *chemcols),
        on="connectivity",
        how="left",
        validate="m:1",
        maintain_order="left",
    )
    u = pl.concat([merged, new], how="diagonal_relaxed")
    vv = v.filter(~pl.col("connectivity").is_in(u["connectivity"].unique().implode())).with_columns(
        (pl.lit("VEITH-") + pl.col("connectivity")).alias("molecule_name"),
        pl.lit("veith").alias("source"),
    )
    u = u.join(
        v.select("connectivity", *maxcols),
        on="connectivity",
        how="left",
        validate="m:1",
        maintain_order="left",
    )
    return pl.concat([u, vv], how="diagonal_relaxed")


def identify(table):
    # Public rows retain author names and distinct stable IDs. All CV masking uses connectivity.
    return table.with_columns(
        pl.col("molecule_name")
        .map_elements(
            lambda n: allnames.get(n, "scp_" + hashlib.sha256(n.encode()).hexdigest()),
            return_dtype=pl.String,
        )
        .alias("compound_id")
    ).with_columns(
        pl.col("compound_id").is_in(compounds["compound_id"].implode()).alias("is_challenge")
    )


u18 = union(base)
# Source joins Emax onto TDI by molecule_name; the Emax-only row is excluded.
tcols = [i + "_pic50_tdi_condition" for i in iso]
ecols = [i + "_emax_vs_pos_ctrl_direct_inhibition" for i in iso]
side = tdi.select("molecule_name", "smiles", *tcols).join(
    emax.select("molecule_name", *ecols),
    on="molecule_name",
    how="left",
    validate="1:1",
    maintain_order="left",
)
expanded = base.join(
    side.drop("smiles"), on="molecule_name", how="left", validate="1:1", maintain_order="left"
)
sideonly = side.filter(
    ~pl.col("molecule_name").is_in(base["molecule_name"].implode())
).with_columns(
    pl.col("smiles").map_elements(key, return_dtype=pl.String).alias("connectivity"),
    pl.lit("challenge_tdi_only").alias("source"),
)
expanded = pl.concat([expanded, sideonly], how="diagonal_relaxed")
assert expanded.height == 6145
uexpanded = union(expanded)
contracts = {
    "challenge8": (base, direct + logfc),
    "union18": (u18, direct + logfc + chemcols + maxcols),
    "d6_isoform6": (
        uexpanded,
        [
            "cyp2d6_pic50_direct_inhibition",
            "cyp2d6_log2fc",
            "cyp2d6_pic50_tdi_condition",
            "cyp2d6_emax_vs_pos_ctrl_direct_inhibition",
            "cyp2d6_pic50_chembl",
            "cyp2d6_max_response",
        ],
    ),
    "d6_single1": (uexpanded, ["cyp2d6_pic50_direct_inhibition"]),
}
counts = []
contract = {}
idrows = []
for name, (table, heads) in contracts.items():
    table = identify(table)
    assert table["molecule_name"].n_unique() == table.height
    assert table["compound_id"].n_unique() == table.height
    assert table["connectivity"].null_count() == 0
    training = table.filter(pl.any_horizontal([pl.col(h).is_finite() for h in heads]))
    # Keep native wide head names and authors' supplied SMILES, not project restandardization.
    table.write_parquet(out / f"{name}-cohort.parquet")
    training.write_parquet(out / f"{name}.parquet")
    stats = {h: training[h].is_finite().sum() for h in heads}
    contract[name] = dict(
        file=f"{name}.parquet",
        heads=heads,
        source_rows=table.height,
        training_rows=training.height,
        label_counts=stats,
    )
    counts.extend(dict(model=name, head=h, n=n) for h, n in stats.items())
    idrows.append(
        table.select(
            "compound_id", "molecule_name", "smiles", "connectivity", "source", "is_challenge"
        )
    )
    print(name, contract[name], flush=True)
# Verify original values and intervals against current project labels.
checks = []
for i in iso:
    original = identify(base).filter(pl.col(i + "_pic50_direct_inhibition").is_finite())
    expected = labels.filter(pl.col("endpoint") == i.upper() + "_pIC50_direct_inhibition")
    joined = original.join(expected, on="compound_id", validate="1:1")
    assert joined.height == expected.height
    for src, dst in [
        (i + "_pic50_direct_inhibition", "value"),
        (i + "_pic50_direct_inhibition_ci_lower", "conf_low"),
        (i + "_pic50_direct_inhibition_ci_upper", "conf_high"),
    ]:
        assert np.allclose(joined[src].to_numpy(), joined[dst].to_numpy(), rtol=0, atol=1e-9), (
            src,
            np.max(np.abs(joined[src].to_numpy() - joined[dst].to_numpy())),
        )
    checks.append(dict(endpoint=i, labels=joined.height, project_values_bounds_match=True))
allchallenge = compounds.join(
    aliases.select("compound_id", "molecule_name", "raw_smiles").unique(
        "compound_id", keep="first"
    ),
    on="compound_id",
    validate="1:1",
)
source_smiles = identify(pl.concat([base, sideonly], how="diagonal_relaxed")).select(
    "compound_id", pl.col("smiles").alias("author_smiles")
)
allchallenge = (
    allchallenge.join(source_smiles, on="compound_id", how="left", validate="1:1")
    .with_columns(pl.coalesce("author_smiles", "raw_smiles").alias("smiles"))
    .with_columns(pl.col("smiles").map_elements(key, return_dtype=pl.String).alias("connectivity"))
)
allchallenge.select("compound_id", "molecule_name", "smiles", "connectivity").write_parquet(
    out / "prediction-compounds.parquet"
)
allchallenge.select("compound_id", "connectivity").join(
    assignments, on="compound_id", how="left", validate="1:1"
).write_parquet(out / "challenge-connectivity-folds.parquet")
testids = set(testmap["compound_id"])
testkeys = set(allchallenge.filter(pl.col("compound_id").is_in(list(testids)))["connectivity"])
assert len(testids) == 750
pubover = []
for name, t in [("chembl37", c), ("veith", v)]:
    pubover.append(
        dict(
            source=name,
            skeletons=t.height,
            overlap_inhibition=t.filter(
                pl.col("connectivity").is_in(base["connectivity"].unique().implode())
            ).height,
            overlap_all_train=t.filter(
                pl.col("connectivity").is_in(
                    allchallenge.filter(~pl.col("compound_id").is_in(list(testids)))["connectivity"]
                    .unique()
                    .implode()
                )
            ).height,
            overlap_test=t.filter(pl.col("connectivity").is_in(list(testkeys))).height,
        )
    )
    t.filter(
        pl.col("connectivity").is_in(allchallenge["connectivity"].unique().implode())
    ).write_parquet(out / f"{name}-challenge-overlaps.parquet")
pl.DataFrame(counts).write_csv(out / "head-counts.csv")
pl.concat(idrows).unique(["compound_id", "smiles"], maintain_order=True).write_parquet(
    out / "identity-provenance.parquet"
)
(out / "contract.json").write_text(
    json.dumps(
        dict(
            models=contract,
            validation_checks=checks,
            public_overlap=pubover,
            source_config=cfg,
            source_skeleton_policy="Author InChiKey block; first public row; shared heads retained",
        ),
        indent=2,
    )
)
