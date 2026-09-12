"""Ligand-only ETKDGv3 geometry and explicitly non-pKa charge proxies."""

import concurrent.futures
import json
import os
import signal
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import polars as pl
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors

TOPO = [
    "formal_charge",
    "positive_n",
    "basic_n_proxy",
    "aromatic_n",
    "aromatic_atoms",
    "aromatic_rings",
    "rotatable_bonds",
    "n_charge_min",
    "n_charge_max",
    "charge_abs_sum",
]
GEOM = [
    "plane_rms",
    "plane_max",
    "npr1",
    "npr2",
    "radius_gyration",
    "asphericity",
    "eccentricity",
    "aromatic_plane_rms",
    "aromatic_center_max_distance",
    "aromatic_normal_abs_dot_mean",
    "basic_n_aromatic_min_distance",
    "basic_n_aromatic_mean_distance",
    "basic_n_exposure_proxy",
]


def alarm(*_):
    raise TimeoutError("geometry timeout")


def calc(task):
    cid, smi, seeds = task
    base = dict(compound_id=cid)
    raw = []
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise ValueError("invalid SMILES")
        # Non-amide aliphatic N and pyridine-like aromatic N: structural proxy only.
        query = Chem.MolFromSmarts("[$([N;X3,X4;!$(N-[C,S,P]=[O,S,N])]),$([nH0;+0])]")
        basic = sorted({x[0] for x in m.GetSubstructMatches(query)})
        aromatic = [a.GetIdx() for a in m.GetAtoms() if a.GetIsAromatic()]
        rings = [
            list(r)
            for r in m.GetRingInfo().AtomRings()
            if all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r)
        ]
        AllChem.ComputeGasteigerCharges(m)
        qs = np.array([float(a.GetProp("_GasteigerCharge")) for a in m.GetAtoms()])
        nq = qs[[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() == 7]]
        base.update(
            zip(
                TOPO,
                [
                    Chem.GetFormalCharge(m),
                    sum(a.GetAtomicNum() == 7 and a.GetFormalCharge() > 0 for a in m.GetAtoms()),
                    len(basic),
                    sum(a.GetAtomicNum() == 7 and a.GetIsAromatic() for a in m.GetAtoms()),
                    len(aromatic),
                    len(rings),
                    rdMolDescriptors.CalcNumRotatableBonds(m),
                    float(nq.min()) if len(nq) else None,
                    float(nq.max()) if len(nq) else None,
                    float(abs(qs).sum()),
                ],
                strict=True,
            )
        )
        for seed in seeds:
            r = dict(compound_id=cid, seed=seed)
            try:
                signal.signal(signal.SIGALRM, alarm)
                signal.alarm(30)
                mh = Chem.AddHs(m)
                p = AllChem.ETKDGv3()
                p.randomSeed = seed
                p.numThreads = 1
                p.maxIterations = 100
                if AllChem.EmbedMolecule(mh, p) != 0:
                    raise ValueError("ETKDG failed")
                xyz = mh.GetConformer().GetPositions()[: m.GetNumAtoms()]
                centered = xyz - xyz.mean(0)
                normal = np.linalg.svd(centered, full_matrices=False)[2][-1]
                distances = abs(centered @ normal)
                r.update(
                    plane_rms=float(np.sqrt(np.mean(distances**2))),
                    plane_max=float(distances.max()),
                    npr1=rdMolDescriptors.CalcNPR1(mh),
                    npr2=rdMolDescriptors.CalcNPR2(mh),
                    radius_gyration=rdMolDescriptors.CalcRadiusOfGyration(mh),
                    asphericity=rdMolDescriptors.CalcAsphericity(mh),
                    eccentricity=rdMolDescriptors.CalcEccentricity(mh),
                )
                if len(aromatic) >= 3:
                    axyz = xyz[aromatic]
                    ac = axyz - axyz.mean(0)
                    r["aromatic_plane_rms"] = float(
                        np.sqrt(np.mean((ac @ np.linalg.svd(ac, full_matrices=False)[2][-1]) ** 2))
                    )
                centers = [xyz[ring].mean(0) for ring in rings]
                normals = [
                    np.linalg.svd(xyz[ring] - center, full_matrices=False)[2][-1]
                    for ring, center in zip(rings, centers, strict=True)
                ]
                pairs = [(i, j) for i in range(len(rings)) for j in range(i)]
                if pairs:
                    r["aromatic_center_max_distance"] = float(
                        max(np.linalg.norm(centers[i] - centers[j]) for i, j in pairs)
                    )
                    r["aromatic_normal_abs_dot_mean"] = float(
                        np.mean([abs(normals[i] @ normals[j]) for i, j in pairs])
                    )
                ds = [
                    np.linalg.norm(xyz[i] - center)
                    for i in basic
                    for ring, center in zip(rings, centers, strict=True)
                    if i not in ring
                ]
                if ds:
                    r.update(
                        basic_n_aromatic_min_distance=float(min(ds)),
                        basic_n_aromatic_mean_distance=float(np.mean(ds)),
                    )
                # Inverse nonbonded heavy-atom crowding within 4 A; not solvent SASA.
                if basic:
                    dmat = np.linalg.norm(xyz[:, None] - xyz[None, :], axis=2)
                    graph = Chem.GetDistanceMatrix(m)
                    r["basic_n_exposure_proxy"] = float(
                        np.mean([1 / (1 + np.sum((dmat[i] < 4) & (graph[i] > 1))) for i in basic])
                    )
                r["status"] = "ok"
                r["error"] = ""
            except Exception as exc:
                r.update(status="failed", error=str(exc))
            finally:
                signal.alarm(0)
            raw.append(r)
    except Exception as exc:
        raw = [dict(compound_id=cid, seed=s, status="failed", error=str(exc)) for s in seeds]
    for r in raw:
        for k in GEOM:
            r.setdefault(k, None)
    for k in TOPO:
        base.setdefault(k, None)
    return base, raw


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    c = pl.read_parquet(next(inp.glob("000-*")))
    tasks = [
        (r["compound_id"], r["standardized_smiles"], cfg["seeds"]) for r in c.iter_rows(named=True)
    ]
    bases, raw = [], []
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg["workers"]) as pool:
        for i, (b, rs) in enumerate(pool.map(calc, tasks), 1):
            bases.append(b)
            raw.extend(rs)
            if i % 500 == 0:
                print(f"{i}/{len(tasks)}", flush=True)
    r = pl.DataFrame(raw, infer_schema_length=None).with_columns(pl.col(GEOM).cast(pl.Float64))
    r.write_parquet(out / "geometry-seeds.parquet")
    b = pl.DataFrame(bases, infer_schema_length=None).with_columns(pl.col(TOPO).cast(pl.Float64))
    features = b.join(
        r.group_by("compound_id").agg(
            *[pl.col(k).mean().alias(k + "_mean") for k in GEOM],
            *[pl.col(k).std(ddof=0).alias(k + "_std") for k in GEOM],
        ),
        on="compound_id",
        validate="1:1",
    ).sort("compound_id")
    features.write_parquet(out / "geometry-features.parquet")
    report = dict(
        config=cfg,
        rdkit=rdBase.rdkitVersion,
        n_compounds=features.height,
        status=r.group_by("status").len().to_dicts(),
        errors=r.filter(pl.col("status") != "ok").group_by("error").len().to_dicts(),
        topology_columns=TOPO,
        mean_columns=[k + "_mean" for k in GEOM],
        std_columns=[k + "_std" for k in GEOM],
        method=(
            "ETKDGv3, explicit hydrogens, no minimization; stored SMILES formal charge, "
            "no pH microstate generation; basic-N SMARTS and Gasteiger are proxies, not pKa; "
            "inverse crowding is not SASA"
        ),
    )
    (out / "feature-report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
