"""Extract atom-level SMARTCyp predictions and aggregate TDI-oriented features."""

import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem

from cyp_challenge.tdi_alerts import compiled_alerts

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
scratch = Path(os.environ["TMPDIR"]) / "smartcyp-features"
scratch.mkdir(parents=True, exist_ok=False)
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def input_path(name: str) -> Path:
    return inp / cfg[name]


jar = input_path("jar")
jar_digest = hashlib.sha256(jar.read_bytes()).hexdigest()
if jar_digest != cfg["jar_sha256"]:
    raise ValueError(f"SMARTCyp JAR hash mismatch: {jar_digest}")

compounds = pl.read_parquet(input_path("compounds")).sort("compound_id")
if compounds["compound_id"].n_unique() != compounds.height:
    raise ValueError("Compound IDs are not unique")

molecules: list[Chem.Mol] = []
for record in compounds.iter_rows(named=True):
    molecule = Chem.MolFromSmiles(record["standardized_smiles"])
    if molecule is None:
        raise ValueError(f"Invalid standardized SMILES for {record['compound_id']}")
    molecules.append(molecule)

compound_ids = compounds["compound_id"].to_list()
smiles = compounds["standardized_smiles"].to_list()
raw_records = []
failures = []
log_lines = []


def run_batch(indices: list[int], tag: str) -> None:
    """Run one batch, recursively isolating structures rejected by the old CDK."""
    batch_dir = scratch / tag
    batch_dir.mkdir()
    sdf_path = batch_dir / "compounds.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    for local_row, global_index in enumerate(indices, start=1):
        molecule = Chem.Mol(molecules[global_index])
        molecule.SetProp("_Name", str(compound_ids[global_index]))
        molecule.SetIntProp("SMARTCyp_input_row", local_row)
        writer.write(molecule)
    writer.close()
    command = ["java", "-jar", str(jar), "-printall", str(sdf_path)]
    result = subprocess.run(
        command,
        cwd=batch_dir,
        text=True,
        capture_output=True,
        timeout=cfg["timeout_seconds"],
        check=False,
    )
    log_lines.append(
        f"{tag}\t{len(indices)}\t{returncode_label(result.returncode)}\t"
        f"{result.stdout.strip()}\t{result.stderr.strip()}"
    )
    csv_files = sorted(batch_dir.glob("SMARTCyp_Results_*.csv"))
    if result.returncode == 0 and len(csv_files) == 1:
        batch = pl.read_csv(csv_files[0], null_values="null")
        for record in batch.iter_rows(named=True):
            local_index = int(record["Molecule"]) - 1
            record["_global_index"] = indices[local_index]
            raw_records.append(record)
        return
    if len(indices) > 1:
        middle = len(indices) // 2
        run_batch(indices[:middle], f"{tag}a")
        run_batch(indices[middle:], f"{tag}b")
        return
    global_index = indices[0]
    failures.append(
        {
            "compound_id": compound_ids[global_index],
            "standardized_smiles": smiles[global_index],
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
    )


def returncode_label(returncode: int) -> str:
    return "ok" if returncode == 0 else f"failed:{returncode}"


chunk_size = cfg["chunk_size"]
for start in range(0, compounds.height, chunk_size):
    stop = min(start + chunk_size, compounds.height)
    run_batch(list(range(start, stop)), f"batch-{start:05d}-{stop:05d}")
    print(f"SMARTCyp attempted {stop}/{compounds.height}", flush=True)

(out / "smartcyp-batch-log.txt").write_text("\n".join(log_lines) + "\n")
failure_table = pl.DataFrame(failures) if failures else pl.DataFrame(
    schema={
        "compound_id": pl.String,
        "standardized_smiles": pl.String,
        "returncode": pl.Int64,
        "stdout": pl.String,
        "stderr": pl.String,
    }
)
failure_table.write_parquet(out / "smartcyp-failures.parquet")
raw = pl.DataFrame(raw_records)
raw.write_csv(out / "smartcyp-raw.csv")
expected_columns = {
    "Molecule",
    "Atom",
    "Ranking",
    "Score",
    "Energy",
    "Relative Span",
    "2D6ranking",
    "2D6score",
    "Span2End",
    "N+Dist",
    "2DSASA",
}
if not expected_columns.issubset(raw.columns):
    raise ValueError(f"Unexpected SMARTCyp columns: {raw.columns}")
alert_queries = compiled_alerts()
alert_atom_sets: list[dict[str, set[int]]] = []
for molecule in molecules:
    alert_atom_sets.append(
        {
            name: {index for match in molecule.GetSubstructMatches(query) for index in match}
            for name, query in alert_queries
        }
    )

atom_rows = []
for record in raw.iter_rows(named=True):
    molecule_index = int(record["_global_index"])
    atom_index = int(str(record["Atom"]).rsplit(".", 1)[1]) - 1
    molecule = molecules[molecule_index]
    if not 0 <= atom_index < molecule.GetNumAtoms():
        raise ValueError(f"Invalid atom index in SMARTCyp output: {record['Atom']}")
    matches = alert_atom_sets[molecule_index]
    hit_names = [name for name, indices in matches.items() if atom_index in indices]
    atom_rows.append(
        {
            "compound_id": compound_ids[molecule_index],
            "standardized_smiles": smiles[molecule_index],
            "atom_index": atom_index,
            "atom_symbol": molecule.GetAtomWithIdx(atom_index).GetSymbol(),
            "ranking_3a4": record["Ranking"],
            "score_3a4": record["Score"],
            "ranking_2d6": record["2D6ranking"],
            "score_2d6": record["2D6score"],
            "energy": record["Energy"],
            "relative_span": record["Relative Span"],
            "span_to_end": record["Span2End"],
            "distance_to_protonated_n": record["N+Dist"],
            "sasa_2d": record["2DSASA"],
            "tdi_alert_atom": bool(hit_names),
            "tdi_alert_names": "|".join(hit_names),
        }
    )
atoms = pl.DataFrame(atom_rows).sort("compound_id", "atom_index")
atoms.write_parquet(out / "smartcyp-atoms.parquet")


def finite(values: list[float | None], ceiling: float = 900.0) -> np.ndarray:
    array = np.asarray([value for value in values if value is not None], dtype=float)
    return array[np.isfinite(array) & (array < ceiling)]


feature_rows = []
for compound_id in compound_ids:
    table = atoms.filter(pl.col("compound_id") == compound_id)
    if table.is_empty():
        continue
    endpoint_values = {}
    row = {
        "compound_id": compound_id,
        "smartcyp_heavy_atoms": table.height,
        "smartcyp_reactive_atoms": int((table["energy"] < 900).sum()),
    }
    for isoform, score_column, rank_column in [
        ("3a4", "score_3a4", "ranking_3a4"),
        ("2d6", "score_2d6", "ranking_2d6"),
    ]:
        scores = finite(table[score_column].to_list())
        if scores.size == 0:
            raise ValueError(f"No reactive SMARTCyp atoms for {compound_id} {isoform}")
        ordered = np.sort(scores)
        for rank in range(3):
            row[f"smartcyp_{isoform}_score_top{rank + 1}"] = float(
                ordered[min(rank, len(ordered) - 1)]
            )
        row[f"smartcyp_{isoform}_top_gap"] = float(
            ordered[min(1, len(ordered) - 1)] - ordered[0]
        )
        row[f"smartcyp_{isoform}_top3_mean"] = float(ordered[:3].mean())
        row[f"smartcyp_{isoform}_score_q25"] = float(np.quantile(scores, 0.25))
        row[f"smartcyp_{isoform}_score_mean"] = float(scores.mean())
        top = table.filter(pl.col(rank_column) == 1)
        if top.height != 1:
            top = table.sort(score_column).head(1)
        top_record = top.row(0, named=True)
        top_atom = int(top_record["atom_index"])
        endpoint_values[isoform] = (float(ordered[0]), top_atom)
        row[f"smartcyp_{isoform}_top_atom"] = top_atom
        row[f"smartcyp_{isoform}_top_energy"] = float(top_record["energy"])
        row[f"smartcyp_{isoform}_top_relative_span"] = float(top_record["relative_span"])
        row[f"smartcyp_{isoform}_top_sasa_2d"] = float(top_record["sasa_2d"])
        row[f"smartcyp_{isoform}_top_is_tdi_alert"] = int(top_record["tdi_alert_atom"])
        ranked = table.sort(score_column).head(3)
        row[f"smartcyp_{isoform}_top3_tdi_alert_count"] = int(
            ranked["tdi_alert_atom"].sum()
        )
    row["smartcyp_min_score_2d6_minus_3a4"] = (
        endpoint_values["2d6"][0] - endpoint_values["3a4"][0]
    )
    row["smartcyp_same_top_atom"] = int(
        endpoint_values["2d6"][1] == endpoint_values["3a4"][1]
    )
    feature_rows.append(row)

successful_features = pl.DataFrame(feature_rows).sort("compound_id")
features = compounds.select("compound_id").join(
    successful_features, on="compound_id", how="left", validate="1:1"
)
if features.height != compounds.height:
    raise ValueError("SMARTCyp aggregate feature row coverage is incomplete")
features.write_parquet(out / "smartcyp-features.parquet")

summary = {
    "source_commit": cfg["source_commit"],
    "jar_sha256": jar_digest,
    "compound_count": compounds.height,
    "successful_compounds": successful_features.height,
    "failed_compounds": failure_table.height,
    "atom_rows": atoms.height,
    "feature_count": len(features.columns) - 1,
    "feature_columns": features.columns[1:],
    "scope": (
        "SMARTCyp generic/3A4 and 2D6 atom scores plus top-SoM interactions "
        "with the fixed TDI alert panel"
    ),
    "limitations": [
        "SMARTCyp 3A4 is the generic reactivity/accessibility score; "
        "2D6 has an isoform correction.",
        "This predicts metabolic sites, not TDI or reactive-metabolite formation directly.",
        "Alert interactions use the existing transparent approximate SMARTS panel.",
    ],
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
