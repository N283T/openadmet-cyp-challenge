"""Transparent RDKit approximations of literature TDI bioactivation alerts.

These SMARTS encode the categories named in the OpenADMET structural-alerts
article. They are an independent, version-controlled implementation: the exact
SMARTS used for the article were not available in its archived source.
"""

from dataclasses import dataclass

import numpy as np
from rdkit import Chem


@dataclass(frozen=True)
class TDIAlert:
    name: str
    smarts: str


TDI_ALERTS = (
    TDIAlert("terminal_alkyne", "[C;H]#[C]"),
    TDIAlert("terminal_alkene", "[C;H2]=[C]"),
    TDIAlert("furan", "[o;r5]1[c;r5][c;r5][c;r5][c;r5]1"),
    TDIAlert("thiophene", "[s;r5]1[c;r5][c;r5][c;r5][c;r5]1"),
    TDIAlert("methylenedioxyphenyl", "c1ccc2OCOc2c1"),
    TDIAlert(
        "primary_secondary_amine",
        "[N;H1,H2;+0;!$(N-[C,S,P]=[O,S,N])]",
    ),
    TDIAlert(
        "hydrazine",
        "[N;!$(N-[C,S,P]=[O,S,N])]-[N;!$(N-[C,S,P]=[O,S,N])]",
    ),
    TDIAlert("aniline_aminophenol", "[N;H1,H2;+0;!$(N-C=O)]-[c]"),
    TDIAlert("thioamide", "[C;X3](=[S])-[N]"),
    TDIAlert("thiourea", "[N]-[C;X3](=[S])-[N]"),
    TDIAlert("cyclopropylamine", "[N]-[C;r3]"),
    TDIAlert(
        "n_dealkylatable_tertiary_amine",
        "[N;X3;H0;+0;!$(N-[C,S,P]=[O,S,N])]-[C;H1,H2,H3]",
    ),
)


def compiled_alerts() -> tuple[tuple[str, Chem.Mol], ...]:
    """Compile and validate the fixed alert panel."""
    result = tuple((alert.name, Chem.MolFromSmarts(alert.smarts)) for alert in TDI_ALERTS)
    invalid = [name for name, query in result if query is None]
    if invalid:
        raise ValueError(f"Invalid TDI SMARTS: {invalid}")
    return result


def alert_matrix(smiles: list[str]) -> tuple[np.ndarray, list[str]]:
    """Return per-molecule match counts plus ``any`` and total-count features."""
    alerts = compiled_alerts()
    molecules = [Chem.MolFromSmiles(value) for value in smiles]
    invalid = [index for index, molecule in enumerate(molecules) if molecule is None]
    if invalid:
        raise ValueError(f"Invalid standardized SMILES at rows: {invalid[:10]}")

    counts = np.asarray(
        [
            [len(molecule.GetSubstructMatches(query, uniquify=True)) for _, query in alerts]
            for molecule in molecules
        ],
        dtype=np.float32,
    )
    any_alert = (counts.sum(axis=1, keepdims=True) > 0).astype(np.float32)
    total = counts.sum(axis=1, keepdims=True)
    names = [name for name, _ in alerts] + ["any_tdi_alert", "tdi_alert_count"]
    return np.column_stack([counts, any_alert, total]), names
