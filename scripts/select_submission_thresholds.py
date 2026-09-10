"""Choose final binary thresholds on labeled OOF only; never use blinded outcomes."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.metrics import matthews_corrcoef

from cyp_challenge.storage import write_json
from cyp_challenge.tracks import ENDPOINTS, validate


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    oof = pl.read_parquet(next(inp.glob("*-oof.parquet")))
    test = pl.read_parquet(next(inp.glob("*-ensemble_probabilities.parquet")))
    scores = []
    selected = {}
    for ep in ENDPOINTS["classification"]:
        data = oof.filter(pl.col("endpoint") == ep)
        y, p = data["y_true"].to_numpy(), data["prediction"].to_numpy()
        assert data["compound_id"].n_unique() == data.height
        for threshold in np.unique(np.append(p, 0.5)):
            score = float(matthews_corrcoef(y, p >= threshold))
            scores.append(dict(endpoint=ep, threshold=float(threshold), selection_mcc=score))
        options = [r for r in scores if r["endpoint"] == ep]
        # Predeclared tie break: closest to 0.5, then higher threshold.
        selected[ep] = max(
            options, key=lambda r: (r["selection_mcc"], -abs(r["threshold"] - 0.5), r["threshold"])
        )
    pl.DataFrame(scores).write_parquet(out / "threshold_search.parquet")
    wide = (
        test.filter(pl.col("endpoint").is_in(ENDPOINTS["classification"]))
        .pivot(on="endpoint", index=["compound_id", "SMILES", "Molecule_Name"], values="prediction")
        .sort("Molecule_Name")
    )
    result = wide.select(
        "SMILES",
        "Molecule_Name",
        *[
            (pl.col(ep) >= selected[ep]["threshold"]).alias(ep)
            for ep in ENDPOINTS["classification"]
        ],
    )
    result.write_parquet(out / "classification.parquet")
    validation = validate(
        out / "classification.parquet",
        "classification",
        next(inp.glob("*-cyp-challenge-TEST-BLINDED.csv")),
    )
    summary = dict(
        config=cfg,
        thresholds=selected,
        positive_counts={ep: int(result[ep].sum()) for ep in ENDPOINTS["classification"]},
        validation=validation,
        warning=(
            "Selection MCC is fitted on the same OOF labels; not an unbiased CV score. "
            "Deployed probabilities average five models, while selection OOF has "
            "one held-out model per row."
        ),
    )
    write_json(out / "summary.json", summary)
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
