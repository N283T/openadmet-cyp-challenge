"""How alike are the OpenADMET CheMeleon checkpoints at one message-passing depth?

Two views for one regression endpoint. Linear CKA between the pooled embeddings of every checkpoint
pair: the 2048 dimensions of different networks do not correspond, so a per-dimension correlation
would be meaningless. Pearson correlation between the checkpoints' OOF readout predictions and
residuals on the incumbent's cohort, with the incumbent included on the prediction side.
"""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from regression_depth_super_ensemble import KEYS, collect

READOUTS = ("lgbm", "tabpfn")


def load_embedding(path):
    frame = pl.read_parquet(path).sort("compound_id")
    columns = [c for c in frame.columns if c.startswith("emb_")]
    x = frame.select(columns).to_numpy().astype(np.float32)
    return frame["compound_id"].to_list(), x - x.mean(axis=0)


def linear_cka(embeddings):
    names = list(embeddings)
    norms = {n: float(np.linalg.norm(embeddings[n].T @ embeddings[n])) for n in names}
    matrix = np.eye(len(names))
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            b = names[j]
            cross = float(np.linalg.norm(embeddings[a].T @ embeddings[b])) ** 2
            matrix[i, j] = matrix[j, i] = cross / (norms[a] * norms[b])
    return names, matrix


def write_matrix(path, names, matrix):
    frame = pl.DataFrame({"checkpoint": names, **{n: matrix[:, k] for k, n in enumerate(names)}})
    frame.write_csv(path)


def off_diagonal_mean(matrix):
    mask = ~np.eye(len(matrix), dtype=bool)
    return float(matrix[mask].mean())


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    depth, endpoint = cfg["depth"], cfg["endpoint"]

    ids, embeddings = None, {}
    for path in sorted(inp.glob("emb-*.parquet")):
        compound_ids, x = load_embedding(path)
        assert ids is None or compound_ids == ids, "checkpoints cover different compounds"
        ids = compound_ids
        embeddings[path.stem.removeprefix("emb-")] = x
    names, cka = linear_cka(embeddings)
    write_matrix(out / "cka.csv", names, cka)

    control = cfg["control"]
    frame = pl.read_parquet(inp / control["file"])
    for column, value in control["filter"].items():
        frame = frame.filter(pl.col(column) == value)
    table = (
        frame.filter(pl.col("endpoint") == endpoint)
        .sort("compound_id")
        .select(*KEYS, "y_true", pl.col("prediction").alias("incumbent"))
    )
    assert table.height == control["rows"]
    members = collect(Path(cfg["runs_root"]), endpoint)
    for name in names:
        for readout in READOUTS:
            oof = members[f"{name}:{depth}:{readout}"][0]
            part = oof.select(*KEYS, pl.col("prediction").alias(f"{name}:{readout}"))
            table = table.join(part, on=KEYS, how="left", validate="1:1")
    assert table.null_count().sum_horizontal().item() == 0
    y = table["y_true"].to_numpy()
    predictions = {"incumbent": table["incumbent"].to_numpy()}
    for name in names:
        predictions[name] = table.select([f"{name}:{r}" for r in READOUTS]).to_numpy().mean(axis=1)
    labels = list(predictions)
    stacked = np.column_stack([predictions[k] for k in labels])
    residuals = y[:, None] - stacked
    write_matrix(out / "prediction-corr.csv", labels, np.corrcoef(stacked.T))
    residual_corr = np.corrcoef(residuals.T)
    write_matrix(out / "residual-corr.csv", labels, residual_corr)

    within = {
        name: float(
            np.corrcoef(
                y - table[f"{name}:lgbm"].to_numpy(), y - table[f"{name}:tabpfn"].to_numpy()
            )[0, 1]
        )
        for name in names
    }
    pairs = residual_corr[1:, 1:]
    summary = {
        "depth": depth,
        "checkpoints": len(names),
        "cka_mean_off_diagonal": off_diagonal_mean(cka),
        "cka_to_chemeleon": {
            n: float(cka[names.index("chemeleon"), k]) for k, n in enumerate(names)
        },
        "cka_mean_to_others": {
            n: float((cka[k].sum() - 1) / (len(names) - 1)) for k, n in enumerate(names)
        },
        "residual_corr_mean_off_diagonal": off_diagonal_mean(pairs),
        "residual_corr_with_incumbent": {
            n: float(residual_corr[0, 1 + k]) for k, n in enumerate(names)
        },
        "residual_corr_lgbm_vs_tabpfn_within_checkpoint": within,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
