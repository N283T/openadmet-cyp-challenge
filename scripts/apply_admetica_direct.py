"""Apply unchanged Admetica CYP classifiers and compare scores with challenge pIC50."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers, nn
from chemprop.data import BatchMolGraph
from rdkit import Chem
from scipy.stats import pearsonr, spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def read(name):
    return pl.read_parquet(next(inp.glob("*-" + name + ".parquet")))


allowed = [
    nn.MeanAggregation,
    nn.BondMessagePassing,
    nn.BinaryClassificationFFN,
    nn.Activation,
    (dict, "lightning.fabric.utilities.data.AttributeDict"),
    (nn.BinaryAccuracy, "chemprop.nn.metrics.BinaryAccuracyMetric"),
    (nn.BinaryAUROC, "chemprop.nn.metrics.BinaryAUROCMetric"),
    (np._core.multiarray.scalar, "numpy.core.multiarray.scalar"),
    (np._core.multiarray._reconstruct, "numpy.core.multiarray._reconstruct"),
    np.ndarray,
    np.dtype,
    np.dtypes.Float64DType,
]
torch.set_num_threads(4)
torch.backends.cuda.matmul.allow_tf32 = False
compounds = read("compounds").sort("compound_id")
labels = read("measurements")
a = read("parent-assignments")
reference = read("aligned-oof")
test_mapping = read("test_mapping")
f = featurizers.SimpleMoleculeMolGraphFeaturizer()
graphs = [f(Chem.MolFromSmiles(s)) for s in compounds["standardized_smiles"]]
all_predictions, metrics, checks = [], [], []
for path in sorted(inp.glob("*.ckpt")):
    enzyme = path.name.split("-", 1)[1].split("-inhibitor")[0].upper()
    ep = enzyme + "_pIC50_direct_inhibition"
    with torch.serialization.safe_globals(allowed):
        ckpt = torch.load(path, map_location="cpu", weights_only=True)
    hp = ckpt["hyper_parameters"]
    print(enzyme, repr(hp), flush=True)
    mp = nn.BondMessagePassing(**{k: v for k, v in hp["message_passing"].items() if k != "cls"})
    predictor = nn.BinaryClassificationFFN(
        **{k: v for k, v in hp["predictor"].items() if k != "cls"}
    )
    agg = nn.MeanAggregation()
    bn = torch.nn.BatchNorm1d(mp.output_dim) if hp["batch_norm"] else torch.nn.Identity()
    for prefix, module in [("message_passing", mp), ("predictor", predictor), ("bn", bn)]:
        state = {
            k.removeprefix(prefix + "."): v
            for k, v in ckpt["state_dict"].items()
            if k.startswith(prefix + ".")
        }
        module.load_state_dict(state, strict=True)
        module.cuda().eval()
    ps, ls = [], []
    with torch.inference_mode():
        for start in range(0, len(graphs), 128):
            graph = BatchMolGraph(graphs[start : start + 128])
            graph.to("cuda")
            emb = bn(agg(mp(graph), graph.batch))
            logits = predictor.train_step(emb)
            probs = predictor(emb)
            torch.testing.assert_close(probs, torch.sigmoid(logits))
            ps.extend(probs.cpu().numpy().ravel().tolist())
            ls.extend(logits.cpu().numpy().ravel().tolist())
    pred = compounds.select("compound_id").with_columns(
        pl.lit(ep).alias("endpoint"), pl.Series("probability", ps), pl.Series("logit", ls)
    )
    assert pred["probability"].is_between(0, 1).all() and pred["logit"].is_finite().all()
    joined = (
        labels.filter(pl.col("endpoint") == ep)
        .join(pred, on=["compound_id", "endpoint"], validate="1:1")
        .join(a, on="compound_id", validate="1:1")
    )
    for subset in ["all_labeled", "matched_oof", "neighbors"]:
        h = joined if subset == "all_labeled" else joined.filter(~pl.col("is_anchor"))
        if subset == "neighbors":
            h = h.filter(pl.col("is_neighbor"))
        for score in ["probability", "logit"]:
            metrics.append(
                dict(
                    endpoint=ep,
                    subset=subset,
                    score=score,
                    n=h.height,
                    pearson=float(pearsonr(h["value"], h[score]).statistic),
                    spearman=float(spearmanr(h["value"], h[score]).statistic),
                )
            )
    matched = reference.filter(pl.col("endpoint") == ep).join(
        pred, on=["compound_id", "endpoint"], validate="1:1"
    )
    for name in ["fp_mordred_lgbm", "emb_tabpfn", "pred_log2fc_tabpfn"]:
        metrics.append(
            dict(
                endpoint=ep,
                subset="correlation_with_" + name,
                score="probability",
                n=matched.height,
                pearson=float(pearsonr(matched[name], matched["probability"]).statistic),
                spearman=float(spearmanr(matched[name], matched["probability"]).statistic),
            )
        )
    checks.append(
        dict(
            endpoint=ep,
            batch_norm=hp["batch_norm"],
            rows=pred.height,
            probability_min=min(ps),
            probability_max=max(ps),
            strict_encoder_and_head_state=True,
        )
    )
    all_predictions.append(pred)
pred = pl.concat(all_predictions)
pred = pred.join(
    a.select("compound_id", "fold"), on="compound_id", how="left", validate="m:1"
).with_columns(pl.col("fold").fill_null(-1), pl.lit("external_pretrained_direct").alias("split_id"))
pred.write_parquet(out / "direct-predictions.parquet")
pred.join(
    test_mapping.select(
        "compound_id",
        pl.col("molecule_name").alias("Molecule_Name"),
        pl.col("raw_smiles").alias("SMILES"),
    ),
    on="compound_id",
    how="inner",
    validate="m:1",
).write_parquet(out / "test-probabilities.parquet")
pl.DataFrame(metrics).write_csv(out / "metrics.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            checks=checks,
            training="none",
            prediction="unchanged pretrained binary inhibition probability and pre-sigmoid logit",
            evaluation="association with continuous pIC50; no probability-to-pIC50 conversion",
        ),
        indent=2,
    )
)
print(pl.DataFrame(metrics).filter(pl.col("subset") == "matched_oof"))
