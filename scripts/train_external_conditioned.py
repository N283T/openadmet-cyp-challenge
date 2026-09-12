"""External Chemprop600 with shared CYP/measurement outputs and assay calibration.

Reuses the assay trainer's data split, curation, objective/source-balanced sampler,
losses, checkpoint selection, and export. Replaces independent assay MLP outputs
with shared biological outputs and small learned assay-specific affine terms.
"""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
import train_external_assay as trainer
from chemprop import nn as cnn


class ConditionedAssayModel(torch.nn.Module):
    def __init__(self, n_assays, dropout):
        super().__init__()
        cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
        out = Path(os.environ["CYP_OUTPUT_DIR"])
        observations = trainer.read_input(Path(os.environ["CYP_INPUT_DIR"]), "observations.parquet")
        tasks = json.loads((out / "tasks.json").read_text())
        assays = sorted({task["assay_id"] for task in tasks})
        assert len(assays) == n_assays
        grouping = cfg.get("conditioning_group", "cyp_measurement")
        if grouping not in {"cyp", "cyp_measurement"}:
            raise ValueError(f"Unknown conditioning_group: {grouping}")
        mode = cfg.get("conditioning_variant", "affine")
        if mode not in {"affine", "offset"}:
            raise ValueError(f"Unknown conditioning_variant: {mode}")
        metadata = (
            observations.filter(pl.col("assay_id").is_in(assays))
            .select("assay_id", "isoform", "measurement_type")
            .unique()
        )
        assert metadata["assay_id"].n_unique() == metadata.height
        assay_group = {
            row["assay_id"]: (
                row["isoform"]
                if grouping == "cyp"
                else row["isoform"] + "/" + row["measurement_type"]
            )
            for row in metadata.iter_rows(named=True)
        }
        groups = sorted(set(assay_group.values()))
        lookup = {key: i for i, key in enumerate(groups)}
        group_indices = torch.tensor([lookup[assay_group[assay]] for assay in assays])
        self.register_buffer("group_indices", group_indices)
        self.mp = cnn.BondMessagePassing(d_h=600, depth=5)
        self.aggregation = cnn.MeanAggregation()
        self.shared_head = torch.nn.Sequential(
            torch.nn.Linear(600, 600),
            torch.nn.ReLU(),
            torch.nn.Dropout(dropout),
            torch.nn.Linear(600, len(groups) * 2),
        )
        self.offset = torch.nn.Parameter(torch.zeros(n_assays, 2))
        self.log_scale = torch.nn.Parameter(
            torch.zeros(n_assays, 2), requires_grad=mode == "affine"
        )
        (out / "conditioning.json").write_text(
            json.dumps(
                dict(
                    grouping=grouping,
                    variant=mode,
                    groups=groups,
                    assay_groups={assay: assay_group[assay] for assay in assays},
                    assays=assays,
                    output_order="assay-major: potency then binary logit",
                    positive_scale="exp(2*tanh(log_scale)); initialized to 1",
                    scale_range=[float(np.exp(-2)), float(np.exp(2))],
                    initial_offset=0.0,
                    parameter_count=sum(p.numel() for p in self.parameters()),
                    trainable_parameters=sum(
                        p.numel() for p in self.parameters() if p.requires_grad
                    ),
                ),
                indent=2,
            )
        )

    def embed(self, graph):
        return self.aggregation(self.mp(graph), graph.batch)

    def forward(self, graph):
        shared = self.shared_head(self.embed(graph)).reshape(
            -1, self.shared_head[-1].out_features // 2, 2
        )
        expanded = shared[:, self.group_indices]
        calibrated = expanded * torch.exp(2 * torch.tanh(self.log_scale)) + self.offset
        return calibrated.flatten(1)


if __name__ == "__main__":
    trainer.AssayModel = ConditionedAssayModel
    trainer.main()
