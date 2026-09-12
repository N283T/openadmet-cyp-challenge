"""Optional medium-sized architecture using the existing external assay trainer.

Architecture comparison only: the external sources, tasks, losses and training-only
scaling remain this project's protocol, not a reproduction of Jeremy's dataset.
"""

import json
import os
from pathlib import Path

import torch
import train_external_assay as trainer
from chemprop import nn as cnn


class MediumAssayModel(torch.nn.Module):
    def __init__(self, n_assays, dropout):
        super().__init__()
        cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
        assert cfg.get("architecture") == "medium"
        assert cfg.get("width", 600) == 600 and cfg.get("depth", 4) == 4
        assert cfg.get("aggregation", "norm") == "norm"
        assert cfg.get("aggregation_norm", 100.0) == 100.0
        assert cfg.get("ffn_hidden_layers", 3) == 3
        assert not cfg.get("conditioning_group")
        self.mp = cnn.BondMessagePassing(d_h=600, depth=4)
        self.aggregation = cnn.NormAggregation(norm=100.0)
        layers = []
        for _ in range(3):
            layers.extend([torch.nn.Linear(600, 600), torch.nn.ReLU(), torch.nn.Dropout(dropout)])
        layers.append(torch.nn.Linear(600, n_assays * 2))
        self.head = torch.nn.Sequential(*layers)
        assert len([m for m in self.head if isinstance(m, torch.nn.Linear)]) == 4
        assert self.mp.output_dim == 600
        metadata = dict(
            architecture="medium",
            width=600,
            depth=4,
            message_passing="BondMessagePassing",
            aggregation="NormAggregation",
            aggregation_norm=100.0,
            ffn_hidden_layers=3,
            ffn_hidden_width=600,
            activation="ReLU",
            dropout=dropout,
            assay_outputs=n_assays,
            output_dimensions=n_assays * 2,
            output_order="assay-major potency then binary logit",
            parameter_count=sum(p.numel() for p in self.parameters()),
            state_key_groups=["mp", "head"],
            protocol="Existing external assay objectives, split, scaler and sampler unchanged",
            scope="Architecture comparison only; no claim to reproduce Jeremy data or training",
        )
        (Path(os.environ["CYP_OUTPUT_DIR"]) / "architecture.json").write_text(
            json.dumps(metadata, indent=2)
        )

    def embed(self, graph):
        return self.aggregation(self.mp(graph), graph.batch)

    def forward(self, graph):
        return self.head(self.embed(graph))


if __name__ == "__main__":
    trainer.AssayModel = MediumAssayModel
    trainer.main()
