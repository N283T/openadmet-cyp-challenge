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
        assert cfg.get("architecture") == "owned_width_grid"
        width = int(cfg["width"])
        depth = int(cfg.get("depth", 4))
        aggregation = cfg.get("aggregation", "norm")
        aggregation_norm = float(cfg.get("aggregation_norm", 100.0))
        hidden_layers = int(cfg.get("ffn_hidden_layers", 3))
        assert width in {256, 384, 448, 512, 600, 768, 896, 1024}
        assert depth in {2, 3, 4, 5, 6}
        assert aggregation in {"norm", "mean", "sum"}
        assert hidden_layers in {1, 2, 3, 4}
        assert not cfg.get("conditioning_group")
        self.mp = cnn.BondMessagePassing(d_h=width, depth=depth)
        if aggregation == "norm":
            self.aggregation = cnn.NormAggregation(norm=aggregation_norm)
        elif aggregation == "mean":
            self.aggregation = cnn.MeanAggregation()
        else:
            self.aggregation = cnn.SumAggregation()
        layers = []
        for _ in range(hidden_layers):
            layers.extend(
                [torch.nn.Linear(width, width), torch.nn.ReLU(), torch.nn.Dropout(dropout)]
            )
        layers.append(torch.nn.Linear(width, n_assays * 2))
        self.head = torch.nn.Sequential(*layers)
        assert len([m for m in self.head if isinstance(m, torch.nn.Linear)]) == hidden_layers + 1
        assert self.mp.output_dim == width
        metadata = dict(
            architecture="owned_width_grid",
            width=width,
            depth=depth,
            message_passing="BondMessagePassing",
            aggregation=type(self.aggregation).__name__,
            aggregation_norm=aggregation_norm if aggregation == "norm" else None,
            ffn_hidden_layers=hidden_layers,
            ffn_hidden_width=width,
            activation="ReLU",
            dropout=dropout,
            assay_outputs=n_assays,
            output_dimensions=n_assays * 2,
            output_order="assay-major potency then binary logit",
            parameter_count=sum(p.numel() for p in self.parameters()),
            state_key_groups=["mp", "head"],
            protocol="Existing external assay objectives, split, scaler and sampler unchanged",
            scope="Project-owned encoder width/seed grid; no Jeremy checkpoint or training data",
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
