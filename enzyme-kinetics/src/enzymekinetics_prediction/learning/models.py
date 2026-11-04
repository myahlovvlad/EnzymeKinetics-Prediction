"""CPU neural regressors; radial attention is not a molecular interaction GNN."""

import torch
from torch import nn


class KineticsMLP(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.regressor = nn.Sequential(
            nn.Linear(input_dim, 128),
            nn.SiLU(),
            nn.Dropout(0.15),
            nn.Linear(128, 64),
            nn.SiLU(),
            nn.Linear(64, 1),
        )

    def forward(self, features, geometry):
        return self.regressor(features).squeeze(-1)


class PhysicsAdd3DGraphNet(nn.Module):
    """Supervised adaptation of the original radial softmax-pooling prototype.

    Uses 32 CA radial quantiles rather than variable-length atoms. The output is
    log10 kinetics, not an uncalibrated docking probability or a physical energy.
    No edges, atom types, charges, QM barriers or binding pose are inferred here.
    """

    def __init__(self, input_dim):
        super().__init__()
        self.node_encoder = nn.Sequential(nn.Linear(4, 32), nn.Tanh())
        self.node_score = nn.Linear(32, 1)
        self.regressor = nn.Sequential(
            nn.Linear(input_dim + 32, 128),
            nn.SiLU(),
            nn.Dropout(0.15),
            nn.Linear(128, 64),
            nn.SiLU(),
            nn.Linear(64, 1),
        )

    def forward(self, features, geometry):
        r = geometry.unsqueeze(-1)
        node_features = torch.cat(
            [r, r.square(), torch.exp(-r), torch.ones_like(r)], dim=-1
        )
        hidden = self.node_encoder(node_features)
        weights = torch.softmax(-self.node_score(hidden), dim=1)
        pooled = (weights * hidden).sum(dim=1)
        return self.regressor(torch.cat([features, pooled], dim=1)).squeeze(-1)
