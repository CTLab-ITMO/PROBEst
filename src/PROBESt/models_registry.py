import torch
import torch.nn as nn
from torch.nn import functional as F

# ---------- Shallow small net ----------
class ShallowNet(nn.Module):
    def __init__(self, input_size):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, 32),
            nn.ReLU(),
            nn.Linear(32, 1),
            nn.Sigmoid()
        )
    def forward(self, x):
        return self.net(x)

# ---------- Wide fully-connected net ----------
class WideNet(nn.Module):
    def __init__(self, input_size):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, 256),
            nn.ReLU(),
            nn.Linear(256, 256),
            nn.ReLU(),
            nn.Linear(256, 1),
            nn.Sigmoid()
        )
    def forward(self, x):
        return self.net(x)

# ---------- Residual MLP ----------
class ResidualBlock(nn.Module):
    def __init__(self, width):
        super().__init__()
        self.fc = nn.Linear(width, width)
        self.bn = nn.BatchNorm1d(width)
    def forward(self, x):
        return F.relu(self.bn(self.fc(x)) + x)

class ResidualNet(nn.Module):
    def __init__(self, input_size):
        super().__init__()
        self.input_layer = nn.Linear(input_size, 128)
        self.block1 = ResidualBlock(128)
        self.block2 = ResidualBlock(128)
        self.output = nn.Linear(128, 1)
    def forward(self, x):
        x = F.relu(self.input_layer(x))
        x = self.block1(x)
        x = self.block2(x)
        return torch.sigmoid(self.output(x))

# ---------- GAIL-style discriminator ----------
class GAILDiscriminator(nn.Module):
    def __init__(self, input_size):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, 256),
            nn.LeakyReLU(0.2),
            nn.Linear(256, 128),
            nn.LeakyReLU(0.2),
            nn.Linear(128, 1),
            nn.Sigmoid() 
        )
    def forward(self, x):
        return self.net(x)

# ---------- TabTransformer ----------
class TabTransformer(nn.Module):
    def __init__(self, input_size, n_heads=4, depth=3):
        super().__init__()
        self.embedding = nn.Linear(input_size, 64)
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=64, nhead=n_heads, batch_first=True
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=depth)
        self.fc_out = nn.Linear(64, 1)

    def forward(self, x):
        x = self.embedding(x).unsqueeze(1)  
        x = self.transformer(x)
        return torch.sigmoid(self.fc_out(x[:, 0]))
