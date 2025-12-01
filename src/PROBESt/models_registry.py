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
    def __init__(self, input_size, hidden1=256, hidden2=128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, hidden1),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden1, hidden2),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden2, 1),
            nn.Sigmoid() 
        )
    def forward(self, x):
        return self.net(x)

# ---------- GAIL variations for architecture search ----------
class GAILDeep(nn.Module):
    """Deeper GAIL with 3 hidden layers"""
    def __init__(self, input_size, hidden1=256, hidden2=128, hidden3=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, hidden1),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden1, hidden2),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden2, hidden3),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden3, 1),
            nn.Sigmoid()
        )
    def forward(self, x):
        return self.net(x)

class GAILWide(nn.Module):
    """Wider GAIL with larger hidden layers"""
    def __init__(self, input_size, hidden1=512, hidden2=256):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, hidden1),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden1, hidden2),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden2, 1),
            nn.Sigmoid()
        )
    def forward(self, x):
        return self.net(x)

class GAILNarrow(nn.Module):
    """Narrower GAIL with smaller hidden layers"""
    def __init__(self, input_size, hidden1=128, hidden2=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, hidden1),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden1, hidden2),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden2, 1),
            nn.Sigmoid()
        )
    def forward(self, x):
        return self.net(x)

class GAILWithDropout(nn.Module):
    """GAIL with dropout for regularization"""
    def __init__(self, input_size, hidden1=256, hidden2=128, dropout=0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, hidden1),
            nn.LeakyReLU(0.2),
            nn.Dropout(dropout),
            nn.Linear(hidden1, hidden2),
            nn.LeakyReLU(0.2),
            nn.Dropout(dropout),
            nn.Linear(hidden2, 1),
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
