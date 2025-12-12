# MIT License
#
# Copyright (c) 2025 CTLab-ITMO
#
# Authors: Daniil Smutin, Aleksandr Serdiukov, Vitalii Dravgelis, Artem Ivanov,
# Aleksei Zabashta, Sergey Muravyov, and the CTLab-ITMO university team.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


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
        # Use LayerNorm instead of BatchNorm to handle batch size 1
        self.ln = nn.LayerNorm(width)
    def forward(self, x):
        return F.relu(self.ln(self.fc(x)) + x)

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

# ---------- GAIL_Wide variations for further architecture search ----------
class GAILWideDeep(nn.Module):
    """Wide and deep GAIL with 3 layers"""
    def __init__(self, input_size, hidden1=512, hidden2=256, hidden3=128):
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

class GAILWideDropout(nn.Module):
    """Wide GAIL with dropout regularization"""
    def __init__(self, input_size, hidden1=512, hidden2=256, dropout=0.3):
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

class GAILWideBatchNorm(nn.Module):
    """Wide GAIL with layer normalization (works with any batch size)"""
    def __init__(self, input_size, hidden1=512, hidden2=256):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, hidden1),
            nn.LayerNorm(hidden1),  # Use LayerNorm instead of BatchNorm for batch size 1 compatibility
            nn.LeakyReLU(0.2),
            nn.Linear(hidden1, hidden2),
            nn.LayerNorm(hidden2),  # Use LayerNorm instead of BatchNorm for batch size 1 compatibility
            nn.LeakyReLU(0.2),
            nn.Linear(hidden2, 1),
            nn.Sigmoid()
        )
    def forward(self, x):
        return self.net(x)

class GAILWideExtra(nn.Module):
    """Extra wide GAIL with even larger layers"""
    def __init__(self, input_size, hidden1=768, hidden2=384):
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

class GAILWideBalanced(nn.Module):
    """Wide GAIL with balanced layer sizes"""
    def __init__(self, input_size, hidden1=512, hidden2=256, hidden3=128):
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
