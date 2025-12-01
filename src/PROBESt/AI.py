import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import torch
import torch.nn as nn
from typing import Tuple, Any, Dict

class BaseAIModel:
    def __init__(self):
        self.model = None
        self.scaler = StandardScaler()
        
    def preprocess_data(self, X: pd.DataFrame) -> np.ndarray:
        """Preprocess the input data by scaling numerical features."""
        return self.scaler.fit_transform(X)
    
    def train(self, X: pd.DataFrame, y: pd.Series) -> None:
        """Train the model."""
        raise NotImplementedError
        
    def predict(self, X: pd.DataFrame) -> np.ndarray:
        """Make predictions."""
        raise NotImplementedError

class LogisticRegressionModel(BaseAIModel):
    def __init__(self, **kwargs):
        super().__init__()
        self.model = LogisticRegression(**kwargs)
        
    def train(self, X: pd.DataFrame, y: pd.Series) -> None:
        """Train the logistic regression model."""
        X_scaled = self.preprocess_data(X)
        self.model.fit(X_scaled, y)
        
    def predict(self, X: pd.DataFrame) -> np.ndarray:
        """Make predictions using the trained model."""
        X_scaled = self.scaler.transform(X)
        return self.model.predict(X_scaled)
    
    def predict_proba(self, X: pd.DataFrame) -> np.ndarray:
        """Get probability predictions."""
        X_scaled = self.scaler.transform(X)
        return self.model.predict_proba(X_scaled)

class SimplePerceptron(nn.Module):
    def __init__(self, input_size: int):
        super().__init__()
        self.layer1 = nn.Linear(input_size, 64)
        self.layer2 = nn.Linear(64, 32)
        self.layer3 = nn.Linear(32, 1)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
        
    def forward(self, x):
        x = self.relu(self.layer1(x))
        x = self.relu(self.layer2(x))
        x = self.sigmoid(self.layer3(x))
        return x

class PerceptronModel(BaseAIModel):
    def __init__(self, input_size: int, learning_rate: float = 0.001):
        super().__init__()
        self.model = SimplePerceptron(input_size)
        self.learning_rate = learning_rate
        self.criterion = nn.BCELoss()
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=learning_rate)
        
    def train(self, X: pd.DataFrame, y: pd.Series, epochs: int = 100) -> None:
        """Train the perceptron model."""
        X_scaled = self.preprocess_data(X)
        X_tensor = torch.FloatTensor(X_scaled)
        y_tensor = torch.FloatTensor(y.values).reshape(-1, 1)
        
        for epoch in range(epochs):
            self.optimizer.zero_grad()
            outputs = self.model(X_tensor)
            loss = self.criterion(outputs, y_tensor)
            loss.backward()
            self.optimizer.step()
            
    def predict(self, X: pd.DataFrame) -> np.ndarray:
        """Make predictions using the trained model."""
        X_scaled = self.scaler.transform(X)
        X_tensor = torch.FloatTensor(X_scaled)
        with torch.no_grad():
            predictions = self.model(X_tensor)
        return (predictions.numpy() > 0.5).astype(int)
    
    def predict_proba(self, X: pd.DataFrame) -> np.ndarray:
        """Get probability predictions."""
        X_scaled = self.scaler.transform(X)
        X_tensor = torch.FloatTensor(X_scaled)
        with torch.no_grad():
            predictions = self.model(X_tensor)
        return predictions.numpy()

class DeepNeuralNetwork(nn.Module):
    def __init__(self, input_size: int, dropout_rate: float = 0.3):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(input_size, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            
            nn.Linear(128, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            
            nn.Linear(64, 32),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            
            nn.Linear(32, 16),
            nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            
            nn.Linear(16, 1),
            nn.Sigmoid()
        )
        
    def forward(self, x):
        return self.network(x)
class TorchClassifier(BaseAIModel):
    def __init__(self, model: nn.Module, learning_rate=0.001, weight_pos=1.0):
        super().__init__()
        self.model = model
        self.learning_rate = learning_rate
        
        pos_weight = torch.tensor([weight_pos])
        self.criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        
        self.optimizer = torch.optim.AdamW(self.model.parameters(), lr=learning_rate)
        
    def train(self, X, y, epochs=100, batch_size=32):
        X_scaled = self.preprocess_data(X)
        X_tensor = torch.FloatTensor(X_scaled)
        y_tensor = torch.FloatTensor(y.values).reshape(-1, 1)
        
        dataset = torch.utils.data.TensorDataset(X_tensor, y_tensor)
        loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True)
        
        for e in range(epochs):
            total_loss = 0
            for bx, by in loader:
                self.optimizer.zero_grad()
                logits = self.model(bx)
                loss = self.criterion(logits, by)
                loss.backward()
                self.optimizer.step()
                total_loss += loss.item()
            
            if e % 20 == 0:
                print(f"Epoch {e}: loss = {total_loss:.4f}")

    def predict(self, X):
        X_scaled = self.scaler.transform(X)
        X_tensor = torch.FloatTensor(X_scaled)
        with torch.no_grad():
            preds = torch.sigmoid(self.model(X_tensor)).numpy()
        return (preds > 0.5).astype(int)

    def predict_proba(self, X):
        X_scaled = self.scaler.transform(X)
        X_tensor = torch.FloatTensor(X_scaled)
        with torch.no_grad():
            preds = torch.sigmoid(self.model(X_tensor)).numpy()
        return preds

class DeepNeuralNetworkModel(BaseAIModel):
    def __init__(self, input_size: int, learning_rate: float = 0.001, dropout_rate: float = 0.3):
        super().__init__()
        self.model = DeepNeuralNetwork(input_size, dropout_rate)
        self.learning_rate = learning_rate
        self.criterion = nn.BCELoss()
        self.optimizer = torch.optim.AdamW(
            self.model.parameters(),
            lr=learning_rate,
            weight_decay=0.01
        )
        self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            self.optimizer,
            mode='min',
            factor=0.5,
            patience=5
        )
        
    def train(self, X: pd.DataFrame, y: pd.Series, epochs: int = 200, batch_size: int = 32) -> None:
        """Train the deep neural network model."""
        X_scaled = self.preprocess_data(X)
        X_tensor = torch.FloatTensor(X_scaled)
        y_tensor = torch.FloatTensor(y.values).reshape(-1, 1)
        
        # Create data loader for batch training
        dataset = torch.utils.data.TensorDataset(X_tensor, y_tensor)
        dataloader = torch.utils.data.DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=True
        )
        
        best_loss = float('inf')
        patience = 10
        patience_counter = 0
        
        for epoch in range(epochs):
            self.model.train()
            epoch_loss = 0
            for batch_X, batch_y in dataloader:
                self.optimizer.zero_grad()
                outputs = self.model(batch_X)
                loss = self.criterion(outputs, batch_y)
                loss.backward()
                self.optimizer.step()
                epoch_loss += loss.item()
            
            epoch_loss /= len(dataloader)
            self.scheduler.step(epoch_loss)
            
            # Early stopping
            if epoch_loss < best_loss:
                best_loss = epoch_loss
                patience_counter = 0
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    print(f"Early stopping at epoch {epoch}")
                    break
            
    def predict(self, X: pd.DataFrame) -> np.ndarray:
        """Make predictions using the trained model."""
        X_scaled = self.scaler.transform(X)
        X_tensor = torch.FloatTensor(X_scaled)
        self.model.eval()
        with torch.no_grad():
            predictions = self.model(X_tensor)
        return (predictions.numpy() > 0.5).astype(int)
    
    def predict_proba(self, X: pd.DataFrame) -> np.ndarray:
        """Get probability predictions."""
        X_scaled = self.scaler.transform(X)
        X_tensor = torch.FloatTensor(X_scaled)
        self.model.eval()
        with torch.no_grad():
            predictions = self.model(X_tensor)
        return predictions.numpy()