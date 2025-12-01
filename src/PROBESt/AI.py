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
            nn.LayerNorm(128),  # Use LayerNorm instead of BatchNorm for batch size 1 compatibility
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            
            nn.Linear(128, 64),
            nn.LayerNorm(64),  # Use LayerNorm instead of BatchNorm for batch size 1 compatibility
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            
            nn.Linear(64, 32),
            nn.LayerNorm(32),  # Use LayerNorm instead of BatchNorm for batch size 1 compatibility
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            
            nn.Linear(32, 16),
            nn.LayerNorm(16),  # Use LayerNorm instead of BatchNorm for batch size 1 compatibility
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
        
        # Track learning curves
        self.train_losses = []
        self.val_losses = []
        self.val_metrics_history = {'f1': [], 'accuracy': [], 'recall': [], 'precision': []}
        
    def train(self, X, y, epochs=100, batch_size=32, val_data=None, track_curves=False):
        X_scaled = self.preprocess_data(X)
        X_tensor = torch.FloatTensor(X_scaled)
        y_tensor = torch.FloatTensor(y.values).reshape(-1, 1)
        
        dataset = torch.utils.data.TensorDataset(X_tensor, y_tensor)
        # Drop last batch if it's smaller than batch_size to avoid BatchNorm issues
        loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=True)
        
        # Prepare validation data if provided
        val_X_tensor = None
        val_y_tensor = None
        if val_data is not None and track_curves:
            val_X = val_data.drop(columns=['type'])
            val_y = val_data['type']
            val_X_scaled = self.scaler.transform(val_X)
            val_X_tensor = torch.FloatTensor(val_X_scaled)
            val_y_tensor = torch.FloatTensor(val_y.values).reshape(-1, 1)
        
        # Reset tracking if tracking curves
        if track_curves:
            self.train_losses = []
            self.val_losses = []
            self.val_metrics_history = {'f1': [], 'accuracy': [], 'recall': [], 'precision': []}
        
        for e in range(epochs):
            self.model.train()
            total_loss = 0
            for bx, by in loader:
                self.optimizer.zero_grad()
                logits = self.model(bx)
                loss = self.criterion(logits, by)
                loss.backward()
                self.optimizer.step()
                total_loss += loss.item()
            
            avg_loss = total_loss / len(loader)
            if track_curves:
                self.train_losses.append(avg_loss)
            
            # Evaluate on validation set if provided
            if val_X_tensor is not None and track_curves:
                self.model.eval()
                with torch.no_grad():
                    val_logits = self.model(val_X_tensor)
                    val_loss = self.criterion(val_logits, val_y_tensor).item()
                    self.val_losses.append(val_loss)
                    
                    # Calculate metrics
                    val_pred_proba = torch.sigmoid(val_logits).numpy()
                    val_pred = (val_pred_proba > 0.5).astype(int)
                    val_y_np = val_y_tensor.numpy().flatten()
                    
                    from sklearn.metrics import f1_score, accuracy_score, recall_score, precision_score
                    self.val_metrics_history['f1'].append(f1_score(val_y_np, val_pred))
                    self.val_metrics_history['accuracy'].append(accuracy_score(val_y_np, val_pred))
                    self.val_metrics_history['recall'].append(recall_score(val_y_np, val_pred))
                    self.val_metrics_history['precision'].append(precision_score(val_y_np, val_pred, zero_division=0))
            
            if e % 20 == 0:
                print(f"Epoch {e}: loss = {avg_loss:.4f}")

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
            shuffle=True,
            drop_last=True  # Drop last incomplete batch to avoid BatchNorm issues
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