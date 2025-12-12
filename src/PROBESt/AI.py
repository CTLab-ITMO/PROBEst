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


import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
import torch
import torch.nn as nn
from typing import Tuple, Any, Dict, Optional

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


def detect_architecture_from_state_dict(state_dict: Dict, input_size: int) -> str:
    """
    Detect model architecture from state_dict shapes.
    
    Args:
        state_dict: Model state dictionary
        input_size: Input size of the model
    
    Returns:
        str: Detected architecture name
    """
    # Get the first layer weight shape to determine architecture
    if 'model_state_dict' in state_dict:
        first_layer_key = 'model_state_dict.net.0.weight'
        if first_layer_key not in state_dict['model_state_dict']:
            first_layer_key = 'model_state_dict.net.0.weight'
        weights = state_dict['model_state_dict']
    else:
        weights = state_dict
    
    # Find first linear layer
    first_layer_shape = None
    for key in weights.keys():
        if 'net.0.weight' in key or '0.weight' in key:
            first_layer_shape = weights[key].shape
            break
    
    if first_layer_shape is None:
        raise ValueError("Could not detect architecture from state_dict")
    
    hidden1 = first_layer_shape[0]
    
    # Find second layer to determine hidden2
    hidden2 = None
    for key in weights.keys():
        if 'net.2.weight' in key or '2.weight' in key:
            hidden2 = weights[key].shape[0]
            break
    
    # Find third layer to determine if it's a 3-layer architecture
    has_third_layer = False
    for key in weights.keys():
        if 'net.4.weight' in key or '4.weight' in key:
            if 'net.6.weight' in str(weights.keys()) or '6.weight' in str(weights.keys()):
                has_third_layer = True
            break
    
    # Map to architecture based on layer sizes
    # GAILDiscriminator: [256, 128]
    # GAILWide: [512, 256]
    # GAILWideExtra: [768, 384]
    # GAILNarrow: [128, 64]
    # GAILDeep: [256, 128, 64]
    # GAILWideDeep: [512, 256, 128]
    # GAILWideBalanced: [512, 256, 128]
    
    if hidden1 == 768 and hidden2 == 384:
        return "GAILWideExtra"
    elif hidden1 == 512 and hidden2 == 256:
        if has_third_layer:
            # Check if it's WideDeep or WideBalanced by checking for dropout/norm
            has_dropout = any('dropout' in str(k).lower() for k in weights.keys())
            has_norm = any('norm' in str(k).lower() for k in weights.keys())
            if has_dropout:
                return "GAILWideDropout"
            elif has_norm:
                return "GAILWideBatchNorm"
            else:
                return "GAILWideDeep"  # Default to WideDeep for 3-layer
        else:
            return "GAILWide"
    elif hidden1 == 256 and hidden2 == 128:
        if has_third_layer:
            return "GAILDeep"
        else:
            return "GAILDiscriminator"
    elif hidden1 == 128 and hidden2 == 64:
        return "GAILNarrow"
    else:
        # Try to match based on common patterns
        if hidden1 == 256 and hidden2 == 128:
            return "GAILDiscriminator"
        elif hidden1 >= 512:
            return "GAILWide"
        else:
            return "GAILDiscriminator"  # Default fallback


def load_torch_classifier(model_path: str, input_size: int, model_architecture: Optional[str] = None) -> TorchClassifier:
    """
    Load a saved TorchClassifier model from file.
    
    Args:
        model_path (str): Path to saved model file (.pt)
        input_size (int): Input size (number of features)
        model_architecture (Optional[str]): Model architecture name. If None, will be auto-detected.
            Options: "GAILDiscriminator", "GAILWide", "GAILDeep", "GAILNarrow", 
            "GAILWithDropout", "GAILWideDeep", "GAILWideDropout", 
            "GAILWideBatchNorm", "GAILWideExtra", "GAILWideBalanced"
    
    Returns:
        TorchClassifier: Loaded and initialized model
    """
    from PROBESt.models_registry import (
        GAILDiscriminator, GAILWide, GAILDeep, GAILNarrow, GAILWithDropout,
        GAILWideDeep, GAILWideDropout, GAILWideBatchNorm, GAILWideExtra, GAILWideBalanced
    )
    
    # Map architecture names to classes
    architecture_map = {
        "GAILDiscriminator": GAILDiscriminator,
        "GAILWide": GAILWide,
        "GAILDeep": GAILDeep,
        "GAILNarrow": GAILNarrow,
        "GAILWithDropout": GAILWithDropout,
        "GAILWideDeep": GAILWideDeep,
        "GAILWideDropout": GAILWideDropout,
        "GAILWideBatchNorm": GAILWideBatchNorm,
        "GAILWideExtra": GAILWideExtra,
        "GAILWideBalanced": GAILWideBalanced,
    }
    
    # Load saved model data (weights_only=False to allow loading StandardScaler)
    checkpoint = torch.load(model_path, map_location='cpu', weights_only=False)
    
    # Auto-detect architecture if not provided
    if model_architecture is None:
        # First try to get from checkpoint
        if 'model_architecture' in checkpoint:
            model_architecture = checkpoint['model_architecture']
            #print(f"Using saved architecture: {model_architecture}")
        else:
            # Auto-detect from state_dict
            state_dict = checkpoint.get('model_state_dict', checkpoint)
            model_architecture = detect_architecture_from_state_dict(state_dict, input_size)
            #print(f"Auto-detected architecture: {model_architecture}")
    
    if model_architecture not in architecture_map:
        raise ValueError(f"Unknown model architecture: {model_architecture}")
    
    # Create model instance
    model_class = architecture_map[model_architecture]
    torch_model = model_class(input_size)
    
    # Load state dict
    if 'model_state_dict' in checkpoint:
        torch_model.load_state_dict(checkpoint['model_state_dict'])
    else:
        # Fallback: try loading directly
        torch_model.load_state_dict(checkpoint)
    
    # Create TorchClassifier wrapper
    weight_pos = checkpoint.get('weight_pos', 5.0)
    if weight_pos is None:
        weight_pos = 5.0
    
    classifier = TorchClassifier(torch_model, 
                                learning_rate=checkpoint.get('learning_rate', 0.001),
                                weight_pos=weight_pos)
    
    # Load scaler
    if 'scaler' in checkpoint:
        classifier.scaler = checkpoint['scaler']
    
    return classifier


def apply_ai_filtration(
    blast_df: pd.DataFrame,
    model_path: str,
    input_size: int,
    model_architecture: str = "GAILDiscriminator",
    fasta_file: Optional[str] = None
) -> pd.DataFrame:
    """
    Apply AI model to BLAST results and calculate sum scores per probe.
    
    Args:
        blast_df (pd.DataFrame): BLAST output DataFrame (should be extended with parameters)
        model_path (str): Path to saved AI model
        input_size (int): Input size for the model
        model_architecture (str): Model architecture name
        fasta_file (Optional[str]): Path to FASTA file (if blast_df needs extension)
        
    Returns:
        pd.DataFrame: DataFrame with AI scores added, including sum_score per probe
    """
    from PROBESt.misc import extend_blast_output_with_parameters
    
    # Extend BLAST output if needed
    if fasta_file and 'sseq' not in blast_df.columns:
        blast_df = extend_blast_output_with_parameters(blast_df, fasta_file)
    
    # Load model
    model = load_torch_classifier(model_path, input_size, model_architecture)
    
    # Prepare features (same columns as training data, excluding 'type')
    # Note: Training data drops non-numeric columns, so we only use numeric features
    numeric_feature_columns = [
        'Formamide', 'GCcontent', 'Lengthnt', 'evalue', 'mismatches', 'length',
        'bitscore', 'identity', 'score', 'hairpin_prob',
        'dimer_DNA', 'dimer_DNA_flank', 'dimer_probe', 'dimer_probe_DNA'
    ]
    
    # Create feature DataFrame with only numeric columns
    X = pd.DataFrame()
    
    # Add numeric features
    for col in numeric_feature_columns:
        if col in blast_df.columns:
            X[col] = pd.to_numeric(blast_df[col], errors='coerce').fillna(0.0)
        else:
            X[col] = 0.0
    
    # Ensure correct order
    X = X[numeric_feature_columns]
    
    # Get predictions (probabilities)
    predictions = model.predict_proba(X)
    
    # Add AI score to DataFrame
    blast_df['ai_score'] = predictions.flatten()
    
    # Calculate sum score per probe (group by qseqid)
    if 'qseqid' in blast_df.columns:
        probe_scores = blast_df.groupby('qseqid')['ai_score'].sum().reset_index()
        probe_scores.columns = ['qseqid', 'sum_score']
        blast_df = blast_df.merge(probe_scores, on='qseqid', how='left')
    else:
        # Fallback: use first column as probe identifier
        probe_id_col = blast_df.columns[0]
        probe_scores = blast_df.groupby(probe_id_col)['ai_score'].sum().reset_index()
        probe_scores.columns = [probe_id_col, 'sum_score']
        blast_df = blast_df.merge(probe_scores, on=probe_id_col, how='left')
    
    return blast_df


def apply_ai_filtration_to_blast_file(
    positive_hits_path: str,
    model_path: str,
    fasta_file: str,
    input_size: int = 14,
    model_architecture: Optional[str] = None,
    threshold_method: str = "median"
) -> bool:
    """
    Apply AI filtration to BLAST output file and filter probes based on sum scores.
    
    This function:
    1. Reads BLAST output from positive_hits_path
    2. Extends it with calculated parameters
    3. Applies AI model to get scores
    4. Filters probes with lower sum scores (better probes)
    5. Writes filtered results back to positive_hits_path
    
    Args:
        positive_hits_path (str): Path to BLAST positive hits file
        model_path (str): Path to saved AI model file
        fasta_file (str): Path to FASTA file containing probe sequences
        input_size (int): Input size for the model (default: 14)
        model_architecture (str): Model architecture name (default: "GAILDiscriminator")
        threshold_method (str): Method for threshold calculation - "median" or "mean" (default: "median")
    
    Returns:
        bool: True if filtration was successful, False otherwise
    """
    import os
    
    try:
        # Check if files exist
        if not os.path.exists(positive_hits_path) or os.path.getsize(positive_hits_path) == 0:
            print("Warning: positive_hits.tsv is empty or doesn't exist, skipping AI filtration")
            return False
        
        if not os.path.exists(model_path):
            print(f"AI model not found at {model_path}, skipping AI filtration")
            return False
        
        if not os.path.exists(fasta_file):
            print(f"FASTA file not found at {fasta_file}, skipping AI filtration")
            return False
        
        print("Applying AI filtration...")
        
        # Read BLAST output
        # Try tab separator first, then space
        try:
            blast_df = pd.read_table(positive_hits_path, sep='\t', header=None)
        except:
            blast_df = pd.read_table(positive_hits_path, sep=' ', header=None)
        
        if blast_df.empty:
            print("Warning: BLAST output is empty, skipping AI filtration")
            return False
        
        # Apply AI filtration
        blast_df_with_scores = apply_ai_filtration(
            blast_df=blast_df,
            model_path=model_path,
            input_size=input_size,
            model_architecture=model_architecture,
            fasta_file=fasta_file
        )
        
        # Filter probes with lower sum scores (better probes)
        if 'sum_score' not in blast_df_with_scores.columns:
            print("Warning: Could not calculate sum scores, skipping AI filtration")
            return False
        
        # Calculate threshold
        if threshold_method == "median":
            threshold = blast_df_with_scores['sum_score'].median()
        elif threshold_method == "mean":
            threshold = blast_df_with_scores['sum_score'].mean()
        else:
            print(f"Warning: Unknown threshold method '{threshold_method}', using median")
            threshold = blast_df_with_scores['sum_score'].median()
        
        # Filter to keep only rows from probes with sum_score <= threshold
        good_probes = blast_df_with_scores[
            blast_df_with_scores['sum_score'] <= threshold
        ]['qseqid'].unique()
        
        # Filter BLAST results to only include good probes
        blast_df_filtered = blast_df_with_scores[
            blast_df_with_scores['qseqid'].isin(good_probes)
        ]
        
        # Keep only original BLAST columns for probe_check (first 7 columns)
        # Select columns: qseqid, sseqid, evalue, sstart, send, ppos, mismatch
        # Get original number of columns from the input
        original_num_cols = len(blast_df.columns)
        if len(blast_df_filtered.columns) >= 7:
            blast_df_output = blast_df_filtered.iloc[:, :7].copy()
        else:
            blast_df_output = blast_df_filtered.iloc[:, :min(original_num_cols, len(blast_df_filtered.columns))].copy()
        
        # Write filtered BLAST output (space-separated to match original format)
        blast_df_output.to_csv(
            positive_hits_path, 
            sep=' ', 
            header=False, 
            index=False
        )
        
        total_probes = len(blast_df_with_scores['qseqid'].unique())
        kept_probes = len(good_probes)
        print(f"AI filtration: Kept {kept_probes} probes out of {total_probes} total")
        
        return True
        
    except Exception as e:
        print(f"Warning: AI filtration failed: {e}")
        print("Continuing with standard filtration...")
        import traceback
        traceback.print_exc()
        return False