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


import os
import pandas as pd
import numpy as np
from typing import Tuple, Any, Dict
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    accuracy_score, f1_score, recall_score, precision_score,
    roc_curve, auc, confusion_matrix
)
from sklearn.model_selection import train_test_split

def train_filtration_AI(
    model: Any,
    data: pd.DataFrame,
    target_col: str = 'type',
    test_size: float = 0.2,
    random_state: int = 42
) -> Tuple[Any, Dict[str, float]]:
    """
    Train and evaluate an AI model for probe filtration.
    
    Args:
        model: AI model instance with train() and predict() methods
        data: Input DataFrame
        target_col: Name of the target column
        test_size: Proportion of data to use for testing
        random_state: Random seed for reproducibility
        
    Returns:
        Tuple of (trained_model, metrics_dict)
    """
    # Prepare data
    X = data.drop(columns=[target_col])
    y = data[target_col]
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state
    )
    
    # Train model
    model.train(X_train, y_train)
    
    # Evaluate
    y_pred = model.predict(X_test)
    y_pred_proba = model.predict_proba(X_test)
    
    metrics = {
        'accuracy': accuracy_score(y_test, y_pred),
        'f1': f1_score(y_test, y_pred),
        'recall': recall_score(y_test, y_pred),
        'precision': precision_score(y_test, y_pred, zero_division=0)
    }
    
    return model, metrics

def validate_filtration_AI(
    model: Any,
    data: pd.DataFrame,
    target_col: str = 'type',
    output_dir: str = 'tests_outs',
    output_name: str = 'validation_plots.png'
) -> Dict[str, float]:
    """
    Validate the trained AI model and generate performance plots.
    
    Args:
        model: Trained AI model
        data: Validation dataset
        target_col: Name of the target column
        output_dir: Directory to save plots
        
    Returns:
        Dictionary of validation metrics
    """
    os.makedirs(output_dir, exist_ok=True)
    
    X = data.drop(columns=[target_col])
    y = data[target_col]
    
    # Get predictions
    y_pred = model.predict(X)
    y_pred_proba = model.predict_proba(X)
    
    # Handle probability predictions for ROC curve
    if y_pred_proba.ndim > 1 and y_pred_proba.shape[1] > 1:
        y_pred_proba = y_pred_proba[:, 1]  # Use probability of positive class for binary classification
    else:
        y_pred_proba = y_pred_proba.flatten()  # Flatten single probability array
    
    # Calculate metrics
    metrics = {
        'accuracy': accuracy_score(y, y_pred),
        'f1': f1_score(y, y_pred),
        'recall': recall_score(y, y_pred),
        'precision': precision_score(y, y_pred, zero_division=0)
    }
    
    # Generate plots
    plt.figure(figsize=(15, 10))
    
    # ROC curve
    plt.subplot(2, 2, 1)
    fpr, tpr, _ = roc_curve(y, y_pred_proba)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, label=f'ROC curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend()
    
    # Confusion matrix
    plt.subplot(2, 2, 2)
    cm = confusion_matrix(y, y_pred)
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.title('Confusion Matrix')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    
    # Metrics bar plot
    plt.subplot(2, 2, 3)
    metrics_values = list(metrics.values())
    metrics_names = list(metrics.keys())
    plt.bar(metrics_names, metrics_values)
    plt.title('Model Metrics')
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, output_name))
    plt.close()
    
    return metrics

def apply_filtration_AI(
    model: Any,
    data: pd.DataFrame,
    target_col: str = 'type',
    threshold: float = 0.5
) -> pd.DataFrame:
    """
    Apply the trained AI model to filter probes.
    
    Args:
        model: Trained AI model
        data: Input DataFrame
        target_col: Name of the target column
        threshold: Probability threshold for classification
        
    Returns:
        DataFrame with predictions
    """
    # Drop target column if it exists
    X = data.drop(columns=[target_col]) if target_col in data.columns else data
    
    predictions = model.predict(X)
    probabilities = model.predict_proba(X)
    
    # Handle different probability output formats
    if probabilities.ndim > 1 and probabilities.shape[1] > 1:
        # For models that return probabilities for both classes (e.g., LogisticRegression)
        proba = probabilities[:, 1]  # Use probability of positive class
    else:
        # For models that return single probability (e.g., DNN)
        proba = probabilities.flatten()
    
    result = data.copy()
    result['predicted_type'] = predictions
    result['prediction_probability'] = proba
    
    return result

def plot_combined_roc_curves(
    models_dict: Dict[str, Tuple[Any, pd.DataFrame]],
    output_dir: str = 'tests_outs',
    output_name: str = 'combined_roc_curves.png'
) -> None:
    """
    Plot ROC curves for all models on a single plot.
    
    Args:
        models_dict: Dictionary mapping model names to (trained_model, validation_data) tuples
        output_dir: Directory to save plots
        output_name: Name of output file
    """
    os.makedirs(output_dir, exist_ok=True)
    
    plt.figure(figsize=(10, 8))
    
    for model_name, (model, val_data) in models_dict.items():
        X = val_data.drop(columns=['type'])
        y = val_data['type']
        
        y_pred_proba = model.predict_proba(X)
        
        # Handle probability predictions for ROC curve
        if y_pred_proba.ndim > 1 and y_pred_proba.shape[1] > 1:
            y_pred_proba = y_pred_proba[:, 1]
        else:
            y_pred_proba = y_pred_proba.flatten()
        
        fpr, tpr, _ = roc_curve(y, y_pred_proba)
        roc_auc = auc(fpr, tpr)
        
        plt.plot(fpr, tpr, label=f'{model_name} (AUC = {roc_auc:.3f})', linewidth=2)
    
    plt.plot([0, 1], [0, 1], 'k--', label='Random Classifier', linewidth=1)
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title('ROC Curves Comparison', fontsize=14, fontweight='bold')
    plt.legend(loc='lower right', fontsize=10)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, output_name), dpi=300, bbox_inches='tight')
    plt.close()

def plot_combined_metrics(
    results_dict: Dict[str, Dict[str, float]],
    output_dir: str = 'tests_outs',
    output_name: str = 'combined_metrics.png'
) -> None:
    """
    Plot combined metrics (f1, recall, accuracy) for all models.
    
    Args:
        results_dict: Dictionary mapping model names to their metrics dictionaries
        output_dir: Directory to save plots
        output_name: Name of output file
    """
    os.makedirs(output_dir, exist_ok=True)
    
    model_names = list(results_dict.keys())
    metrics_to_plot = ['f1', 'recall', 'accuracy']
    
    # Prepare data for plotting
    metric_values = {metric: [results_dict[name][metric] for name in model_names] 
                     for metric in metrics_to_plot}
    
    # Create grouped bar chart
    x = np.arange(len(model_names))
    width = 0.25
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    for i, metric in enumerate(metrics_to_plot):
        offset = (i - 1) * width
        bars = ax.bar(x + offset, metric_values[metric], width, label=metric.upper(), alpha=0.8)
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{height:.3f}',
                   ha='center', va='bottom', fontsize=8)
    
    ax.set_xlabel('Model', fontsize=12)
    ax.set_ylabel('Score', fontsize=12)
    ax.set_title('Model Performance Metrics Comparison', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(model_names, rotation=45, ha='right')
    ax.set_ylim([0, 1.1])
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, output_name), dpi=300, bbox_inches='tight')
    plt.close()

def plot_learning_curves(
    model: Any,
    output_dir: str = 'tests_outs',
    output_name: str = 'learning_curves.png'
) -> None:
    """
    Plot learning curves (loss and metrics over epochs) for a trained model.
    
    Args:
        model: Trained model with train_losses, val_losses, and val_metrics_history attributes
        output_dir: Directory to save plots
        output_name: Name of output file
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if not hasattr(model, 'train_losses') or len(model.train_losses) == 0:
        print("Warning: Model does not have learning curve data. Skipping plot.")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    epochs = range(1, len(model.train_losses) + 1)
    
    # Plot training and validation loss
    ax1 = axes[0, 0]
    ax1.plot(epochs, model.train_losses, 'b-', label='Training Loss', linewidth=2)
    if hasattr(model, 'val_losses') and len(model.val_losses) > 0:
        ax1.plot(epochs, model.val_losses, 'r-', label='Validation Loss', linewidth=2)
    ax1.set_xlabel('Epoch', fontsize=11)
    ax1.set_ylabel('Loss', fontsize=11)
    ax1.set_title('Loss Curves', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # Plot F1 score
    ax2 = axes[0, 1]
    if hasattr(model, 'val_metrics_history') and len(model.val_metrics_history.get('f1', [])) > 0:
        ax2.plot(epochs, model.val_metrics_history['f1'], 'g-', label='F1 Score', linewidth=2)
        ax2.set_xlabel('Epoch', fontsize=11)
        ax2.set_ylabel('F1 Score', fontsize=11)
        ax2.set_title('F1 Score Over Epochs', fontsize=12, fontweight='bold')
        ax2.legend()
        ax2.grid(alpha=0.3)
        ax2.set_ylim([0, 1])
    
    # Plot Accuracy
    ax3 = axes[1, 0]
    if hasattr(model, 'val_metrics_history') and len(model.val_metrics_history.get('accuracy', [])) > 0:
        ax3.plot(epochs, model.val_metrics_history['accuracy'], 'm-', label='Accuracy', linewidth=2)
        ax3.set_xlabel('Epoch', fontsize=11)
        ax3.set_ylabel('Accuracy', fontsize=11)
        ax3.set_title('Accuracy Over Epochs', fontsize=12, fontweight='bold')
        ax3.legend()
        ax3.grid(alpha=0.3)
        ax3.set_ylim([0, 1])
    
    # Plot Recall and Precision
    ax4 = axes[1, 1]
    if hasattr(model, 'val_metrics_history') and len(model.val_metrics_history.get('recall', [])) > 0:
        ax4.plot(epochs, model.val_metrics_history['recall'], 'c-', label='Recall', linewidth=2)
        if len(model.val_metrics_history.get('precision', [])) > 0:
            ax4.plot(epochs, model.val_metrics_history['precision'], 'orange', label='Precision', linewidth=2)
        ax4.set_xlabel('Epoch', fontsize=11)
        ax4.set_ylabel('Score', fontsize=11)
        ax4.set_title('Recall and Precision Over Epochs', fontsize=12, fontweight='bold')
        ax4.legend()
        ax4.grid(alpha=0.3)
        ax4.set_ylim([0, 1])
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, output_name), dpi=300, bbox_inches='tight')
    plt.close()

def plot_gail_architecture_comparison(
    results_dict: Dict[str, Dict[str, float]],
    output_dir: str = 'tests_outs',
    output_name: str = 'gail_architecture_comparison.png'
) -> None:
    """
    Plot comparison of all GAIL architecture variants.
    
    Args:
        results_dict: Dictionary mapping GAIL variant names to their metrics dictionaries
        output_dir: Directory to save plots
        output_name: Name of output file
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Filter to only GAIL variants
    gail_results = {k: v for k, v in results_dict.items() 
                   if k.startswith("GAIL") or "GAIL" in k}
    
    if len(gail_results) == 0:
        print("Warning: No GAIL architectures found for comparison.")
        return
    
    model_names = list(gail_results.keys())
    metrics_to_plot = ['f1', 'recall', 'accuracy', 'precision']
    
    # Prepare data for plotting
    metric_values = {metric: [gail_results[name].get(metric, 0) for name in model_names] 
                     for metric in metrics_to_plot}
    
    # Create grouped bar chart
    x = np.arange(len(model_names))
    width = 0.2
    
    fig, ax = plt.subplots(figsize=(16, 8))
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    for i, metric in enumerate(metrics_to_plot):
        offset = (i - 1.5) * width
        bars = ax.bar(x + offset, metric_values[metric], width, 
                     label=metric.upper(), alpha=0.8, color=colors[i])
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{height:.3f}',
                       ha='center', va='bottom', fontsize=7)
    
    ax.set_xlabel('GAIL Architecture Variant', fontsize=12)
    ax.set_ylabel('Score', fontsize=12)
    ax.set_title('GAIL Architecture Comparison', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(model_names, rotation=45, ha='right', fontsize=9)
    ax.set_ylim([0, 1.1])
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, output_name), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Also create a heatmap-style comparison
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create matrix for heatmap
    heatmap_data = np.array([[gail_results[name].get(metric, 0) 
                             for metric in metrics_to_plot] 
                            for name in model_names])
    
    im = ax.imshow(heatmap_data, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)
    
    # Set ticks and labels
    ax.set_xticks(np.arange(len(metrics_to_plot)))
    ax.set_yticks(np.arange(len(model_names)))
    ax.set_xticklabels([m.upper() for m in metrics_to_plot])
    ax.set_yticklabels(model_names)
    
    # Rotate x-axis labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    
    # Add text annotations
    for i in range(len(model_names)):
        for j in range(len(metrics_to_plot)):
            text = ax.text(j, i, f'{heatmap_data[i, j]:.3f}',
                          ha="center", va="center", color="black", fontsize=8)
    
    ax.set_title("GAIL Architecture Performance Heatmap", fontsize=14, fontweight='bold', pad=20)
    plt.colorbar(im, ax=ax, label='Score')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, output_name.replace('.png', '_heatmap.png')), 
               dpi=300, bbox_inches='tight')
    plt.close() 