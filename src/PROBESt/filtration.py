import os
import pandas as pd
from typing import Tuple, Any, Dict
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    accuracy_score, f1_score, recall_score, r2_score,
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
        'r2': r2_score(y_test, y_pred)
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
        'r2': r2_score(y, y_pred)
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