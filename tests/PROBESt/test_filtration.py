import pytest
import pandas as pd
import numpy as np
import os
from PROBESt.AI import LogisticRegressionModel
from PROBESt.filtration import train_filtration_AI, validate_filtration_AI, apply_filtration_AI

@pytest.fixture
def sample_data():
    # Create sample data for testing
    np.random.seed(42)
    n_samples = 100
    n_features = 5
    
    data = pd.DataFrame(
        np.random.randn(n_samples, n_features),
        columns=[f'feature_{i}' for i in range(n_features)]
    )
    data['type'] = np.random.randint(0, 2, n_samples)
    
    return data

def test_train_filtration_AI(sample_data):
    model = LogisticRegressionModel()
    trained_model, metrics = train_filtration_AI(model, sample_data)
    
    assert all(metric in metrics for metric in ['accuracy', 'f1', 'recall', 'precision'])
    assert all(0 <= value <= 1 for value in [metrics['accuracy'], metrics['f1'], metrics['recall'], metrics['precision']])

def test_validate_filtration_AI(sample_data):
    model = LogisticRegressionModel()
    trained_model, _ = train_filtration_AI(model, sample_data)
    
    metrics = validate_filtration_AI(trained_model, sample_data)
    
    assert all(metric in metrics for metric in ['accuracy', 'f1', 'recall', 'precision'])
    assert all(0 <= value <= 1 for value in [metrics['accuracy'], metrics['f1'], metrics['recall'], metrics['precision']])
    assert os.path.exists('tests_outs/validation_plots.png')

def test_apply_filtration_AI(sample_data):
    model = LogisticRegressionModel()
    trained_model, _ = train_filtration_AI(model, sample_data)
    
    result = apply_filtration_AI(trained_model, sample_data)
    
    assert 'predicted_type' in result.columns
    assert 'prediction_probability' in result.columns
    assert len(result) == len(sample_data)
    assert all(pred in [0, 1] for pred in result['predicted_type'])
    assert all(0 <= prob <= 1 for prob in result['prediction_probability']) 