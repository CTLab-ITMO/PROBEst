import pytest
import pandas as pd
import numpy as np
from PROBESt.AI import LogisticRegressionModel, PerceptronModel, DeepNeuralNetworkModel
import torch

@pytest.fixture
def sample_data():
    # Create sample data for testing
    np.random.seed(42)
    n_samples = 100
    n_features = 5
    
    X = pd.DataFrame(
        np.random.randn(n_samples, n_features),
        columns=[f'feature_{i}' for i in range(n_features)]
    )
    y = pd.Series(np.random.randint(0, 2, n_samples))
    
    return X, y

def test_logistic_regression(sample_data):
    X, y = sample_data
    model = LogisticRegressionModel()
    
    # Test training
    model.train(X, y)
    
    # Test prediction
    predictions = model.predict(X)
    assert len(predictions) == len(y)
    assert all(pred in [0, 1] for pred in predictions)
    
    # Test probability prediction
    proba = model.predict_proba(X)
    assert proba.shape == (len(y), 2)
    assert all(0 <= p <= 1 for p in proba.flatten())

def test_perceptron(sample_data):
    X, y = sample_data
    model = PerceptronModel(input_size=X.shape[1])
    
    # Test training
    model.train(X, y, epochs=10)
    
    # Test prediction
    predictions = model.predict(X)
    assert len(predictions) == len(y)
    assert all(pred in [0, 1] for pred in predictions)
    
    # Test probability prediction
    proba = model.predict_proba(X)
    assert proba.shape == (len(y), 1)
    assert all(0 <= p <= 1 for p in proba.flatten())

def test_deep_neural_network(sample_data):
    X, y = sample_data
    model = DeepNeuralNetworkModel(input_size=X.shape[1])
    
    # Test training
    model.train(X, y, epochs=10, batch_size=16)
    
    # Test prediction
    predictions = model.predict(X)
    assert len(predictions) == len(y)
    assert all(pred in [0, 1] for pred in predictions)
    
    # Test probability prediction
    proba = model.predict_proba(X)
    assert proba.shape == (len(y), 1)
    assert all(0 <= p <= 1 for p in proba.flatten())
    
    # Test model evaluation mode
    model.model.eval()
    with torch.no_grad():
        eval_predictions = model.predict(X)
    assert len(eval_predictions) == len(y)
    assert all(pred in [0, 1] for pred in eval_predictions) 