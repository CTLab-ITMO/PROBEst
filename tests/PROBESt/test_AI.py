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


def test_calculate_gc_content():
    from PROBESt.misc import calculate_gc_content
    
    # Test with GC-rich sequence
    assert calculate_gc_content("GCGCGC") == 100.0
    # Test with AT-rich sequence
    assert calculate_gc_content("ATATAT") == 0.0
    # Test with mixed sequence
    assert abs(calculate_gc_content("ATGC") - 50.0) < 0.1
    # Test with empty sequence
    assert calculate_gc_content("") == 0.0


def test_extend_blast_output_with_parameters(tmp_path):
    from PROBESt.misc import extend_blast_output_with_parameters
    import pandas as pd
    
    # Create a test FASTA file
    fasta_file = tmp_path / "test.fa"
    with open(fasta_file, 'w') as f:
        f.write(">probe1\nATGCATGC\n")
        f.write(">probe2\nGCGCGCGC\n")
    
    # Create a test BLAST DataFrame
    blast_df = pd.DataFrame({
        0: ['probe1', 'probe2'],
        1: ['seq1', 'seq2'],
        2: [0.001, 0.01],
        3: [1, 10],
        4: [8, 17],
        5: [100, 90],
        6: [0, 1]
    })
    
    # Extend with parameters
    extended_df = extend_blast_output_with_parameters(blast_df, str(fasta_file))
    
    # Check that new columns were added
    assert 'sseq' in extended_df.columns
    assert 'GCcontent' in extended_df.columns
    assert 'Lengthnt' in extended_df.columns
    assert 'hairpin_prob' in extended_df.columns
    
    # Check that sequences were loaded
    assert extended_df.loc[0, 'sseq'] == 'ATGCATGC'
    assert extended_df.loc[1, 'sseq'] == 'GCGCGCGC'
    
    # Check GC content calculation
    assert abs(extended_df.loc[0, 'GCcontent'] - 50.0) < 0.1
    assert extended_df.loc[1, 'GCcontent'] == 100.0 