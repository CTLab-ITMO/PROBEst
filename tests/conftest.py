"""Common pytest configurations and fixtures."""

import pytest
import os
import sys
from pathlib import Path

# Add the project root directory to the Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

@pytest.fixture(autouse=True)
def setup_test_env():
    """Setup test environment."""
    # Set test-specific environment variables
    os.environ['ENTREZ_EMAIL'] = 'test@example.com'
    
    # Create temporary directory for test files
    test_dir = project_root / 'tests' / 'test_data'
    test_dir.mkdir(exist_ok=True)
    
    yield
    
    # Cleanup after tests
    for file in test_dir.glob('*'):
        if file.is_file():
            file.unlink()
    test_dir.rmdir()

@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return project_root / 'tests' / 'test_data'
