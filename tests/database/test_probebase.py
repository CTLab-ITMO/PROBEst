"""Tests for ProbeBase database parsing."""

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
from bs4 import BeautifulSoup
from scripts.databases.probeBase import parse_probebase_page

@pytest.fixture
def mock_response():
    """Create a mock response object."""
    mock = MagicMock()
    mock.status_code = 200
    return mock

@patch('requests.get')
def test_response_problem(mock_get, mock_response):
    """Test page with response problem."""
    # Setup mock response
    mock_response.content = b"<html><body>Error page</body></html>"
    mock_get.return_value = mock_response
    
    data = "https://probebase.csb.univie.ac.at/pb_report/probe"
    resp = parse_probebase_page(data)
    assert isinstance(resp, Warning)

@patch('requests.get')
def test_response_empty(mock_get, mock_response):
    """Test page with empty table."""
    # Setup mock response with empty table
    mock_response.content = b"""
    <html>
        <body>
            <table>
                <tr><th>Header</th></tr>
                <tr><td>No data</td></tr>
            </table>
        </body>
    </html>
    """
    mock_get.return_value = mock_response
    
    data = "https://probebase.csb.univie.ac.at/pb_report/probe/1"
    resp = parse_probebase_page(data)
    assert isinstance(resp, Warning)

@patch('requests.get')
def test_response_content(mock_get, mock_response):
    """Test page with content."""
    # Setup mock response with valid table data
    mock_response.content = b"""
    <html>
        <body>
            <table>
                <tr><th>Probe ID</th><th>Sequence</th></tr>
                <tr><td>1</td><td>ATGC</td></tr>
            </table>
        </body>
    </html>
    """
    mock_get.return_value = mock_response
    
    data = "https://probebase.csb.univie.ac.at/pb_report/probe/2"
    resp = parse_probebase_page(data)
    assert isinstance(resp, pd.Series)