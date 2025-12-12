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


"""Tests for ProbeBase database parsing."""

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
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