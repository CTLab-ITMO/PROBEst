"""Tests for genome operations module."""

import pytest
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from unittest.mock import patch, MagicMock, mock_open
from PROBESt.genome_operations import genome_fetch, genome_blastn, genome_parse

@pytest.fixture
def mock_genome_file(tmp_path):
    """Create a mock genome file for testing."""
    genome_content = (
        ">mock_genome\n"
        "ATGCATGCATGCATGCATGCATGCATGCATGC"
        "GCATGCATGCATGCATGCATGCATGCATGCAT"
    )
    file_path = tmp_path / "mock_species.fasta"
    with open(file_path, "w") as f:
        f.write(genome_content)
    return str(file_path)

@pytest.fixture
def mock_probe():
    """Create a mock probe sequence."""
    return "ATGCATGCATGC"

@pytest.fixture
def mock_entrez_response():
    """Mock Entrez response."""
    return {
        "IdList": ["12345"],
        "Count": "1",
        "RetMax": "1",
        "RetStart": "0"
    }

@pytest.fixture
def mock_genome_record():
    """Mock genome record."""
    return SeqRecord(
        Seq("ATGCATGCATGCATGCATGCATGCATGCATGCGCATGCATGCATGCATGCATGCATGCATGCAT"),
        id="mock_genome",
        description="Mock genome for testing"
    )

@patch('Bio.Entrez.esearch')
@patch('Bio.Entrez.efetch')
@patch('Bio.Entrez.read')
def test_genome_fetch(mock_read, mock_efetch, mock_esearch, mock_entrez_response, tmp_path):
    """Test genome fetching functionality."""
    # Setup mocks
    mock_read.return_value = mock_entrez_response
    
    mock_efetch.return_value = MagicMock()
    mock_efetch.return_value.__enter__.return_value = MagicMock()
    mock_efetch.return_value.__enter__.return_value.read.return_value = ">mock_genome\nATGCATGC"
    
    # Test the function
    result = genome_fetch("mock_species")
    
    # Verify results
    assert os.path.exists(result)
    assert result.endswith("mock_species.fasta")
    
    # Cleanup
    if os.path.exists(result):
        os.unlink(result)

@patch('os.system')
@patch('builtins.open', new_callable=mock_open)
@patch('Bio.Blast.NCBIXML.parse')
def test_genome_blastn(mock_parse, mock_file, mock_system, mock_genome_file, mock_probe):
    """Test BLAST search functionality."""
    # Setup mock system call
    mock_system.return_value = 0
    
    # Create mock BLAST output
    mock_blast_record = MagicMock()
    mock_blast_record.alignments = [MagicMock()]
    mock_blast_record.alignments[0].hsps = [MagicMock()]
    mock_blast_record.alignments[0].hsps[0].sbjct_start = 1
    mock_blast_record.alignments[0].hsps[0].sbjct_end = 12
    mock_parse.return_value = [mock_blast_record]
    
    # Create mock genome file
    mock_file.return_value.__enter__.return_value.read.return_value = ">mock_genome\nATGCATGC"
    
    result = genome_blastn(mock_genome_file, mock_probe)
    
    assert isinstance(result, str)
    assert len(result) >= len(mock_probe)

@patch('PROBESt.genome_operations.genome_fetch')
@patch('PROBESt.genome_operations.genome_blastn')
def test_genome_parse_single_species(mock_blastn, mock_fetch, mock_probe, tmp_path):
    """Test genome parsing with a single species."""
    # Setup mocks
    mock_fetch.return_value = str(tmp_path / "mock_file.fasta")
    mock_blastn.return_value = "ATGCATGCATGCATGCATGC"
    
    # Test the function
    species = ["Escherichia coli"]
    result = genome_parse(species, mock_probe, extend=5)
    
    # Verify results
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ['species', 'probe_id', 'species_dna_string', 
                                  'extend', 'probe_dna_string']
    assert len(result) == 1
    assert result.iloc[0]['species'] == species[0]
    assert result.iloc[0]['probe_dna_string'] == mock_probe
    assert result.iloc[0]['extend'] == 5

@patch('PROBESt.genome_operations.genome_fetch')
@patch('PROBESt.genome_operations.genome_blastn')
def test_genome_parse_output_file(mock_blastn, mock_fetch, mock_probe, tmp_path):
    """Test if genome parse creates output CSV file."""
    # Setup mocks
    mock_fetch.return_value = str(tmp_path / "mock_file.fasta")
    mock_blastn.return_value = "ATGCATGCATGCATGCATGC"
    
    # Test the function
    species = ["Escherichia coli"]
    probe = "ATGCATGCATGC"
    output_file = tmp_path / "genome_parse_results.csv"
    
    result = genome_parse(species, probe)
    
    # Verify results
    assert os.path.exists(output_file)
    df = pd.read_csv(output_file)
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ['species', 'probe_id', 'species_dna_string', 
                              'extend', 'probe_dna_string']

@patch('PROBESt.genome_operations.genome_fetch')
def test_invalid_species(mock_fetch):
    """Test handling of invalid species names."""
    # Setup mock to raise exception
    mock_fetch.side_effect = ValueError("No genome found")
    
    species = ["ThisIsNotARealSpecies12345"]
    probe = "ATGCATGCATGC"
    result = genome_parse(species, probe)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 0  # Should return empty DataFrame for invalid species 