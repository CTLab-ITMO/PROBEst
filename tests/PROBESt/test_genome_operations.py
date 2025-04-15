"""Tests for genome operations module."""

import pytest
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from unittest.mock import patch, MagicMock, mock_open
from io import StringIO, BytesIO
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
def mock_blast_xml():
    """Mock BLAST XML output."""
    return b"""<?xml version="1.0"?>
<BlastOutput>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>100</Hsp_bit-score>
              <Hsp_score>50</Hsp_score>
              <Hsp_evalue>1e-10</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>12</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>12</Hsp_hit-to>
              <Hsp_identity>12</Hsp_identity>
              <Hsp_align-len>12</Hsp_align-len>
              <Hsp_qseq>ATGCATGCATGC</Hsp_qseq>
              <Hsp_hseq>ATGCATGCATGC</Hsp_hseq>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""

@pytest.fixture
def mock_entrez_xml():
    """Mock Entrez XML response."""
    return b"""<?xml version="1.0"?>
<eSearchResult>
    <Count>1</Count>
    <RetMax>1</RetMax>
    <RetStart>0</RetStart>
    <IdList>
        <Id>12345</Id>
    </IdList>
</eSearchResult>
"""

@patch('Bio.Entrez.esearch')
@patch('Bio.Entrez.efetch')
@patch('Bio.Entrez.read')
@patch('builtins.open', new_callable=mock_open)
def test_genome_fetch(mock_file, mock_read, mock_efetch, mock_esearch, mock_entrez_xml):
    """Test genome fetching functionality."""
    # Setup mocks
    mock_esearch.return_value = BytesIO(mock_entrez_xml)
    mock_read.return_value = {"IdList": ["12345"]}
    mock_efetch.return_value = StringIO(">mock_genome\nATGCATGC")
    
    # Test the function
    result = genome_fetch("mock_species")
    
    # Verify results
    assert result == "mock_species.fasta"
    
    # Verify file was written correctly
    mock_file.assert_called_once_with("mock_species.fasta", "w")
    mock_file.return_value.__enter__.return_value.write.assert_called_once_with(">mock_genome\nATGCATGC")

@patch('PROBESt.genome_operations.NcbiblastnCommandline')
@patch('builtins.open', new_callable=mock_open)
@patch('Bio.SeqIO.read')
def test_genome_blastn(mock_seqio, mock_file, mock_blastn, mock_genome_file, mock_blast_xml):
    """Test BLAST search functionality."""
    # Setup mock BLAST command
    mock_blastn.return_value = MagicMock()
    mock_blastn.return_value.return_value = ("", "")  # stdout, stderr
    
    # Setup mock file for XML reading
    mock_file.return_value.__enter__.return_value.read.return_value = mock_blast_xml
    
    # Setup mock SeqIO.read
    mock_seqio.return_value = SeqRecord(
        Seq("ATGCATGCATGCATGCATGC"),
        id="mock_genome"
    )
    
    # Test the function
    result = genome_blastn(mock_genome_file, "ATGCATGCATGC")
    
    # Verify results
    assert isinstance(result, str)
    assert len(result) > 0
    assert "ATGC" in result
    
    # Verify BLAST command was called correctly
    mock_blastn.assert_called_once()
    args, kwargs = mock_blastn.call_args
    assert kwargs['outfmt'] == '5'  # XML format

@patch('PROBESt.genome_operations.genome_fetch')
@patch('PROBESt.genome_operations.genome_blastn')
def test_genome_parse_single_species(mock_blastn, mock_fetch):
    """Test parsing BLAST results into DataFrame."""
    # Setup mocks
    mock_fetch.return_value = "mock_species.fasta"
    mock_blastn.return_value = "ATGCATGCATGC"
    
    # Test the function
    species = ["mock_species"]
    probe = "ATGCATGCATGC"
    df = genome_parse(species, probe)
    
    # Verify DataFrame structure and content
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert list(df.columns) == [
        'species', 'probe_id', 'species_dna_string', 'extend', 'probe_dna_string'
    ]
    assert df['species'].iloc[0] == "mock_species"
    assert df['probe_dna_string'].iloc[0] == probe

@patch('PROBESt.genome_operations.genome_fetch')
@patch('PROBESt.genome_operations.genome_blastn')
def test_genome_parse_output_file(mock_blastn, mock_fetch):
    """Test creating output file from BLAST results."""
    # Setup mocks
    mock_fetch.return_value = "mock_species.fasta"
    mock_blastn.return_value = "ATGCATGCATGC"
    
    # Test the function
    species = ["mock_species"]
    probe = "ATGCATGCATGC"
    df = genome_parse(species, probe)
    
    # Verify DataFrame structure and content
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1
    assert list(df.columns) == [
        'species', 'probe_id', 'species_dna_string', 'extend', 'probe_dna_string'
    ]
    
    # Verify CSV file was created
    assert os.path.exists("genome_parse_results.csv")
    
    # Cleanup
    if os.path.exists("genome_parse_results.csv"):
        os.unlink("genome_parse_results.csv")

@patch('Bio.Entrez.esearch')
def test_invalid_species(mock_esearch):
    """Test handling of invalid species."""
    # Setup mock
    mock_esearch.return_value = BytesIO(b"""<?xml version="1.0"?>
<eSearchResult>
    <Count>0</Count>
    <RetMax>0</RetMax>
    <RetStart>0</RetStart>
    <IdList>
    </IdList>
</eSearchResult>
""")
    
    # Test the function
    with pytest.raises(ValueError, match="No genome found for species: invalid_species"):
        genome_fetch("invalid_species") 