"""Tests for modeling module."""

import pytest
import os
import subprocess
import pandas as pd
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from PROBESt.modeling import (
    prepare_blast_db,
    get_probes_from_output,
    run_blast_mapping,
    extract_extended_sequence,
    extract_best_hits,
    run_modeling,
    to_RNA,
    reverse_complement_RNA
)


def test_to_RNA():
    """Test DNA to RNA conversion."""
    assert to_RNA("ATGC") == "AUGC"
    assert to_RNA("ATCGATCG") == "AUCGAUCG"
    assert to_RNA("TTTT") == "UUUU"
    assert to_RNA("atgc") == "AUGC"  # Test lowercase conversion
    assert to_RNA("ATGCatgc") == "AUGCAUGC"  # Test mixed case
    assert to_RNA("") == ""  # Test empty string


def test_reverse_complement_RNA():
    """Test RNA reverse complement conversion."""
    assert reverse_complement_RNA("AUGC") == "GCAU"
    assert reverse_complement_RNA("AUCG") == "CGAU"
    assert reverse_complement_RNA("UUUU") == "AAAA"
    assert reverse_complement_RNA("AUGCAUGC") == "GCAUGCAU"
    # Test that reverse complement of reverse complement gives original
    seq = "AUGCAUGC"
    assert reverse_complement_RNA(reverse_complement_RNA(seq)) == seq


@pytest.fixture
def tmp_dir(tmp_path):
    """Create a temporary directory for test files."""
    return str(tmp_path)


@pytest.fixture
def input_fasta(tmp_dir):
    """Create a test input FASTA file."""
    fasta_path = os.path.join(tmp_dir, "input.fasta")
    with open(fasta_path, "w") as f:
        f.write(">seq1\n")
        f.write("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n")
        f.write(">seq2\n")
        f.write("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n")
    return fasta_path


@pytest.fixture
def output_fasta(tmp_dir):
    """Create a test output FASTA file with probes."""
    fasta_path = os.path.join(tmp_dir, "output.fa")
    with open(fasta_path, "w") as f:
        f.write(">H5_probe1\n")
        f.write("ATGCATGCATGC\n")
        f.write(">H3_probe2\n")
        f.write("GCTAGCTAGCTA\n")
    return fasta_path


@pytest.fixture
def mock_args():
    """Create a mock args object."""
    args = MagicMock()
    args.blastn = "blastn"
    args.threads = "1"
    args.word_size = "7"
    args.reward = "3"
    args.penalty = "-3"
    args.gapopen = "6"
    args.gapextend = "3"
    args.evalue = "1"
    args.output = "/tmp/test_output"
    return args


@pytest.fixture
def blast_output_content():
    """Mock BLAST output content."""
    return """probe1\tseq1\t100.0\t12\t0\t0\t1\t12\t1\t12\t1e-10\t24.0\tATGCATGCATGC\tATGCATGCATGC
probe2\tseq2\t100.0\t12\t0\t0\t1\t12\t1\t12\t1e-10\t24.0\tGCTAGCTAGCTA\tGCTAGCTAGCTA
probe1\tseq1\t95.0\t12\t1\t0\t1\t12\t5\t16\t1e-8\t22.0\tATGCATGCATGC\tATGCATGCATGC"""


@patch('subprocess.run')
def test_prepare_blast_db(mock_subprocess, tmp_dir, input_fasta):
    """Test BLAST database preparation."""
    mock_subprocess.return_value = MagicMock(returncode=0)
    
    db_path = os.path.join(tmp_dir, "test_db")
    result = prepare_blast_db(input_fasta, db_path)
    
    assert result == db_path
    mock_subprocess.assert_called_once()
    args, kwargs = mock_subprocess.call_args
    assert "makeblastdb" in args[0]
    assert input_fasta in args[0]
    assert db_path in args[0]


def test_get_probes_from_output(output_fasta):
    """Test probe extraction from output FASTA."""
    probes = get_probes_from_output(output_fasta)
    
    assert len(probes) == 2
    assert "probe1" in probes
    assert "probe2" in probes
    # Probes should be converted to RNA (T -> U)
    assert probes["probe1"] == "AUGCAUGCAUGC"  # ATGCATGCATGC -> AUGCAUGCAUGC
    assert probes["probe2"] == "GCUAGCUAGCUA"  # GCTAGCTAGCTA -> GCUAGCUAGCUA


def test_get_probes_from_output_missing_file():
    """Test probe extraction with missing file."""
    with pytest.raises(FileNotFoundError):
        get_probes_from_output("/nonexistent/file.fa")


@patch('subprocess.run')
@patch('tempfile.NamedTemporaryFile')
def test_run_blast_mapping(mock_tempfile, mock_subprocess, mock_args, tmp_dir, input_fasta):
    """Test BLAST mapping execution."""
    # Setup mocks
    mock_file = MagicMock()
    mock_file.name = os.path.join(tmp_dir, "tmp_probes.fa")
    mock_file.__enter__ = MagicMock(return_value=mock_file)
    mock_file.__exit__ = MagicMock(return_value=None)
    mock_tempfile.return_value = mock_file
    
    mock_subprocess.return_value = MagicMock(returncode=0)
    
    probes = {"probe1": "ATGCATGCATGC", "probe2": "GCTAGCTAGCTA"}
    db_path = os.path.join(tmp_dir, "test_db")
    output_file = os.path.join(tmp_dir, "blast_output.tsv")
    
    result = run_blast_mapping(probes, db_path, mock_args, output_file)
    
    assert result == output_file
    mock_subprocess.assert_called_once()
    args, kwargs = mock_subprocess.call_args
    assert "blastn" in args[0] or args[0][0] == "blastn"


def test_extract_extended_sequence_forward(input_fasta):
    """Test extended sequence extraction (forward strand)."""
    seq = extract_extended_sequence(input_fasta, "seq1", 1, 12, extend=10)
    
    assert seq is not None
    assert len(seq) >= 12  # Should have at least the original sequence
    assert "ATGCATGCATGC" in seq


def test_extract_extended_sequence_reverse(input_fasta):
    """Test extended sequence extraction (reverse strand)."""
    # For reverse complement, send < sstart
    seq = extract_extended_sequence(input_fasta, "seq1", 12, 1, extend=10)
    
    assert seq is not None
    assert len(seq) >= 12


def test_extract_extended_sequence_not_found(input_fasta):
    """Test extended sequence extraction with non-existent sequence."""
    seq = extract_extended_sequence(input_fasta, "nonexistent", 1, 12, extend=10)
    
    assert seq is None


def test_extract_best_hits(tmp_dir, input_fasta, blast_output_content):
    """Test best hit extraction from BLAST results."""
    blast_output = os.path.join(tmp_dir, "blast_output.tsv")
    with open(blast_output, "w") as f:
        f.write(blast_output_content)
    
    results = extract_best_hits(blast_output, input_fasta, extend_nuc=10)
    
    assert isinstance(results, pd.DataFrame)
    assert len(results) == 2  # Should have one row per unique probe
    expected_columns = [
        "probe_seq", "target_seq", "target_seq_reverse_complement",
        "vienna_rna_mfe", "dna_dna_duplex_dg", "dna_rna_duplex_dg", "rna_rna_duplex_dg"
    ]
    for col in expected_columns:
        assert col in results.columns
    assert "probe1" in results["probe_seq"].values or "AUGCAUGCAUGC" in results["probe_seq"].values
    # Verify reverse complement column is present and not empty
    assert all(results["target_seq_reverse_complement"].str.len() > 0)


def test_extract_best_hits_empty_file(tmp_dir, input_fasta):
    """Test best hit extraction with empty BLAST output."""
    blast_output = os.path.join(tmp_dir, "empty_blast.tsv")
    with open(blast_output, "w") as f:
        pass  # Empty file
    
    results = extract_best_hits(blast_output, input_fasta, extend_nuc=10)
    
    assert isinstance(results, pd.DataFrame)
    assert len(results) == 0
    expected_columns = [
        "probe_seq", "target_seq", "target_seq_reverse_complement",
        "vienna_rna_mfe", "dna_dna_duplex_dg", "dna_rna_duplex_dg", "rna_rna_duplex_dg"
    ]
    assert list(results.columns) == expected_columns


def test_extract_best_hits_nonexistent_file(tmp_dir, input_fasta):
    """Test best hit extraction with non-existent file."""
    blast_output = os.path.join(tmp_dir, "nonexistent.tsv")
    
    results = extract_best_hits(blast_output, input_fasta, extend_nuc=10)
    
    assert isinstance(results, pd.DataFrame)
    assert len(results) == 0
    expected_columns = [
        "probe_seq", "target_seq", "target_seq_reverse_complement",
        "vienna_rna_mfe", "dna_dna_duplex_dg", "dna_rna_duplex_dg", "rna_rna_duplex_dg"
    ]
    assert list(results.columns) == expected_columns


@patch('PROBESt.modeling.create_visualizations_from_table')
@patch('PROBESt.modeling.prepare_blast_db')
@patch('PROBESt.modeling.get_probes_from_output')
@patch('PROBESt.modeling.run_blast_mapping')
@patch('PROBESt.modeling.extract_best_hits')
@patch('pandas.DataFrame.to_csv')
def test_run_modeling(
    mock_to_csv,
    mock_extract_best_hits,
    mock_run_blast_mapping,
    mock_get_probes,
    mock_prepare_db,
    mock_create_visualizations,
    mock_args,
    tmp_dir,
    input_fasta,
    output_fasta
):
    """Test the complete modeling pipeline."""
    # Setup mocks
    mock_prepare_db.return_value = os.path.join(tmp_dir, "test_db")
    mock_get_probes.return_value = {"probe1": "ATGCATGCATGC", "probe2": "GCTAGCTAGCTA"}
    mock_run_blast_mapping.return_value = os.path.join(tmp_dir, "blast_output.tsv")
    mock_create_visualizations.return_value = os.path.join(tmp_dir, "visualizations")
    
    mock_df = pd.DataFrame({
        "probe_seq": ["ATGCATGCATGC", "GCTAGCTAGCTA"],
        "target_seq": ["ATGCATGCATGCATGCATGCATGC", "GCTAGCTAGCTAGCTAGCTAGCTA"]
    })
    mock_extract_best_hits.return_value = mock_df
    
    output_table = os.path.join(tmp_dir, "modeling_results.tsv")
    result = run_modeling(mock_args, input_fasta, output_fasta, output_table)
    
    assert result == output_table
    mock_prepare_db.assert_called_once()
    mock_get_probes.assert_called_once_with(output_fasta)
    mock_run_blast_mapping.assert_called_once()
    mock_extract_best_hits.assert_called_once()
    mock_to_csv.assert_called_once()
    mock_create_visualizations.assert_called_once()


def test_integration_modeling(tmp_dir, input_fasta, output_fasta, mock_args):
    """Integration test for modeling module (requires BLAST tools)."""
    # Skip if BLAST tools are not available
    try:
        subprocess.run(["makeblastdb", "-version"], capture_output=True, check=True)
        subprocess.run(["blastn", "-version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        pytest.skip("BLAST tools not available")
    
    # Create output directory
    os.makedirs(mock_args.output, exist_ok=True)
    os.makedirs(os.path.join(mock_args.output, ".tmp"), exist_ok=True)
    
    # Prepare database
    db_path = os.path.join(tmp_dir, "test_db")
    prepare_blast_db(input_fasta, db_path)
    
    # Get probes
    probes = get_probes_from_output(output_fasta)
    
    # Run BLAST
    blast_output = os.path.join(tmp_dir, "blast_output.tsv")
    run_blast_mapping(probes, db_path, mock_args, blast_output)
    
    # Extract best hits
    if os.path.exists(blast_output) and os.path.getsize(blast_output) > 0:
        results = extract_best_hits(blast_output, input_fasta, extend_nuc=10)
        
        assert isinstance(results, pd.DataFrame)
        assert len(results) > 0
        assert "probe_seq" in results.columns
        assert "target_seq" in results.columns

