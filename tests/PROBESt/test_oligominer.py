import pytest
import os
import sys
import tempfile
from Bio import SeqIO
from src.PROBESt.oligominer import fastq_to_fasta
from src.PROBESt.args import arguments_parse

@pytest.fixture
def mock_args(monkeypatch):
    """Mock command line arguments for testing."""
    test_args = [
        "script.py",
        "-i", "input.fasta",
        "-tb", "true_base",
        "-fb", "false_base",
        "-c", "contig_table",
        "-o", "output"
    ]
    monkeypatch.setattr(sys, 'argv', test_args)
    return arguments_parse()

@pytest.fixture
def test_dir(tmp_path):
    """Create a temporary test directory with test files."""
    # Create a sample FASTQ file
    # Note: Quality scores must match sequence length exactly
    fastq_file = tmp_path / "test_probes.fastq"
    with open(fastq_file, "w") as f:
        # seq1:0-25 = 25 chars, quality = 25 chars
        seq1_25 = "ATGCATGCATGCATGCATGCATGCA"
        f.write(f"@seq1:0-25\n{seq1_25}\n+\n{'~' * len(seq1_25)}\n")
        # seq1:5-30 = 23 chars (seq1 is only 28 chars, so 5-30 becomes 5-28 = 23), quality = 23 chars
        seq1_23 = "TGCATGCATGCATGCATGCATGC"
        f.write(f"@seq1:5-30\n{seq1_23}\n+\n{'~' * len(seq1_23)}\n")
        # seq2:0-25 = 25 chars, quality = 25 chars
        seq2_25 = "CGTAACGTAACGTAACGTAACGTAA"
        f.write(f"@seq2:0-25\n{seq2_25}\n+\n{'~' * len(seq2_25)}\n")
    
    return {
        "fastq_file": str(fastq_file),
        "tmp_path": tmp_path
    }

def test_fastq_to_fasta(test_dir):
    """Test the fastq_to_fasta function."""
    output_fasta = os.path.join(test_dir["tmp_path"], "output.fasta")
    fastq_to_fasta(test_dir["fastq_file"], output_fasta)
    
    # Verify output FASTA file was created
    assert os.path.exists(output_fasta), "Output FASTA file was not created"
    
    # Read and verify contents
    records = list(SeqIO.parse(output_fasta, "fasta"))
    assert len(records) == 3, f"Expected 3 sequences, got {len(records)}"
    
    # Check first sequence
    assert records[0].id == "seq1:0-25"
    assert str(records[0].seq) == "ATGCATGCATGCATGCATGCATGCA"
    
    # Check second sequence
    assert records[1].id == "seq1:5-30"
    assert str(records[1].seq) == "TGCATGCATGCATGCATGCATGC"
    
    # Check third sequence
    assert records[2].id == "seq2:0-25"
    assert str(records[2].seq) == "CGTAACGTAACGTAACGTAACGTAA"

def test_fastq_to_fasta_empty_file(tmp_path):
    """Test converting an empty FASTQ file."""
    # Create empty FASTQ file
    fastq_file = tmp_path / "empty.fastq"
    fastq_file.touch()
    
    output_fasta = tmp_path / "output.fasta"
    fastq_to_fasta(str(fastq_file), str(output_fasta))
    
    # Output file should exist but be empty or have minimal content
    assert os.path.exists(output_fasta), "Output FASTA file should be created even for empty input"

def test_fastq_to_fasta_nonexistent_file(tmp_path):
    """Test converting a non-existent FASTQ file."""
    output_fasta = tmp_path / "output.fasta"
    
    # Should raise FileNotFoundError
    with pytest.raises(FileNotFoundError):
        fastq_to_fasta("/nonexistent/file.fastq", str(output_fasta))
