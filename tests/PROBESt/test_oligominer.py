import pytest
import os
import sys
import tempfile
import shutil
from Bio import SeqIO
from src.PROBESt.oligominer import parse_oligominer_bed, extract_sequence_from_fasta, oligominer2fasta
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
    # Create a sample FASTA file
    fasta_file = tmp_path / "test_sequence.fasta"
    with open(fasta_file, "w") as f:
        f.write(">seq1\nATGCATGCATGCATGCATGCATGCATGC\n>seq2\nCGTAACGTAACGTAACGTAACGTAACG\n")
    
    # Create a sample BED file (OligoMiner format)
    bed_file = tmp_path / "probes.bed"
    with open(bed_file, "w") as f:
        # BED format: chrom, start, end, name, score, strand
        f.write("seq1\t0\t25\tprobe1\t100\t+\n")
        f.write("seq1\t5\t30\tprobe2\t100\t+\n")
        f.write("seq2\t0\t25\tprobe3\t100\t+\n")
    
    return {
        "fasta_file": str(fasta_file),
        "bed_file": str(bed_file),
        "tmp_path": tmp_path
    }

def test_parse_oligominer_bed(test_dir):
    """Test the parse_oligominer_bed function."""
    # Parse BED file without FASTA (sequences will be empty)
    result = parse_oligominer_bed(test_dir["bed_file"])
    
    # Expected output (without sequences)
    expected_output = [
        ("seq1", "0_25", "PROBE", ""),
        ("seq1", "5_30", "PROBE", ""),
        ("seq2", "0_25", "PROBE", ""),
    ]
    
    # Verify the parsed output matches (sequences will be empty)
    assert len(result) == len(expected_output)
    for i, (exp, res) in enumerate(zip(expected_output, result)):
        assert res[0] == exp[0], f"Sequence ID mismatch at index {i}"
        assert res[1] == exp[1], f"Probe number mismatch at index {i}"
        assert res[2] == exp[2], f"Probe type mismatch at index {i}"

def test_parse_oligominer_bed_with_fasta(test_dir):
    """Test the parse_oligominer_bed function with FASTA file for sequence extraction."""
    # Parse BED file with FASTA
    result = parse_oligominer_bed(test_dir["bed_file"], test_dir["fasta_file"])
    
    # Verify sequences are extracted
    assert len(result) == 3
    # seq1 is "ATGCATGCATGCATGCATGCATGCATGC" (28 chars), [0:25] gives first 25 chars
    assert result[0][3] == "ATGCATGCATGCATGCATGCATGCA"  # seq1[0:25] = 25 chars
    assert len(result[0][3]) == 25, f"Expected 25 chars, got {len(result[0][3])}"
    # seq1[5:30] gives chars 5-28 (seq1 is only 28 chars) = 23 chars
    assert result[1][3] == "TGCATGCATGCATGCATGCATGC"  # seq1[5:30] = 23 chars (limited by seq length)
    assert len(result[1][3]) == 23, f"Expected 23 chars, got {len(result[1][3])}"
    # seq2 is "CGTAACGTAACGTAACGTAACGTAACG" (27 chars), [0:25] gives first 25 chars
    assert result[2][3] == "CGTAACGTAACGTAACGTAACGTAA"  # seq2[0:25] = 25 chars
    assert len(result[2][3]) == 25, f"Expected 25 chars, got {len(result[2][3])}"

def test_extract_sequence_from_fasta(test_dir):
    """Test the extract_sequence_from_fasta function."""
    # Extract sequence from FASTA
    result = extract_sequence_from_fasta(test_dir["fasta_file"], "seq1", 0, 25)
    
    # Expected output: seq1 is "ATGCATGCATGCATGCATGCATGCATGC" (28 chars), [0:25] = 25 chars
    expected = "ATGCATGCATGCATGCATGCATGCA"
    
    assert result == expected, f"Extracted sequence does not match. Got '{result}' (len={len(result)}), expected '{expected}' (len={len(expected)})"
    assert len(result) == 25, f"Expected 25 characters, got {len(result)}"
    
    # Test with different coordinates: seq1[5:30] = 23 chars (seq1 is only 28 chars, so [5:30] becomes [5:28])
    result2 = extract_sequence_from_fasta(test_dir["fasta_file"], "seq1", 5, 30)
    expected2 = "TGCATGCATGCATGCATGCATGC"
    assert result2 == expected2, f"Extracted sequence with offset does not match. Got '{result2}' (len={len(result2)}), expected '{expected2}' (len={len(expected2)})"
    assert len(result2) == 23, f"Expected 23 characters (seq1 is only 28 chars), got {len(result2)}"

def test_oligominer2fasta(test_dir, mock_args):
    """Test the oligominer2fasta function."""
    # Set add_set to None
    mock_args.add_set = None
    
    # Create probes list
    probes = [
        ("seq1", "0_25", "PROBE", "ATGCATGCATGCATGCATGCATGCA"),
        ("seq2", "0_25", "PROBE", "CGTAACGTAACGTAACGTAACGTAA"),
    ]
    
    # Write probes to FASTA
    oligominer2fasta(mock_args, str(test_dir["tmp_path"]), probes, test_dir["fasta_file"])
    
    # Verify output FASTA file was created
    output_fasta = os.path.join(test_dir["tmp_path"], "output.fa")
    assert os.path.exists(output_fasta), "Output FASTA file was not created"
    
    # Read and verify contents
    with open(output_fasta, "r") as f:
        content = f.read()
        assert ">seq1_0_25_PROBE" in content
        assert "ATGCATGCATGCATGCATGCATGCA" in content
        assert ">seq2_0_25_PROBE" in content
        assert "CGTAACGTAACGTAACGTAACGTAA" in content  # seq2[0:25] = 25 chars

def test_oligominer2fasta_with_add_set(test_dir, mock_args):
    """Test the oligominer2fasta function with additional sequences."""
    # Create additional sequences file
    add_set_file = test_dir["tmp_path"] / "additional.fasta"
    with open(add_set_file, "w") as f:
        f.write(">extra_seq\nAAAAATTTTT\n")
    
    mock_args.add_set = [str(add_set_file)]
    
    # Create probes list
    probes = [
        ("seq1", "0_25", "PROBE", "ATGCATGCATGCATGCATGCATGC"),
    ]
    
    # Write probes to FASTA
    oligominer2fasta(mock_args, str(test_dir["tmp_path"]), probes, test_dir["fasta_file"])
    
    # Verify output FASTA file contains both original and additional sequences
    output_fasta = os.path.join(test_dir["tmp_path"], "output.fa")
    with open(output_fasta, "r") as f:
        content = f.read()
        assert ">seq1_0_25_PROBE" in content
        assert ">extra_seq" in content
        assert "AAAAATTTTT" in content

def test_parse_oligominer_bed_empty_file(tmp_path):
    """Test parsing an empty BED file."""
    # Create empty BED file
    bed_file = tmp_path / "empty.bed"
    bed_file.touch()
    
    result = parse_oligominer_bed(str(bed_file))
    assert result == [], "Empty BED file should return empty list"

def test_parse_oligominer_bed_nonexistent_file():
    """Test parsing a non-existent BED file."""
    result = parse_oligominer_bed("/nonexistent/path/to/file.bed")
    assert result == [], "Non-existent BED file should return empty list"

def test_extract_sequence_from_fasta_nonexistent():
    """Test extracting sequence from non-existent FASTA file."""
    import pytest
    # The function will raise FileNotFoundError when trying to parse non-existent file
    # We should handle this gracefully or let it raise
    try:
        result = extract_sequence_from_fasta("/nonexistent/file.fasta", "seq1", 0, 25)
        # If it doesn't raise, it should return empty string
        assert result == "", "Non-existent FASTA should return empty string or raise error"
    except FileNotFoundError:
        # This is also acceptable behavior
        pass

def test_extract_sequence_from_fasta_invalid_coords(test_dir):
    """Test extracting sequence with invalid coordinates."""
    # Test with coordinates beyond sequence length
    result = extract_sequence_from_fasta(test_dir["fasta_file"], "seq1", 0, 1000)
    # Should return available sequence (up to end)
    assert len(result) <= 30, "Should not exceed sequence length"
    
    # Test with invalid sequence ID
    result2 = extract_sequence_from_fasta(test_dir["fasta_file"], "nonexistent", 0, 25)
    assert result2 == "", "Invalid sequence ID should return empty string"


