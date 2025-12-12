"""Tests for the de-degeneration module."""

import pytest
import os
import random
from src.PROBESt.evolution import dedegenerate_position
from src.PROBESt.dedegeneration import (
    count_degenerate_nucleotides,
    read_fasta_sequences
)


def test_count_degenerate_nucleotides():
    """Test counting degenerate nucleotides in sequences."""
    # Test with no degenerate nucleotides
    assert count_degenerate_nucleotides("ATGC") == 0
    assert count_degenerate_nucleotides("ATGCATGC") == 0
    
    # Test with degenerate nucleotides
    assert count_degenerate_nucleotides("ATNG") == 1
    assert count_degenerate_nucleotides("NRY") == 3
    assert count_degenerate_nucleotides("ATNCN") == 2
    
    # Test with mixed case
    assert count_degenerate_nucleotides("AtNc") == 1
    assert count_degenerate_nucleotides("nry") == 3
    
    # Test with all degenerate types
    assert count_degenerate_nucleotides("NWSMKRYBDHVN") == 12
    
    # Test with empty string
    assert count_degenerate_nucleotides("") == 0


def test_dedegenerate_position():
    """Test de-degenerating individual nucleotide positions."""
    random.seed(42)  # Set seed for reproducibility
    
    # Test with non-degenerate nucleotide (should return as-is)
    assert dedegenerate_position("A", 1.0) == "A"
    assert dedegenerate_position("T", 1.0) == "T"
    assert dedegenerate_position("G", 1.0) == "G"
    assert dedegenerate_position("C", 1.0) == "C"
    
    # Test with degenerate nucleotide at rate 1.0 (should always de-degenerate)
    result_n = dedegenerate_position("N", 1.0)
    assert result_n in ["A", "T", "G", "C"]
    
    result_r = dedegenerate_position("R", 1.0)
    assert result_r in ["A", "G"]
    
    result_y = dedegenerate_position("Y", 1.0)
    assert result_y in ["C", "T"]
    
    # Test with degenerate nucleotide at rate 0.0 (should never de-degenerate)
    assert dedegenerate_position("N", 0.0) == "N"
    assert dedegenerate_position("R", 0.0) == "R"
    
    # Test case insensitivity
    result_n_lower = dedegenerate_position("n", 1.0)
    assert result_n_lower in ["A", "T", "G", "C"]
    
    # Test all degenerate codes
    degenerate_codes = {
        'N': ['A', 'T', 'G', 'C'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G']
    }
    
    for code, possible_values in degenerate_codes.items():
        result = dedegenerate_position(code, 1.0)
        assert result in possible_values, f"Failed for {code}: got {result}, expected one of {possible_values}"


def test_read_fasta_sequences(tmp_path):
    """Test reading sequences from FASTA files."""
    # Test FISH algorithm
    fasta_fish = tmp_path / "test_fish.fasta"
    with open(fasta_fish, "w") as f:
        f.write(">seq1\nATGC\n>seq2\nCGTA\n>seq3\nNNNN\n")
    
    seqs = read_fasta_sequences(str(fasta_fish), algorithm='FISH')
    assert len(seqs) == 3
    assert seqs["seq1"] == "ATGC"
    assert seqs["seq2"] == "CGTA"
    assert seqs["seq3"] == "NNNN"
    
    # Test primer algorithm
    fasta_primer = tmp_path / "test_primer.fasta"
    with open(fasta_primer, "w") as f:
        f.write(">seq1_LEFT\nATGC\n>seq1_RIGHT\nCGTA\n")
    
    seqs = read_fasta_sequences(str(fasta_primer), algorithm='primer')
    assert len(seqs) == 2
    assert seqs["seq1_LEFT"] == "ATGC"
    assert seqs["seq1_RIGHT"] == "CGTA"


def test_read_fasta_sequences_empty(tmp_path):
    """Test reading empty FASTA file."""
    fasta_empty = tmp_path / "empty.fasta"
    with open(fasta_empty, "w") as f:
        f.write("")
    
    seqs = read_fasta_sequences(str(fasta_empty), algorithm='FISH')
    assert len(seqs) == 0


def test_dedegenerate_sequence_integration():
    """Test de-degenerating a full sequence."""
    random.seed(42)
    
    # Create a sequence with degenerate nucleotides
    seq = "ATNGC"
    
    # De-degenerate with rate 1.0 (all degenerate should be replaced)
    result = "".join([dedegenerate_position(_, dedegen_rate=1.0) for _ in seq])
    
    # Should have same length
    assert len(result) == len(seq)
    
    # First three should be unchanged (non-degenerate)
    assert result[0] == "A"
    assert result[1] == "T"
    assert result[4] == "C"
    
    # Fourth position was N, should be replaced with A, T, G, or C
    assert result[3] in ["A", "T", "G", "C"]


def test_count_vs_dedegenerate():
    """Test that de-degeneration reduces degenerate nucleotide count."""
    random.seed(42)
    
    # Sequence with many degenerate nucleotides
    seq = "NNNYYRRR"
    initial_count = count_degenerate_nucleotides(seq)
    assert initial_count == 8
    
    # De-degenerate with rate 1.0
    dedegened = "".join([dedegenerate_position(_, dedegen_rate=1.0) for _ in seq])
    final_count = count_degenerate_nucleotides(dedegened)
    
    # Should have fewer or equal degenerate nucleotides
    assert final_count <= initial_count
    # With rate 1.0, should have 0
    assert final_count == 0


def test_dedegenerate_position_probabilistic():
    """Test that de-degeneration rate works probabilistically."""
    random.seed(42)
    
    # Run many times with rate 0.5
    results = []
    for _ in range(100):
        result = dedegenerate_position("N", 0.5)
        results.append(result)
    
    # Some should be de-degenerated (ATGC) and some should remain N
    degen_count = sum(1 for r in results if r == "N")
    dedegened_count = sum(1 for r in results if r in ["A", "T", "G", "C"])
    
    # Both should happen (probabilistic)
    assert degen_count > 0
    assert dedegened_count > 0
    assert degen_count + dedegened_count == 100

