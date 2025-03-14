import pytest
from Bio.SeqUtils import MeltingTemp as mt
from PROBESt.probe_alignment_profiler import (
    calculate_gc_content,
    calculate_melting_temperature,
    process_sequence_pair
)


def test_gc_content_basic():
    """Test GC content calculation for basic sequences."""
    assert calculate_gc_content("GCGC") == 100.0
    assert calculate_gc_content("ATAT") == 0.0
    assert calculate_gc_content("GATC") == 50.0
    assert calculate_gc_content("") == 0.0


def test_gc_content_iupac():
    """Test GC content calculation with IUPAC ambiguous bases."""
    # S is always G or C
    assert calculate_gc_content("SSSS") == 100.0

    # R (A/G) and Y (C/T) are 50% GC
    assert calculate_gc_content("RRRR") == 50.0
    assert calculate_gc_content("YYYY") == 50.0

    # K (G/T) and M (A/C) are 50% GC
    assert calculate_gc_content("KKKK") == 50.0
    assert calculate_gc_content("MMMM") == 50.0

    # B (C/G/T) and V (A/C/G) are 66.7% GC
    assert pytest.approx(calculate_gc_content("BBBB"), 0.1) == 66.7
    assert pytest.approx(calculate_gc_content("VVVV"), 0.1) == 66.7

    # D (A/G/T) and H (A/C/T) are 33.3% GC
    assert pytest.approx(calculate_gc_content("DDDD"), 0.1) == 33.3
    assert pytest.approx(calculate_gc_content("HHHH"), 0.1) == 33.3

    # N (any base) is 50% GC
    assert calculate_gc_content("NNNN") == 50.0


def test_gc_content_mixed():
    """Test GC content calculation with mixed regular and IUPAC bases."""
    # "GC" + "SS" (all GC) = 100%
    assert calculate_gc_content("GCSS") == 100.0

    # "AT" (0% GC) + "RY" (50% GC each) = 25%
    assert calculate_gc_content("ATRY") == 25.0

    # Complex mixture
    # "GC" (100%) + "AT" (0%) + "NN" (50%) = 50%
    assert calculate_gc_content("GCATNN") == 50.0


def test_melting_temperature_basic():
    """Test melting temperature calculation using Biopython's nearest-neighbor method."""
    # Test with some standard sequences
    seq1 = "GCGCGC"
    seq2 = "ATATATA"
    seq3 = "GCATGCAT"

    # Using pytest.approx because Tm calculations can have small floating-point differences
    assert pytest.approx(calculate_melting_temperature(seq1), 0.01) == mt.Tm_NN(seq1, strict=False)
    assert pytest.approx(calculate_melting_temperature(seq2), 0.01) == mt.Tm_NN(seq2, strict=False)
    assert pytest.approx(calculate_melting_temperature(seq3), 0.01) == mt.Tm_NN(seq3, strict=False)


def test_melting_temperature_empty():
    """Test melting temperature calculation for empty sequence."""
    assert calculate_melting_temperature("") == 0


def test_melting_temperature_iupac():
    """Test melting temperature calculation with IUPAC ambiguous bases."""
    # Test that the function handles IUPAC codes without raising exceptions
    ambiguous_seq = "GCRYSWKMBDHVN"
    try:
        tm = calculate_melting_temperature(ambiguous_seq)
        assert isinstance(tm, float)
    except Exception as e:
        pytest.fail(f"calculate_melting_temperature failed with IUPAC sequence: {str(e)}")


def test_process_sequence_pair():
    """Test processing of sequence pairs."""
    probe = "GCTAGCAT"
    target = "GCTGGCAT"

    result = process_sequence_pair(probe, target)

    # Test presence of all expected keys
    expected_keys = {
        'probe_sequence', 'target_sequence',
        'probe_gc_content', 'target_gc_content',
        'probe_tm', 'target_tm',
        'probe_length', 'target_length'
    }
    assert set(result.keys()) == expected_keys

    # Test specific values
    assert result['probe_sequence'] == probe
    assert result['target_sequence'] == target
    assert result['probe_gc_content'] == 50.0
    assert result['target_gc_content'] == 62.5
    assert result['probe_length'] == 8
    assert result['target_length'] == 8

    # Test that Tm values are reasonable (should be positive numbers for these sequences)
    assert result['probe_tm'] > 0
    assert result['target_tm'] > 0

    # Test that Tm values match Biopython's calculations
    assert pytest.approx(result['probe_tm'], 0.01) == mt.Tm_NN(probe, strict=False)
    assert pytest.approx(result['target_tm'], 0.01) == mt.Tm_NN(target, strict=False)
