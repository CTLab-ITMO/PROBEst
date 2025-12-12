import pytest
from src.PROBESt.rna_structure import calculate_hairpin_prob, calculate_dimer_G, get_reverse_complement

def test_calculate_hairpin_prob():
    # Test with a simple hairpin sequence
    seq = "GGGAAACCC"
    prob = calculate_hairpin_prob(seq)
    assert isinstance(prob, float)
    assert 0 <= prob <= 1
    
    # Test with a more complex sequence
    seq = "AUGCAUGCAUGC"
    prob = calculate_hairpin_prob(seq)
    assert isinstance(prob, float)
    assert 0 <= prob <= 1

def test_get_reverse_complement():
    # Test RNA reverse complement
    seq = "AUGCAUGC"
    rc = get_reverse_complement(seq, "RNA")
    assert rc == "GCAUGCAU"
    
    # Test DNA reverse complement
    seq = "ATGCATGC"
    rc = get_reverse_complement(seq, "DNA")
    assert rc == "GCATGCAT"
    
    # Test with invalid characters
    with pytest.raises(ValueError):
        get_reverse_complement("ATGCXYZ", "DNA")
    
    # Test with empty sequence
    with pytest.raises(ValueError):
        get_reverse_complement("", "DNA")

def test_calculate_dimer_G():
    # Test RNA-RNA dimer
    seq1 = "AUGCAUGCAUGC"
    seq2 = "UACGUACGUACG"
    energy = calculate_dimer_G(seq1, seq2, type1="RNA", type2="RNA")
    assert isinstance(energy, float)
    
    # Test DNA-DNA dimer
    seq1 = "ATGCATGCATGC"
    seq2 = "TACGTACGTACG"
    energy = calculate_dimer_G(seq1, seq2, type1="DNA", type2="DNA")
    assert isinstance(energy, float)
    
    # Test RNA-DNA dimer
    seq1 = "AUGCAUGCAUGC"
    seq2 = "TACGTACGTACG"
    energy = calculate_dimer_G(seq1, seq2, type1="RNA", type2="DNA")
    assert isinstance(energy, float)
    
    # Test with None string2 (reverse complement)
    seq1 = "AUGCAUGCAUGC"
    energy = calculate_dimer_G(seq1, type1="RNA")
    assert isinstance(energy, float)
    
    # Test with None string2 and type2
    seq1 = "ATGCATGCATGC"
    energy = calculate_dimer_G(seq1, type1="DNA")
    assert isinstance(energy, float)

def test_input_validation():
    # Test with invalid characters
    with pytest.raises(ValueError):
        calculate_hairpin_prob("ATGCXYZ")
    
    with pytest.raises(ValueError):
        calculate_dimer_G("ATGCXYZ", "TACGTAC", type1="DNA", type2="DNA")
    
    # Test with empty sequences
    with pytest.raises(ValueError, match="Sequence cannot be empty"):
        calculate_hairpin_prob("")
    
    with pytest.raises(ValueError, match="Sequence cannot be empty"):
        calculate_dimer_G("", "TACGTAC")