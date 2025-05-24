import RNA
from typing import Union, Tuple

def validate_sequence(sequence: str, seq_type: str = "RNA") -> None:
    """
    Validate sequence for RNA/DNA characters.
    
    Args:
        sequence (str): Sequence to validate
        seq_type (str): Type of sequence ("RNA" or "DNA")
        
    Raises:
        ValueError: If sequence contains invalid characters or is empty
    """
    if not sequence:
        raise ValueError("Sequence cannot be empty")
    
    valid_chars = set("ATGCU" if seq_type == "RNA" else "ATGC")
    invalid_chars = set(sequence.upper()) - valid_chars
    
    if invalid_chars:
        raise ValueError(f"Invalid characters in {seq_type} sequence: {invalid_chars}")

def calculate_hairpin_prob(sequence: str) -> float:
    """
    Calculate the probability of hairpin formation in an RNA sequence.
    
    Args:
        sequence (str): RNA sequence (ATGCU)
        
    Returns:
        float: Probability of the most stable structure
        
    Raises:
        ValueError: If sequence contains invalid characters or is empty
    """
    # Validate input
    validate_sequence(sequence, "RNA")
    
    # Convert sequence to uppercase
    seq = sequence.upper()
    
    # Create model details
    md = RNA.md()
    md.noLP = 1  # No lonely pairs
    
    # Calculate partition function
    fc = RNA.fold_compound(seq, md)
    (ss, mfe) = fc.mfe()
    
    # Calculate structure probability
    fc.pf()
    prob = fc.pr_structure(ss)
    
    return prob

def calculate_dimer_G(string1: str, string2: str, type1: str = "RNA", type2: str = "DNA") -> float:
    """
    Calculate the Gibbs free energy of dimer formation between two sequences.
    
    Args:
        string1 (str): First sequence
        string2 (str): Second sequence
        type1 (str): Type of first sequence ("RNA" or "DNA")
        type2 (str): Type of second sequence ("RNA" or "DNA")
        
    Returns:
        float: Gibbs free energy in kcal/mol
        
    Raises:
        ValueError: If sequences contain invalid characters or are empty
    """
    # Validate inputs
    validate_sequence(string1, type1)
    validate_sequence(string2, type2)
    
    # Convert sequences to uppercase
    seq1 = string1.upper()
    seq2 = string2.upper()
    
    # Convert DNA to RNA if needed
    if type1 == "DNA":
        seq1 = seq1.replace("T", "U")
    if type2 == "DNA":
        seq2 = seq2.replace("T", "U")
    
    # Combine sequences with '&' for RNAcofold
    combined_seq = f"{seq1}&{seq2}"
    
    # Calculate cofolding energy
    (ss, energy) = RNA.cofold(combined_seq)
    
    return energy 