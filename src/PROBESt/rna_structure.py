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


import RNA
from typing import Union, Tuple, Optional

def validate_sequence(sequence: str, seq_type: str = "RNA") -> None:
    """
    Validate sequence for RNA/DNA characters and type consistency.
    
    Args:
        sequence (str): Sequence to validate
        seq_type (str): Type of sequence ("RNA" or "DNA")
        
    Raises:
        ValueError: If sequence contains invalid characters, is empty, or has type mismatches
    """
    if not sequence:
        raise ValueError("Sequence cannot be empty")
    
    # Check for type mismatches
    seq = sequence.upper()
    
    # Check for other invalid characters
    valid_chars = set("ATGCUIN")
    invalid_chars = set(seq) - valid_chars
    
    if invalid_chars:
        raise ValueError(f"Invalid characters in {seq_type} sequence: {invalid_chars}")

def get_reverse_complement(sequence: str, seq_type: str = "RNA") -> str:
    """
    Get the reverse complement of a sequence.
    
    Args:
        sequence (str): Input sequence
        seq_type (str): Type of sequence ("RNA" or "DNA")
        
    Returns:
        str: Reverse complement sequence
        
    Raises:
        ValueError: If sequence contains invalid characters or is empty
    """
    # Validate input
    validate_sequence(sequence, seq_type)
    
    # Define complement mapping
    complement_map = {
        'A': 'U' if seq_type == "RNA" else 'T',
        'T': 'A',
        'U': 'A',
        'G': 'C',
        'C': 'G'
    }
    
    # Convert to uppercase and get complement
    seq = sequence.upper()
    complement = ''.join(complement_map[base] for base in seq)
    
    # Reverse the sequence
    return complement[::-1]

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

def calculate_dimer_G(string1: str, string2: Optional[str] = None, type1: str = "RNA", type2: Optional[str] = None) -> float:
    """
    Calculate the Gibbs free energy of dimer formation between two sequences.
    If string2 is None, calculates dimer with reverse complement of string1.
    
    Args:
        string1 (str): First sequence
        string2 (str, optional): Second sequence. If None, uses reverse complement of string1
        type1 (str): Type of first sequence ("RNA" or "DNA")
        type2 (str, optional): Type of second sequence ("RNA" or "DNA"). If None, uses type1
        
    Returns:
        float: Gibbs free energy in kcal/mol
        
    Raises:
        ValueError: If sequences contain invalid characters or is empty
    """
    # Validate first sequence
    validate_sequence(string1, type1)
    
    # Handle reverse complement case
    if string2 is None:
        string2 = get_reverse_complement(string1, type1)
        type2 = type1
    elif type2 is None:
        type2 = type1
    
    # Validate second sequence
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