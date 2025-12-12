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


import pytest
import pandas as pd
import numpy as np
import os
import tempfile
from PROBESt.misc import (
    calculate_gc_content, 
    get_sequence_from_fasta,
    load_sequences_from_fasta,
    extend_blast_output_with_parameters
)


def test_calculate_gc_content():
    """Test GC content calculation."""
    # Test with GC-rich sequence
    assert calculate_gc_content("GCGCGC") == 100.0
    
    # Test with AT-rich sequence
    assert calculate_gc_content("ATATAT") == 0.0
    
    # Test with mixed sequence
    assert abs(calculate_gc_content("ATGC") - 50.0) < 0.1
    
    # Test with empty sequence
    assert calculate_gc_content("") == 0.0
    
    # Test with RNA sequence
    assert abs(calculate_gc_content("AUGCAUGC") - 50.0) < 0.1


def test_get_sequence_from_fasta(tmp_path):
    """Test sequence extraction from FASTA file."""
    # Create a test FASTA file
    fasta_file = tmp_path / "test.fa"
    with open(fasta_file, 'w') as f:
        f.write(">probe1\nATGCATGC\n")
        f.write(">probe2\nGCGCGCGC\n")
        f.write(">H1_probe3\nTTTTTTTT\n")
        f.write(">probe4_LEFT\nAAAA\n")
    
    # Test normal probe name
    seq = get_sequence_from_fasta("probe1", str(fasta_file))
    assert seq == "ATGCATGC"
    
    # Test probe with H prefix
    seq = get_sequence_from_fasta("probe3", str(fasta_file))
    assert seq == "TTTTTTTT"
    
    # Test probe with LEFT suffix
    seq = get_sequence_from_fasta("probe4", str(fasta_file))
    assert seq == "AAAA"
    
    # Test non-existent probe
    seq = get_sequence_from_fasta("nonexistent", str(fasta_file))
    assert seq is None


def test_load_sequences_from_fasta(tmp_path):
    """Test loading all sequences from FASTA file."""
    # Create a test FASTA file
    fasta_file = tmp_path / "test.fa"
    with open(fasta_file, 'w') as f:
        f.write(">probe1\nATGCATGC\n")
        f.write(">probe2\nGCGCGCGC\n")
        f.write(">H1_probe3\nTTTTTTTT\n")
    
    sequences = load_sequences_from_fasta(str(fasta_file))
    
    # Check that sequences were loaded
    assert 'probe1' in sequences
    assert 'probe2' in sequences
    assert 'H1_probe3' in sequences
    assert 'probe3' in sequences  # Should also be available without H prefix
    
    assert sequences['probe1'] == "ATGCATGC"
    assert sequences['probe2'] == "GCGCGCGC"
    assert sequences['probe3'] == "TTTTTTTT"


def test_extend_blast_output_with_parameters(tmp_path):
    """Test extending BLAST output with calculated parameters."""
    # Create a test FASTA file
    fasta_file = tmp_path / "test.fa"
    with open(fasta_file, 'w') as f:
        f.write(">probe1\nATGCATGC\n")
        f.write(">probe2\nGCGCGCGC\n")
    
    # Create a test BLAST DataFrame (format: qseqid sseqid evalue sstart send ppos mismatch)
    blast_df = pd.DataFrame({
        0: ['probe1', 'probe2'],
        1: ['seq1', 'seq2'],
        2: [0.001, 0.01],  # evalue
        3: [1, 10],        # sstart
        4: [8, 17],        # send
        5: [100, 90],      # ppos
        6: [0, 1]          # mismatch
    })
    
    # Extend with parameters
    extended_df = extend_blast_output_with_parameters(blast_df, str(fasta_file))
    
    # Check that new columns were added
    required_columns = ['sseq', 'qseq', 'GCcontent', 'Lengthnt', 'hairpin_prob',
                       'dimer_DNA', 'dimer_DNA_flank', 'dimer_probe', 'dimer_probe_DNA',
                       'evalue', 'mismatches', 'length', 'identity', 'Formamide']
    for col in required_columns:
        assert col in extended_df.columns, f"Column {col} not found"
    
    # Check that sequences were loaded
    assert extended_df.loc[0, 'sseq'] == 'ATGCATGC'
    assert extended_df.loc[1, 'sseq'] == 'GCGCGCGC'
    
    # Check GC content calculation
    assert abs(extended_df.loc[0, 'GCcontent'] - 50.0) < 0.1
    assert extended_df.loc[1, 'GCcontent'] == 100.0
    
    # Check length calculation
    assert extended_df.loc[0, 'Lengthnt'] == 8
    assert extended_df.loc[1, 'Lengthnt'] == 8
    
    # Check that mismatches column was created
    assert 'mismatches' in extended_df.columns
    
    # Check that length was calculated from sstart and send
    assert extended_df.loc[0, 'length'] == 8  # send - sstart + 1 = 8 - 1 + 1 = 8
    assert extended_df.loc[1, 'length'] == 8  # 17 - 10 + 1 = 8
