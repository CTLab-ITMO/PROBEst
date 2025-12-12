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


"""Module for analyzing BLAST alignments of probe-target sequence pairs.

This module provides functionality to process BLAST alignment results and calculate
various sequence properties for probe-target pairs. It handles both standard DNA
sequences and sequences containing IUPAC ambiguous nucleotide codes.
"""

import pandas as pd
from Bio import SearchIO
from typing import Dict, Any


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of a DNA sequence, handling IUPAC ambiguous bases.

    IUPAC nucleotide codes that contribute to GC content:
    - G, C: 100% GC
    - S: 100% GC (G or C)
    - K, M: 50% GC (K = G/T, M = A/C)
    - B, V: 66.7% GC (B = C/G/T, V = A/C/G)
    - D, H: 33.3% GC (D = A/G/T, H = A/C/T)
    - N: 50% GC (any base)
    - R, Y: 50% GC (R = A/G, Y = C/T)
    """
    sequence = sequence.upper()
    length = len(sequence)
    if not length:
        return 0

    # Full GC bases
    gc_count = sequence.count('G') + sequence.count('C') + sequence.count('S')

    # 50% GC bases
    gc_count += 0.5 * (
        sequence.count('K') + sequence.count('M') +
        sequence.count('N') + sequence.count('R') +
        sequence.count('Y'))

    # 66.7% GC bases
    gc_count += (2/3) * (sequence.count('B') + sequence.count('V'))

    # 33.3% GC bases
    gc_count += (1/3) * (sequence.count('D') + sequence.count('H'))

    return (gc_count / length) * 100


def calculate_melting_temperature(sequence: str) -> float:
    """Calculate melting temperature of a DNA sequence using Biopython's MeltingTemp.

    This uses Biopython's implementation of nearest-neighbor thermodynamics
    which is more accurate than simple rules like Wallace.
    """
    from Bio.SeqUtils import MeltingTemp as mt
    sequence = sequence.upper()
    if not sequence:
        return 0
    return round(mt.Tm_NN(sequence, strict=False), 2)


def process_sequence_pair(probe_seq: str, target_seq: str) -> Dict[str, Any]:
    """Process a pair of sequences (probe and target) and calculate their features.

    This function calculates various properties for both probe and target sequences.
    New feature calculations can be added here.
    """
    features = {
        'probe_sequence': probe_seq,
        'target_sequence': target_seq,

        # Paired properties (calculated for both probe and target)
        'probe_gc_content': calculate_gc_content(probe_seq),
        'target_gc_content': calculate_gc_content(target_seq),
        'probe_tm': calculate_melting_temperature(probe_seq),
        'target_tm': calculate_melting_temperature(target_seq),
        'probe_length': len(probe_seq),
        'target_length': len(target_seq),
        # Single properties

    }
    return features


def profile_alignment(blast_file: str) -> pd.DataFrame:
    """Process BLAST alignment file and extract sequence features.

    Args:
        blast_file: Path to the BLAST alignment file in tabular format

    Returns:
        DataFrame containing sequence pairs and their calculated features
    """

    # Code requires a BLAST alignment file with the following fields:
    # qseqid, sseqid, evalue, sstart, send, ppos, mismatch, qseq, sseq

    custom_fields = ['qseqid', 'sseqid', 'evalue', 'sstart', 'send', 'ppos', 'mismatch', 'qseq', 'sseq']
    blast_records = SearchIO.parse(blast_file, 'blast-tab', fields=custom_fields)
    results = []

    for qresult in blast_records:
        for hit in qresult:
            # Multimappings should be filtered before this step
            # Optionally we can add this filter here too
            for fragment in hit[0].fragments:
                # Basic alignment information
                result = {
                    'query_id': qresult.id,
                    'target_id': hit.id,
                }
                # Calculate sequence-based features
                probe_seq = str(fragment.query.seq)
                target_seq = str(fragment.hit.seq)
                features = process_sequence_pair(probe_seq, target_seq)
                result.update(features)
                results.append(result)

    return pd.DataFrame(results)
