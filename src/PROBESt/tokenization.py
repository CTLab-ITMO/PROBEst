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


"""Tokenization module for DNA sequences.

This module provides functions to tokenize DNA sequences into k-mers,
similar to how DNA language models process sequences.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from collections import Counter


def tokenize_seq(sequence: str, k: int = 3) -> List[str]:
    """Tokenize a DNA sequence into k-mers.
    
    Args:
        sequence: DNA sequence string (e.g., "ATCGATCG")
        k: Size of k-mers (default: 3, producing 3-mers like "ATC", "TCG")
    
    Returns:
        List of k-mer tokens extracted from the sequence.
        For example, "ATCGATCG" with k=3 returns ["ATC", "TCG", "CGA", "GAT", "ATC", "TCG"]
    
    Example:
        >>> tokenize_seq("ATCGATCG", k=3)
        ['ATC', 'TCG', 'CGA', 'GAT', 'ATC', 'TCG']
    """
    if not sequence or pd.isna(sequence):
        return []
    
    sequence = str(sequence).upper().strip()
    if len(sequence) < k:
        return []
    
    tokens = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        # Only include valid DNA k-mers (containing only A, T, C, G)
        if all(base in 'ATCG' for base in kmer):
            tokens.append(kmer)
    
    return tokens


def tokenize_table(input_csv: str, output_csv: Optional[str] = None, 
                  add_tokens: int = 100, k: int = 3, 
                  sequence_columns: Optional[List[str]] = None,
                  drop_original_sequences: bool = True) -> pd.DataFrame:
    """Tokenize DNA sequences in a CSV table and add token count columns.
    
    This function reads a CSV file, tokenizes sequences in specified columns
    (default: 'sseq' and 'qseq'), and adds new columns with k-mer counts.
    The new columns are named like 'sseq_token_TTC', 'qseq_token_AAA', etc.
    
    Args:
        input_csv: Path to input CSV file
        output_csv: Path to output CSV file (if None, appends '_tokenized' to input name)
        add_tokens: Number of top k-mers to add as columns per sequence column (default: 100)
        k: Size of k-mers (default: 3)
        sequence_columns: List of column names to tokenize (default: ['sseq', 'qseq'])
        drop_original_sequences: If True, drop original sequence columns after tokenization (default: True)
    
    Returns:
        DataFrame with original columns plus new token count columns (or without original sequences if dropped)
    
    Example:
        >>> df = tokenize_table('data.csv', add_tokens=50)
        >>> # Adds columns like 'sseq_token_AAA', 'sseq_token_TTC', etc.
    """
    if sequence_columns is None:
        sequence_columns = ['sseq', 'qseq']
    
    # Read the CSV
    df = pd.read_csv(input_csv)
    
    # Collect all k-mers from all sequences to find the most frequent ones
    all_kmers = Counter()
    
    for col in sequence_columns:
        if col not in df.columns:
            print(f"Warning: Column '{col}' not found in CSV. Skipping.")
            continue
        
        for seq in df[col]:
            if pd.notna(seq):
                tokens = tokenize_seq(str(seq), k=k)
                all_kmers.update(tokens)
    
    # Get top k-mers to add as columns
    top_kmers = [kmer for kmer, _ in all_kmers.most_common(add_tokens)]
    
    print(f"Found {len(all_kmers)} unique {k}-mers. Adding top {len(top_kmers)} as columns.")
    
    # Build all token count columns at once to avoid DataFrame fragmentation
    new_columns = {}
    
    for col in sequence_columns:
        if col not in df.columns:
            continue
        
        # Pre-compute tokens for all sequences to avoid repeated computation
        print(f"Tokenizing {col} column...")
        all_tokens = df[col].apply(
            lambda seq: tokenize_seq(str(seq), k=k) if pd.notna(seq) else []
        )
        
        # Count k-mers for each row and each top k-mer
        for kmer in top_kmers:
            col_name = f"{col}_token_{kmer}"
            new_columns[col_name] = all_tokens.apply(
                lambda tokens: tokens.count(kmer) if isinstance(tokens, list) else 0
            )
    
    # Add all new columns at once using pd.concat to avoid fragmentation
    if new_columns:
        new_df = pd.DataFrame(new_columns, index=df.index)
        df = pd.concat([df, new_df], axis=1)
    
    # Count added columns before dropping original sequences
    num_added_columns = len(new_columns) if new_columns else 0
    
    # Optionally drop original sequence columns after tokenization
    if drop_original_sequences:
        for col in sequence_columns:
            if col in df.columns:
                df = df.drop(columns=[col])
                print(f"Dropped original sequence column: {col}")
    
    # Save to output file
    if output_csv is None:
        base_name = input_csv.rsplit('.', 1)[0]
        output_csv = f"{base_name}_tokenized.csv"
    
    df.to_csv(output_csv, index=False)
    print(f"Tokenized table saved to {output_csv}")
    print(f"Added {num_added_columns} new token columns.")
    
    return df

