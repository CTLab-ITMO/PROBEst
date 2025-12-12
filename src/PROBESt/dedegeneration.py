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


"""
De-degeneration module for PROBESt.

This module provides functionality to reduce degenerate nucleotides in probe sequences
while maintaining their specificity. It iteratively tries to replace degenerate nucleotides
with specific ones (ATGCU) and validates the results using BLAST and probe checking.
"""

import os
import subprocess
import numpy as np
import pandas as pd
import re
from typing import Dict, Tuple

# Try to import dedegenerate_position from evolution module
# If not available (package not installed), define it here
try:
    from PROBESt.evolution import dedegenerate_position
except ImportError:
    # Define locally if import fails
    import random
    def dedegenerate_position(x, dedegen_rate):
        """De-degenerates a nucleotide based on a given de-degeneration rate."""
        degenerate_map = {
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
        if x.upper() not in degenerate_map:
            return x
        rstate = random.random()
        if rstate <= dedegen_rate:
            return random.choice(degenerate_map[x.upper()])
        return x

# Define degenerate nucleotide codes
DEGENERATE_CODES = "WSMKRYBDHVN"


def count_degenerate_nucleotides(sequence: str) -> int:
    """
    Count the number of degenerate nucleotides in a sequence.
    
    Args:
        sequence: DNA/RNA sequence string
        
    Returns:
        Number of degenerate nucleotides in the sequence
    """
    return sum(1 for char in sequence.upper() if char in DEGENERATE_CODES)


def read_fasta_sequences(fasta_path: str, algorithm: str = 'FISH') -> Dict[str, str]:
    """
    Read sequences from a FASTA file and return as a dictionary.
    
    Args:
        fasta_path: Path to the FASTA file
        algorithm: Algorithm type ('FISH' or 'primer')
        
    Returns:
        Dictionary mapping sequence names to sequences
    """
    seqs = {}
    fasta = open(fasta_path, "r")
    for iter_line, line in enumerate(fasta):
        if iter_line % 2 == 0:
            line_name = line[1:-1]  # Remove '>' and newline
            if algorithm == 'primer':
                line_clear = re.sub(r'(_LEFT|_RIGHT)$', '', line_name)
            else:
                line_clear = line_name
            keep = True
        else:
            if keep:
                seqs[line_name] = line[:-1]
    fasta.close()
    return seqs


def run_dedegeneration_iteration(
    args,
    input_fasta: str,
    output_dir: str,
    iteration: int,
    blastn_cmd: str,
    probe_check_cmd: str
) -> Tuple[Dict[str, str], pd.Series]:
    """
    Run one iteration of the de-degeneration algorithm.
    
    Args:
        args: Arguments object with all parameters
        input_fasta: Path to input FASTA file
        output_dir: Output directory for this iteration
        iteration: Current iteration number
        blastn_cmd: Pre-built blastn command string
        probe_check_cmd: Pre-built probe_check command string
        
    Returns:
        Tuple of (sequences dictionary, probe_vals Series)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Read input sequences
    seqs = read_fasta_sequences(input_fasta, args.algorithm)
    
    # Generate de-degenerated variants
    if args.dedegeneration_append:
        seqs_dedegen = seqs.copy()
    else:
        seqs_dedegen = dict()
    
    dedegen_rate = getattr(args, 'dedegeneration_rate', 0.5)
    
    for seq_name, init_seq in seqs.items():
        for seq_iter in range(args.dedegeneration_set_size):
            dedegenerated_seq = init_seq
            # Try to de-degenerate at least one position
            attempts = 0
            while dedegenerated_seq == init_seq and attempts < 100:
                dedegenerated_seq = "".join([
                    dedegenerate_position(_, dedegen_rate=dedegen_rate) 
                    for _ in init_seq
                ])
                attempts += 1
            
            # Only add if we actually changed something
            if dedegenerated_seq != init_seq:
                dseq_name = "DD" + str(iteration) + "N" + str(seq_iter) + "_" + seq_name
                seqs_dedegen[dseq_name] = dedegenerated_seq
    
    # Write sequences to FASTA
    fasta_path = os.path.join(output_dir, "output.fa")
    fasta = open(fasta_path, "w")
    for fname in seqs_dedegen.keys():
        fasta.write(">" + fname + "\n" + seqs_dedegen[fname] + "\n")
    fasta.close()
    
    # Merge sequences
    from PROBESt.merge import merge
    merge(
        algo=args.algorithm,
        input=fasta_path,
        output=os.path.join(output_dir, "merged.fa"),
        tmp=os.path.join(output_dir, "fasta_table.tsv"),
        NNN=10,
        script_path=args.script_path
    )
    
    # Run BLASTn against true base
    blastn_iter = blastn_cmd + " -query " + os.path.join(output_dir, "merged.fa")
    blastn_db = blastn_iter + " -db " + args.true_base + \
        " > " + os.path.join(output_dir, "positive_hits.tsv")
    subprocess.run(blastn_db, shell=True)
    
    # Run BLASTn against false bases
    for db_neg in args.false_base:
        blastn_db = blastn_iter + " -db " + db_neg + \
            " >> " + os.path.join(output_dir, "negative_hits.tsv")
        subprocess.run(blastn_db, shell=True)
    
    # Run probe check
    probe_check_iter = probe_check_cmd + \
        " -o " + os.path.join(output_dir, "clear_hits.tsv") + \
        " -r " + os.path.join(output_dir, "probe_check/") + \
        " -t " + os.path.join(output_dir, "positive_hits.tsv") + " " + \
        os.path.join(output_dir, "negative_hits.tsv")
    
    subprocess.run(probe_check_iter, shell=True)
    
    # Read probe check results
    try:
        probe_out = pd.read_table(os.path.join(output_dir, "clear_hits.tsv"),
                                  sep=' ', header=None)
    except:
        raise InterruptedError(
            "Empty file after filtration in de-degeneration iteration, "
            "try to use other probe_check properties and review false databases")
    
    probe_vals = probe_out.iloc[:, 0].value_counts()
    
    return seqs_dedegen, probe_vals


def run_dedegeneration(args, input_fasta: str, output_fasta: str) -> str:
    """
    Run the full de-degeneration algorithm.
    
    This function iteratively tries to reduce degenerate nucleotides in probe sequences
    while maintaining their specificity through BLAST validation.
    
    Args:
        args: Arguments object with all parameters
        input_fasta: Path to input FASTA file (output from evolutionary algorithm)
        output_fasta: Path to output FASTA file
        
    Returns:
        Path to the output FASTA file
    """
    print("\n---- De-degeneration module ----")
    
    # Prepare commands
    from PROBESt.bash_wrappers import blastn_function, probe_check_function
    blastn_cmd = blastn_function(args)
    probe_check_cmd = probe_check_function(args)
    
    # Set up output directory structure
    dedegen_tmp_dir = os.path.join(args.output_tmp, "dedegeneration")
    os.makedirs(dedegen_tmp_dir, exist_ok=True)
    
    # Start with input sequences
    current_fasta = input_fasta
    all_stats = {}
    seqs_selected = {}
    probe_vals = pd.Series(dtype=int)
    
    # Run de-degeneration iterations
    for dedeg_iter in range(1, args.dedegeneration_iterations + 1):
        print(f"\nDe-degeneration iteration {dedeg_iter} ----")
        
        iter_dir = os.path.join(dedegen_tmp_dir, str(dedeg_iter))
        seqs_dedegen, probe_vals = run_dedegeneration_iteration(
            args, current_fasta, iter_dir, dedeg_iter,
            blastn_cmd, probe_check_cmd
        )
        
        # Calculate stats
        if len(probe_vals) > 0:
            max_hits = probe_vals.iloc[0]
            mean_hits = round(sum(probe_vals) / len(probe_vals), 1)
        else:
            max_hits = 0
            mean_hits = 0
        
        all_stats[dedeg_iter] = {
            'max_hits': max_hits,
            'mean_hits': mean_hits
        }
        print(f"Maximum hits: {max_hits}")
        print(f"Mean hits: {mean_hits}")
        
        # Select best sequences based on:
        # 1. They passed probe_check (are in probe_vals)
        # 2. Fewer degenerate nucleotides is better
        probe_list = list(set(probe_vals.index))
        probe_list_hash = [hash(_) for _ in probe_list]
        
        # Read sequences and filter
        seqs_filtered = {}
        fasta = open(os.path.join(iter_dir, "output.fa"), "r")
        for iter_line, line in enumerate(fasta):
            if iter_line % 2 == 0:
                line_name = line[1:-1]
                if args.algorithm == 'primer':
                    line_clear = re.sub(r'(_LEFT|_RIGHT)$', '', line_name)
                else:
                    line_clear = line_name
                
                if np.isin(hash(line_clear), probe_list_hash):
                    keep = True
                else:
                    keep = False
            else:
                if keep:
                    seqs_filtered[line_name] = line[:-1]
        fasta.close()
        
        # Sort sequences by number of degenerate nucleotides (ascending)
        # and by number of hits (descending)
        seqs_with_scores = []
        for seq_name, seq in seqs_filtered.items():
            deg_count = count_degenerate_nucleotides(seq)
            if args.algorithm == 'primer':
                seq_key = re.sub(r'(_LEFT|_RIGHT)$', '', seq_name)
            else:
                seq_key = seq_name
            
            if seq_key in probe_vals.index:
                hits = probe_vals.loc[seq_key]
            else:
                hits = 0
            
            seqs_with_scores.append((seq_name, seq, deg_count, hits))
        
        # Sort: first by degenerate count (ascending), then by hits (descending)
        seqs_with_scores.sort(key=lambda x: (x[2], -x[3]))
        
        # Select top sequences (those with fewest degenerate nucleotides)
        # Keep sequences that passed filtering
        seqs_selected = {}
        for seq_name, seq, deg_count, hits in seqs_with_scores:
            if args.algorithm == 'primer':
                seq_key = re.sub(r'(_LEFT|_RIGHT)$', '', seq_name)
            else:
                seq_key = seq_name
            
            # Only keep if it passed probe check and we don't have too many already
            if seq_key in probe_vals.index:
                # Keep if we haven't seen this sequence key before
                if seq_key not in [re.sub(r'(_LEFT|_RIGHT)$', '', k) if args.algorithm == 'primer' 
                                  else k for k in seqs_selected.keys()]:
                    seqs_selected[seq_name] = seq
        
        # Update current_fasta for next iteration
        if dedeg_iter < args.dedegeneration_iterations:
            current_fasta = os.path.join(iter_dir, "selected.fa")
            fasta = open(current_fasta, "w")
            for fname in seqs_selected.keys():
                fasta.write(">" + fname + "\n" + seqs_selected[fname] + "\n")
            fasta.close()
        
        print(f"Selected {len(seqs_selected)} sequences")
        print("Done")
    
    # If no sequences were selected, use original input
    if len(seqs_selected) == 0:
        print("Warning: No sequences passed de-degeneration filtering. Using original sequences.")
        seqs_selected = read_fasta_sequences(input_fasta, args.algorithm)
        # Create empty probe_vals if needed
        if len(probe_vals) == 0:
            probe_vals = pd.Series(index=list(seqs_selected.keys()), data=0)
    
    # Write final output
    # Use sequences from the last iteration
    fasta = open(output_fasta, "w")
    for fname, seq in seqs_selected.items():
        if args.algorithm == 'primer':
            seq_key = re.sub(r'(_LEFT|_RIGHT)$', '', fname)
        else:
            seq_key = fname
        
        if seq_key in probe_vals.index:
            seqs_hits = str(probe_vals.loc[seq_key])
        else:
            seqs_hits = "0"
        
        deg_count = count_degenerate_nucleotides(seq)
        fasta.write(">H" + seqs_hits + "_D" + str(deg_count) + "_" + fname + "\n" + seq + "\n")
    fasta.close()
    
    print(f"\nDe-degeneration completed. Output written to {output_fasta}")
    print(f"Final sequences: {len(seqs_selected)}")
    
    return output_fasta

