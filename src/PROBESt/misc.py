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


import os
import re
import subprocess
import glob
import pandas as pd
import numpy as np
from typing import List, Dict, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import Namespace


def write_stats(stats: dict, output_dir: str):
    """
    Writes statistics about probe performance to a CSV file.

    Args:
        stats (dict): A dictionary where the keys are iteration numbers (int) and the
                      values are dictionaries containing the following keys:
                      - 'max_hits' (int): The maximum number of hits for the iteration.
                      - 'mean_hits' (float): The mean number of hits for the iteration.
        output_dir (str): The path to the output directory where the CSV file will be saved.

    Returns:
        None: The function writes the statistics to a file and does not return a value.

    Raises:
        OSError: If the output directory cannot be created or the file cannot be written.
    """
    # Create the output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    # Construct the full path to the stats file
    stats_file = os.path.join(output_dir, 'stats.csv')

    try:
        # Write the statistics to the CSV file
        with open(stats_file, 'w') as f:
            # Write the header row
            f.write('iteration,max_hits,mean_hits\n')

            # Write each iteration's statistics
            for iter, stat in stats.items():
                f.write(f"{iter},{stat['max_hits']},{stat['mean_hits']}\n")
    except OSError as e:
        raise OSError(f"Failed to write statistics to {stats_file}: {e}")


def write_fasta(probes, output_file):
    """
    Writes primer sequences to a FASTA file.

    Args:
        probes (list): A list of tuples, where each tuple contains:
                       - sequence_id (str): The ID of the sequence.
                       - primer_num (str): The primer number.
                       - side (str): The type of primer ("LEFT" or "RIGHT").
                       - sequence (str): The primer sequence.
        output_file (str): Path to the output FASTA file.

    Returns:
        None
    """
    try:
        with open(output_file, "w") as fasta:
            for primer in probes:
                sequence_id, primer_num, side, sequence = primer
                header = f">{sequence_id}_{primer_num}_{side}"
                fasta.write(f"{header}\n{sequence}\n")
    except IOError as e:
        print(f"Error writing to file {output_file}: {e}")

def pairing(x):
    """
    Generates a pair of primer names by appending "RIGHT" and "LEFT" to a base name.

    This function takes a primer name (e.g., "primer1_LEFT" or "primer1_RIGHT") and removes
    the "LEFT" or "RIGHT" suffix to extract the base name. It then generates a pair of
    primer names by appending "RIGHT" and "LEFT" to the base name.

    Args:
        x (str): The primer name, which may contain a "LEFT" or "RIGHT" suffix.

    Returns:
        list: A list containing two strings: the base name with "RIGHT" appended and
              the base name with "LEFT" appended.
    """
    # Remove "LEFT" or "RIGHT" from the primer name
    clear_primer = re.sub(r"(RIGHT)|(LEFT)", "", string=x)

    # Return a list with the base name appended by "RIGHT" and "LEFT"
    return [clear_primer + "RIGHT", clear_primer + "LEFT"]


def collect_fasta_files(inputs: List[str]) -> List[str]:
    """
    Collect all FASTA files from input paths (files or directories).
    
    Args:
        inputs: List of input paths. Each can be:
            - A FASTA file (*.fa, *.fna, *.fasta or their .gz versions)
            - A directory containing FASTA files
    
    Returns:
        List of paths to all FASTA files found.
        
    Raises:
        ValueError: If an input path doesn't exist or no FASTA files are found.
    """
    fasta_extensions = ['.fa', '.fna', '.fasta', '.fa.gz', '.fna.gz', '.fasta.gz']
    fasta_files = []
    
    for input_path in inputs:
        if not os.path.exists(input_path):
            raise ValueError(f"Input path does not exist: {input_path}")
        
        if os.path.isfile(input_path):
            # Check if it's a FASTA file
            if any(input_path.endswith(ext) for ext in fasta_extensions):
                fasta_files.append(input_path)
            else:
                raise ValueError(f"Input file is not a recognized FASTA file: {input_path}")
        elif os.path.isdir(input_path):
            # Collect all FASTA files from directory
            dir_files = []
            for ext in fasta_extensions:
                pattern = os.path.join(input_path, f"*{ext}")
                dir_files.extend(glob.glob(pattern))
            
            if not dir_files:
                raise ValueError(f"No FASTA files found in directory: {input_path}")
            
            fasta_files.extend(sorted(dir_files))
        else:
            raise ValueError(f"Input path is neither a file nor a directory: {input_path}")
    
    if not fasta_files:
        raise ValueError("No FASTA files found in any of the provided inputs")
    
    return sorted(fasta_files)


def process_multiple_inputs(
    args: 'Namespace',
    input_fasta_files: List[str],
    output_dir: str,
    initial_generator: str = 'primer3'
) -> str:
    """
    Process multiple input FASTA files separately for initial set generation and merge their outputs.
    
    This function processes each input FASTA file independently:
    1. Creates a temporary directory for each input
    2. Runs uniline_fasta and initial_set_generation for each
    3. Merges individual outputs
    4. Combines all outputs into a single output.fa file
    
    Args:
        args: Arguments object containing all pipeline parameters
        input_fasta_files: List of paths to FASTA files to process
        output_dir: Output directory path (e.g., .tmp/0/)
        initial_generator: Initial generator to use ('primer3' or 'oligominer')
    
    Returns:
        Path to the final merged output.fa file
    
    Raises:
        ValueError: If initial_generator is unknown
        ImportError: If required modules cannot be imported
    """
    # Import here to avoid circular dependencies
    from PROBESt.bash_wrappers import uniline_fasta
    from PROBESt.merge import merge
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each input FASTA separately for initial set generation
    temp_output_dirs = []
    
    for idx, fasta_file in enumerate(input_fasta_files):
        print(f"\nProcessing input {idx + 1}/{len(input_fasta_files)}: {fasta_file}")
        
        # Create temporary directory for this input
        temp_dir = os.path.join(output_dir, f"input_{idx}/")
        os.makedirs(temp_dir, exist_ok=True)
        temp_output_dirs.append(temp_dir)
        
        # Create a temporary args object with this specific input file
        temp_args = type('Args', (), {})()
        for attr in dir(args):
            if not attr.startswith('_'):
                setattr(temp_args, attr, getattr(args, attr))
        temp_args.input = fasta_file
        
        # Make uniline fasta for this input
        uniline_fasta(temp_args, temp_dir)
        print(f"Input fasta parsed: {fasta_file}")
        
        # Template generation for this input
        if initial_generator == "primer3":
            from PROBESt.primer3 import initial_set_generation
            initial_set_generation(temp_args, temp_dir)
        elif initial_generator == "oligominer":
            from PROBESt.oligominer import initial_set_generation
            initial_set_generation(temp_args, temp_dir)
        else:
            raise ValueError(f"Unknown initial generator: {initial_generator}")
        
        # Merge this input's output
        merge(algo=temp_args.algorithm,
              input=os.path.join(temp_dir, "output.fa"),
              output=os.path.join(temp_dir, "merged.fa"),
              tmp=os.path.join(temp_dir, "fasta_table.tsv"),
              NNN=10,
              script_path=temp_args.script_path)
    
    # Merge all outputs to the first directory's output.fa
    print("\nMerging all initial set outputs...")
    final_output_fa = os.path.join(output_dir, "output.fa")
    with open(final_output_fa, 'w') as outfile:
        for temp_dir in temp_output_dirs:
            merged_fa = os.path.join(temp_dir, "merged.fa")
            if os.path.exists(merged_fa) and os.path.getsize(merged_fa) > 0:
                with open(merged_fa, 'r') as infile:
                    outfile.write(infile.read())
            else:
                # If merged.fa doesn't exist, try output.fa
                output_fa = os.path.join(temp_dir, "output.fa")
                if os.path.exists(output_fa) and os.path.getsize(output_fa) > 0:
                    with open(output_fa, 'r') as infile:
                        outfile.write(infile.read())
    
    return final_output_fa


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate GC content percentage of a DNA/RNA sequence.
    
    Args:
        sequence (str): DNA/RNA sequence
        
    Returns:
        float: GC content as percentage (0-100)
    """
    if not sequence or len(sequence) == 0:
        return 0.0
    
    seq_upper = sequence.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    total_valid = sum(seq_upper.count(base) for base in 'ATGCU')
    
    if total_valid == 0:
        return 0.0
    
    return (gc_count / total_valid) * 100


def get_sequence_from_fasta(probe_name: str, fasta_file: str) -> Optional[str]:
    """
    Extract sequence from FASTA file by probe name.
    
    Args:
        probe_name (str): Name of the probe (header without '>')
        fasta_file (str): Path to FASTA file
        
    Returns:
        Optional[str]: Sequence if found, None otherwise
    """
    try:
        with open(fasta_file, 'r') as f:
            current_name = None
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Extract name, handling potential prefixes like "H1_"
                    current_name = line[1:].strip()
                    # Remove hit count prefix if present (e.g., "H1_probe_name" -> "probe_name")
                    if re.match(r'^H\d+_', current_name):
                        current_name = re.sub(r'^H\d+_', '', current_name)
                    # Remove LEFT/RIGHT suffix for primer algorithm
                    if '_LEFT' in current_name or '_RIGHT' in current_name:
                        current_name = re.sub(r'(_LEFT|_RIGHT)$', '', current_name)
                else:
                    if current_name:
                        # Check if this is the probe we're looking for
                        name_to_check = current_name
                        if re.match(r'^H\d+_', name_to_check):
                            name_to_check = re.sub(r'^H\d+_', '', name_to_check)
                        if '_LEFT' in name_to_check or '_RIGHT' in name_to_check:
                            name_to_check = re.sub(r'(_LEFT|_RIGHT)$', '', name_to_check)
                        
                        if name_to_check == probe_name or current_name == probe_name:
                            return line.upper()
        return None
    except Exception as e:
        print(f"Error reading FASTA file {fasta_file}: {e}")
        return None


def load_sequences_from_fasta(fasta_file: str) -> Dict[str, str]:
    """
    Load all sequences from a FASTA file into a dictionary.
    
    Args:
        fasta_file (str): Path to FASTA file
        
    Returns:
        Dict[str, str]: Dictionary mapping probe names to sequences
    """
    sequences = {}
    try:
        with open(fasta_file, 'r') as f:
            current_name = None
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    current_name = line[1:].strip()
                else:
                    if current_name:
                        # Store with original name
                        sequences[current_name] = line.upper()
                        # Also store without H prefix and LEFT/RIGHT suffix for easier lookup
                        name_clean = current_name
                        if re.match(r'^H\d+_', name_clean):
                            name_clean = re.sub(r'^H\d+_', '', name_clean)
                        if '_LEFT' in name_clean or '_RIGHT' in name_clean:
                            name_clean = re.sub(r'(_LEFT|_RIGHT)$', '', name_clean)
                        if name_clean != current_name:
                            sequences[name_clean] = line.upper()
        return sequences
    except Exception as e:
        print(f"Error reading FASTA file {fasta_file}: {e}")
        return {}


def extend_blast_output_with_parameters(
    blast_df: pd.DataFrame,
    fasta_file: str,
    model_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Extend BLAST output DataFrame with calculated parameters needed for AI model.
    
    BLAST format expected: qseqid, sseqid, evalue, sstart, send, ppos, mismatch
    
    Args:
        blast_df (pd.DataFrame): BLAST output DataFrame
        fasta_file (str): Path to FASTA file containing probe sequences
        model_path (Optional[str]): Path to saved AI model (if None, only calculates parameters)
        
    Returns:
        pd.DataFrame: Extended DataFrame with all calculated parameters
    """
    # Load sequences from FASTA
    sequences = load_sequences_from_fasta(fasta_file)
    
    # Rename columns if needed (BLAST output format: qseqid sseqid evalue sstart send ppos mismatch)
    if len(blast_df.columns) >= 7:
        blast_df.columns = ['qseqid', 'sseqid', 'evalue', 'sstart', 'send', 'ppos', 'mismatch'] + list(blast_df.columns[7:])
    
    # Create mismatches column to match training data format (from mismatch column)
    if 'mismatch' in blast_df.columns and 'mismatches' not in blast_df.columns:
        blast_df['mismatches'] = blast_df['mismatch']
    
    # Initialize new columns
    blast_df['sseq'] = None
    blast_df['qseq'] = None
    blast_df['left_flank'] = ''
    blast_df['right_flank'] = ''
    blast_df['Formamide'] = 0.0
    blast_df['GCcontent'] = 0.0
    blast_df['Lengthnt'] = 0
    blast_df['Modifiedversions'] = ''
    blast_df['length'] = blast_df['send'] - blast_df['sstart'] + 1
    blast_df['bitscore'] = 0.0  # Will be calculated if needed
    blast_df['identity'] = 0.0
    blast_df['qseq_aln'] = None
    blast_df['sseq_aln'] = None
    blast_df['score'] = 0.0
    blast_df['hairpin_prob'] = 0.0
    blast_df['dimer_DNA'] = 0.0
    blast_df['dimer_DNA_flank'] = 0.0
    blast_df['dimer_probe'] = 0.0
    blast_df['dimer_probe_DNA'] = 0.0
    
    # Import RNA structure functions
    try:
        from PROBESt.rna_structure import calculate_hairpin_prob, calculate_dimer_G
    except ImportError:
        print("Warning: Could not import RNA structure functions. Some parameters will be set to 0.")
        calculate_hairpin_prob = lambda x: 0.0
        calculate_dimer_G = lambda x, y=None, type1="DNA", type2=None: 0.0
    
    # Process each row
    for idx, row in blast_df.iterrows():
        probe_name = str(row['qseqid'])
        
        # Get sequence from FASTA
        seq = sequences.get(probe_name)
        if seq is None:
            # Try without H prefix
            if re.match(r'^H\d+_', probe_name):
                probe_name_clean = re.sub(r'^H\d+_', '', probe_name)
                seq = sequences.get(probe_name_clean)
            if seq is None:
                # Try with LEFT/RIGHT removal
                probe_name_clean = re.sub(r'(_LEFT|_RIGHT)$', '', probe_name)
                seq = sequences.get(probe_name_clean)
        
        if seq:
            blast_df.at[idx, 'sseq'] = seq
            blast_df.at[idx, 'qseq'] = seq  # For now, use same sequence
            blast_df.at[idx, 'GCcontent'] = calculate_gc_content(seq)
            blast_df.at[idx, 'Lengthnt'] = len(seq)
            
            # Calculate identity from mismatch and length
            mismatch_val = row.get('mismatch', row.get('mismatches', 0))
            if pd.notna(mismatch_val) and pd.notna(row['length']):
                mismatches = int(mismatch_val)
                length = int(row['length'])
                if length > 0:
                    identity = ((length - mismatches) / length) * 100
                    blast_df.at[idx, 'identity'] = identity
            
            # Calculate hairpin probability
            try:
                blast_df.at[idx, 'hairpin_prob'] = calculate_hairpin_prob(seq)
            except:
                blast_df.at[idx, 'hairpin_prob'] = 0.0
            
            # Calculate dimer energies
            try:
                blast_df.at[idx, 'dimer_DNA'] = calculate_dimer_G(seq, type1="DNA", type2="DNA")
            except:
                blast_df.at[idx, 'dimer_DNA'] = 0.0
            
            try:
                # For dimer_DNA_flank, we'd need flanking sequences, but we don't have them
                # So we'll use the same as dimer_DNA for now
                blast_df.at[idx, 'dimer_DNA_flank'] = blast_df.at[idx, 'dimer_DNA']
            except:
                blast_df.at[idx, 'dimer_DNA_flank'] = 0.0
            
            try:
                blast_df.at[idx, 'dimer_probe'] = calculate_dimer_G(seq, type1="RNA", type2="RNA")
            except:
                blast_df.at[idx, 'dimer_probe'] = 0.0
            
            try:
                blast_df.at[idx, 'dimer_probe_DNA'] = calculate_dimer_G(seq, seq, type1="RNA", type2="DNA")
            except:
                blast_df.at[idx, 'dimer_probe_DNA'] = 0.0
            
            # Set alignment sequences (same as original for now)
            blast_df.at[idx, 'qseq_aln'] = seq
            blast_df.at[idx, 'sseq_aln'] = seq
            
            # Calculate score (simplified - could be improved)
            blast_df.at[idx, 'score'] = float(row['evalue']) if pd.notna(row['evalue']) else 0.0
        else:
            # If sequence not found, set defaults
            blast_df.at[idx, 'sseq'] = ''
            blast_df.at[idx, 'qseq'] = ''
            blast_df.at[idx, 'GCcontent'] = 0.0
            blast_df.at[idx, 'Lengthnt'] = 0
    
    # Ensure numeric columns are numeric
    numeric_cols = ['evalue', 'mismatch', 'mismatches', 'length', 'identity', 'GCcontent', 'Lengthnt', 
                    'score', 'hairpin_prob', 'dimer_DNA', 'dimer_DNA_flank', 
                    'dimer_probe', 'dimer_probe_DNA', 'Formamide']
    for col in numeric_cols:
        if col in blast_df.columns:
            blast_df[col] = pd.to_numeric(blast_df[col], errors='coerce').fillna(0.0)
    
    return blast_df
