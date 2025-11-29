#!/usr/bin/env python3

import pandas as pd
import numpy as np
import random
import re
from typing import Tuple, List, Dict
import argparse

def clean_sequence(sequence: str) -> str:
    """Clean sequence by removing non-ACGT characters and converting to uppercase."""
    sequence = sequence.strip().upper()
    sequence = re.sub(r'[^ACGT]', '', sequence)
    return sequence

def calculate_gc_content(sequence: str) -> int:
    """Calculate GC content percentage of a sequence and round to integer."""
    if not sequence:
        return 0
    gc_count = sequence.count('G') + sequence.count('C')
    return round((gc_count / len(sequence)) * 100)

def apply_mutations_sequence(sequence: str, 
                   insertion_rate: float = 0.01,
                   deletion_rate: float = 0.01,
                   mutation_rate: float = 0.1,
                   values: list=[None]) -> str:
    """Apply random mutations to the sequence."""
    nucleotides = ['A', 'C', 'G', 'T']
    result = []
    
    while sequence == ''.join(result):
        for base in sequence:
            # Apply mutations based on probabilities
            if random.random() < mutation_rate:
                # SNP mutation
                new_base = random.choice([n for n in nucleotides if n != base])
                result.append(new_base)
            else:
                result.append(base)
                
            # Insertion after current position
            if random.random() < insertion_rate:
                result.append(random.choice(nucleotides))
                
            # Deletion of current position
            if random.random() < deletion_rate:
                if result:
                    result.pop()
    
    return ''.join(result)

def apply_mutations_experiment(sequence,
                               values: list=[None],
                               insertion_rate: float = 0.01,
                               deletion_rate: float = 0.01,
                               mutation_rate: float = 0.1):
    """Apply random mutations to the experiment parameters."""
    if random.random() < mutation_rate:
        return random.choice(values)
    else:
        return None

def process_probe_data(df: pd.DataFrame) -> List[Dict]:
    """Process the input DataFrame into a list of probe dictionaries."""
    probes = []
    current_probe = {}
    
    for _, row in df.iterrows():
        if pd.isna(row[0]) or row[0] == '':
            continue
            
        if row[0] == 'Test':
            if current_probe:
                probes.append(current_probe)
            current_probe = {}
        else:
            current_probe[row[0]] = row[1]
    
    if current_probe:
        probes.append(current_probe)
    
    return probes

def generate_noisy_probes(input_file: str,
                         output_file: str,
                         mutation_number: int = 1,
                         obligate_mutations: list = ["Sequence"],
                         facultative_mutations: list = ["Formamide [%]","Modified version(s)"],
                         insertion_rate: float = 0.01,
                         deletion_rate: float = 0.01,
                         mutation_rate: float = 0.1,
                         iterations: int = 1) -> None:
    """Generate noisy probe data with sequence mutations and experiment parameter variations.
    Parameters:
        input_file (str): Path to the input CSV file containing probe data.
        output_file (str): Path to the output file where noisy probes will be saved.
        insertion_rate (float, optional): Probability of inserting a nucleotide at each position. Default is 0.01.
        deletion_rate (float, optional): Probability of deleting a nucleotide at each position. Default is 0.01.
        mutation_rate (float, optional): Probability of introducing a SNP at each position. Default is 0.1.
        iterations (int, optional): Number of noisy probe sets to generate. Default is 1.
    """
       
    # Read input probe data
    df = pd.read_csv(input_file, header=None)
    # Prepare
    mutation_list = [*obligate_mutations, *facultative_mutations]
    mutation_functions = {
        "Sequence":apply_mutations_sequence,
        "Formamide [%]":apply_mutations_experiment,
        "Modified version(s)":apply_mutations_experiment
    }
    mutation_values = {
        "Sequence":None,
        "Formamide [%]":df.loc[df[0] == "Formamide [%]", 1].tolist(),
        "Modified version(s)":df.loc[df[0] == "Modified version(s)", 1].tolist()
    }
    
    # Process probes
    probes = process_probe_data(df)
    all_noisy_probes = []
    
    for iteration in range(iterations):
        noisy_probes = []
        
        for probe in probes:
            if 'Sequence' not in probe:
                continue
                
            # Clean sequence
            clean_seq = clean_sequence(probe['Sequence'])
            probe['Sequence'] = clean_seq
            
            # Apply mutations
            for OM in obligate_mutations:
                tmp = mutation_functions[OM]
                probe[OM] = tmp(sequence = probe[OM],
                                insertion_rate = insertion_rate,
                                deletion_rate = deletion_rate,
                                mutation_rate = mutation_rate,
                                values = mutation_values[OM])
                
            if mutation_number < len(mutation_list):
                OMs = random.choices(mutation_list, k=mutation_number)
            else:
                add = ["Sequence"]*(mutation_number - len(mutation_list))
                OMs = [*mutation_list, *add]
            
            for OM in OMs:
                tmp = mutation_functions[OM]
                probe[OM] = tmp(sequence = probe[OM],
                                insertion_rate = insertion_rate,
                                deletion_rate = deletion_rate,
                                mutation_rate = mutation_rate,
                                values = mutation_values[OM])
                
            # Calculate new properties
            new_length = len(probe['Sequence'])
            new_gc_content = calculate_gc_content(probe['Sequence'])
            probe['Length [nt]'] = str(new_length)
            probe['G+C content [%]'] = str(new_gc_content)
        
            # Generate unique ID for this iteration
            if 'Accession no.' in probe:
                probe['Accession no.'] = f"{probe['Accession no.']}_{iteration + 1}"
                
            probe['Mutation number'] = mutation_number
            noisy_probes.append(probe)
        
        all_noisy_probes.extend(noisy_probes)
    
    # Convert to DataFrame format
    rows = []
    for probe in all_noisy_probes:
        for key, value in probe.items():
            rows.append([key, value, len(rows) // len(probe) + 1])
    
    # Create new DataFrame and save to CSV
    noisy_df = pd.DataFrame(rows, columns=['', '1', 'id'])
    noisy_df.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description='Generate noisy probe data')
    parser.add_argument('--input', required=True, help='Input probeBase CSV file')
    parser.add_argument('--output', type=str, default='data/databases/open/probeBase_false.csv', help='Output noisy probeBase CSV file')
    parser.add_argument('--mutation-number', type=int, default=5, help='Maximum number of mutations to apply to each probe')
    parser.add_argument('--insertion-rate', type=float, default=0.01, help='Insertion mutation rate')
    parser.add_argument('--deletion-rate', type=float, default=0.01, help='Deletion mutation rate')
    parser.add_argument('--mutation-rate', type=float, default=0.1, help='SNP mutation rate')
    parser.add_argument('--iterations', type=int, default=10, help='Number of iterations to generate noisy data')
    
    args = parser.parse_args()
    
    for i in range(args.mutation_number):
        generate_noisy_probes(
            args.input,
            args.output,
            i,
            insertion_rate=args.insertion_rate,
            deletion_rate=args.deletion_rate,
            mutation_rate=args.mutation_rate,
            iterations=args.iterations
        )

if __name__ == '__main__':
    main() 
    
# Usage:
# python scripts/databases/generate_noisy_probes.py --input data/databases/open/probeBase.csv --output data/databases/open/probeBase_mutated.csv 