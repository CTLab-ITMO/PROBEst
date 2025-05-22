#!/usr/bin/env python3

import pandas as pd
import numpy as np
import random
import re
from typing import Tuple, List, Dict
import argparse

def load_formamide_probabilities(file_path: str) -> Tuple[List[float], List[float]]:
    """Load formamide probabilities from TSV file and normalize them."""
    # Read the file with explicit float conversion
    df = pd.read_csv(file_path, sep=r'\s+', header=None, names=['formamide', 'probability'])
    
    # Convert to float explicitly and handle any potential NaN values
    df['probability'] = pd.to_numeric(df['probability'], errors='coerce')
    df = df.dropna()  # Remove any rows with NaN values
    
    if df.empty:
        raise ValueError("No valid probability values found in the formamide file")
    
    # Convert probabilities to numpy array and normalize
    probs = df['probability'].values.astype(float)
    
    # Check for any invalid values
    if np.any(np.isnan(probs)) or np.any(np.isinf(probs)):
        raise ValueError("Invalid probability values found (NaN or Inf)")
    
    # Normalize probabilities
    total = np.sum(probs)
    if total <= 0:
        raise ValueError("Sum of probabilities must be positive")
    
    probs = probs / total
    
    # Verify the probabilities sum to 1 (within numerical precision)
    if not np.isclose(np.sum(probs), 1.0, rtol=1e-5):
        raise ValueError("Probabilities do not sum to 1 after normalization")
    
    return df['formamide'].tolist(), probs.tolist()

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

def apply_mutations(sequence: str, 
                   insertion_rate: float = 0.01,
                   deletion_rate: float = 0.01,
                   mutation_rate: float = 0.1) -> str:
    """Apply random mutations to the sequence."""
    nucleotides = ['A', 'C', 'G', 'T']
    result = []
    
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
                         formamide_file: str,
                         insertion_rate: float = 0.01,
                         deletion_rate: float = 0.01,
                         mutation_rate: float = 0.1,
                         iterations: int = 1) -> None:
    """Generate noisy probe data with mutations and formamide variations."""
    
    # Load formamide probabilities
    formamide_values, formamide_probs = load_formamide_probabilities(formamide_file)
    
    # Read input probe data
    df = pd.read_csv(input_file, header=None)
    
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
            
            # Apply mutations
            mutated_seq = apply_mutations(clean_seq, 
                                        insertion_rate,
                                        deletion_rate,
                                        mutation_rate)
            
            # Calculate new properties
            new_length = len(mutated_seq)
            new_gc_content = calculate_gc_content(mutated_seq)
            
            # Randomly select new formamide value based on probabilities
            new_formamide = random.choices(formamide_values, 
                                         weights=formamide_probs, 
                                         k=1)[0]
            
            # Create new probe with updated values
            new_probe = probe.copy()
            new_probe['Sequence'] = mutated_seq
            new_probe['Length [nt]'] = str(new_length)
            new_probe['G+C content [%]'] = str(new_gc_content)
            new_probe['Formamide [%]'] = str(new_formamide)
            
            # Generate unique ID for this iteration
            if 'Accession no.' in new_probe:
                new_probe['Accession no.'] = f"{new_probe['Accession no.']}_{iteration + 1}"
            
            noisy_probes.append(new_probe)
        
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
    parser.add_argument('--formamide', type=str, default='data/databases/open/probeBase_formamide.tsv', help='Formamide probabilities TSV file')
    parser.add_argument('--insertion-rate', type=float, default=0.01, help='Insertion mutation rate')
    parser.add_argument('--deletion-rate', type=float, default=0.01, help='Deletion mutation rate')
    parser.add_argument('--mutation-rate', type=float, default=0.1, help='SNP mutation rate')
    parser.add_argument('--iterations', type=int, default=10, help='Number of iterations to generate noisy data')
    
    args = parser.parse_args()
    
    generate_noisy_probes(
        args.input,
        args.output,
        args.formamide,
        args.insertion_rate,
        args.deletion_rate,
        args.mutation_rate,
        args.iterations
    )

if __name__ == '__main__':
    main() 