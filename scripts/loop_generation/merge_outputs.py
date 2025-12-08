#!/usr/bin/env python3
"""
Merge modeling outputs from all species into a single file with species_name column.
"""

import os
import sys
import pandas as pd
import argparse
from pathlib import Path


def find_modeling_files(output_dir):
    """Find all modeling_results.tsv files in species subdirectories."""
    modeling_files = []
    
    for species_dir in Path(output_dir).iterdir():
        if not species_dir.is_dir():
            continue
        
        modeling_file = species_dir / "modeling_results.tsv"
        if modeling_file.exists():
            modeling_files.append((species_dir.name, modeling_file))
    
    return modeling_files


def merge_modeling_outputs(output_dir, output_file):
    """
    Merge all modeling_results.tsv files and add species_name column.
    
    Args:
        output_dir: Directory containing species subdirectories
        output_file: Output file path for merged results
    """
    modeling_files = find_modeling_files(output_dir)
    
    if not modeling_files:
        print(f"WARNING: No modeling_results.tsv files found in {output_dir}")
        return False
    
    print(f"Found {len(modeling_files)} species with modeling results")
    
    all_data = []
    for species_name, modeling_file in modeling_files:
        try:
            df = pd.read_csv(modeling_file, sep='\t')
            df['species_name'] = species_name
            all_data.append(df)
            print(f"  Loaded {len(df)} rows from {species_name}")
        except Exception as e:
            print(f"  WARNING: Failed to load {modeling_file}: {e}")
            continue
    
    if not all_data:
        print("ERROR: No data could be loaded")
        return False
    
    # Merge all dataframes
    merged_df = pd.concat(all_data, ignore_index=True)
    
    # Reorder columns to put species_name first
    cols = ['species_name'] + [c for c in merged_df.columns if c != 'species_name']
    merged_df = merged_df[cols]
    
    # Save merged output
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged_df.to_csv(output_path, sep='\t', index=False)
    
    print(f"\nMerged {len(merged_df)} total rows from {len(all_data)} species")
    print(f"Output saved to: {output_path}")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Merge modeling outputs from all species"
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="test_out/loop_generation",
        help="Input directory containing species subdirectories (default: test_out/loop_generation)"
    )
    parser.add_argument(
        "--output-file",
        type=str,
        default="test_out/loop_generation/output_pulled.tsv",
        help="Output file for merged results (default: test_out/loop_generation/output_pulled.tsv)"
    )
    
    args = parser.parse_args()
    
    # Get project root
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "../.."))
    
    input_dir = os.path.join(project_root, args.input_dir)
    output_file = os.path.join(project_root, args.output_file)
    
    if not os.path.exists(input_dir):
        print(f"ERROR: Input directory not found: {input_dir}")
        sys.exit(1)
    
    print("=" * 60)
    print("Merging modeling outputs")
    print("=" * 60)
    print(f"Input directory: {input_dir}")
    print(f"Output file: {output_file}")
    print()
    
    success = merge_modeling_outputs(input_dir, output_file)
    
    if success:
        print("\n✓ Merge completed successfully")
        sys.exit(0)
    else:
        print("\n✗ Merge failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
