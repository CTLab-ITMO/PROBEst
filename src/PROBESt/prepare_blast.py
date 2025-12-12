"""
Module for preparing BLAST databases from FASTA directories.

This module provides functions to detect whether input paths are directories
containing FASTA files (instead of pre-built BLAST databases) and to convert
them to BLAST databases using the prep_db.sh script.
"""

import os
import subprocess
import glob
from typing import List, Tuple, Optional


# Supported FASTA file extensions
FASTA_EXTENSIONS = {'.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz'}


def is_fasta_directory(path: str) -> bool:
    """
    Check if the given path is a directory containing FASTA files.
    
    Args:
        path: Path to check.
        
    Returns:
        True if path is a directory containing FASTA files, False otherwise.
    """
    if not os.path.isdir(path):
        return False
    
    # Check for FASTA files in the directory
    for ext in FASTA_EXTENSIONS:
        pattern = os.path.join(path, f"*{ext}")
        if glob.glob(pattern):
            return True
    
    return False


def get_fasta_files(directory: str) -> List[str]:
    """
    Get all FASTA files from a directory.
    
    Args:
        directory: Path to the directory containing FASTA files.
        
    Returns:
        List of paths to FASTA files.
    """
    fasta_files = []
    for ext in FASTA_EXTENSIONS:
        pattern = os.path.join(directory, f"*{ext}")
        fasta_files.extend(glob.glob(pattern))
    return sorted(fasta_files)


def prepare_blast_database(
    fasta_dir: str,
    output_db_path: str,
    contig_table_path: str,
    tmp_dir: Optional[str] = None,
    script_path: Optional[str] = None
) -> str:
    """
    Prepare a BLAST database from a directory containing FASTA files.
    
    Uses the prep_db.sh script to merge FASTA files and create a BLAST database.
    
    Args:
        fasta_dir: Path to directory containing FASTA files.
        output_db_path: Path for the output BLAST database.
        contig_table_path: Path for the contig names output file.
        tmp_dir: Optional temporary directory for intermediate files.
        script_path: Path to the directory containing prep_db.sh script.
        
    Returns:
        Path to the created BLAST database.
        
    Raises:
        ValueError: If the directory contains no FASTA files.
        RuntimeError: If the prep_db.sh script fails.
    """
    fasta_files = get_fasta_files(fasta_dir)
    
    if not fasta_files:
        raise ValueError(f"No FASTA files found in directory: {fasta_dir}")
    
    # Determine script path
    if script_path is None:
        script_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
            "scripts", "generator"
        )
    
    prep_db_script = os.path.join(script_path, "prep_db.sh")
    
    if not os.path.exists(prep_db_script):
        raise FileNotFoundError(f"prep_db.sh script not found at: {prep_db_script}")
    
    # Ensure output directories exist
    os.makedirs(os.path.dirname(output_db_path) or '.', exist_ok=True)
    os.makedirs(os.path.dirname(contig_table_path) or '.', exist_ok=True)
    
    # Build command
    cmd = [
        "bash", prep_db_script,
        "-n", output_db_path,
        "-c", contig_table_path,
    ]
    
    if tmp_dir:
        cmd.extend(["-t", tmp_dir])
    
    cmd.extend(fasta_files)
    
    # Run prep_db.sh
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        raise RuntimeError(
            f"prep_db.sh failed with return code {result.returncode}.\n"
            f"stdout: {result.stdout}\n"
            f"stderr: {result.stderr}"
        )
    
    return output_db_path


def prepare_bases_if_needed(
    true_base: str,
    false_bases: List[str],
    output_dir: str,
    contig_table_path: str,
    tmp_dir: Optional[str] = None,
    script_path: Optional[str] = None
) -> Tuple[str, List[str]]:
    """
    Prepare BLAST databases from FASTA directories if needed.
    
    Checks if true_base and false_bases are directories containing FASTA files.
    If they are, creates BLAST databases from them. Otherwise, returns the
    original paths (assumed to be pre-built BLAST databases).
    
    Args:
        true_base: Path to true base (BLAST database or FASTA directory).
        false_bases: List of paths to false bases (BLAST databases or FASTA directories).
        output_dir: Output directory for created BLAST databases.
        contig_table_path: Path for the contig names output file.
        tmp_dir: Optional temporary directory for intermediate files.
        script_path: Path to the directory containing prep_db.sh script.
        
    Returns:
        Tuple of (processed_true_base, processed_false_bases) paths to BLAST databases.
    """
    blast_db_dir = os.path.join(output_dir, ".blast_db")
    os.makedirs(blast_db_dir, exist_ok=True)
    
    # Process true base
    if is_fasta_directory(true_base):
        dir_name = os.path.basename(true_base.rstrip('/'))
        processed_true_base = os.path.join(blast_db_dir, f"true_{dir_name}")
        true_contig_table = contig_table_path if os.path.dirname(contig_table_path) else os.path.join(blast_db_dir, contig_table_path)
        
        print(f"Preparing BLAST database from FASTA directory: {true_base}")
        prepare_blast_database(
            fasta_dir=true_base,
            output_db_path=processed_true_base,
            contig_table_path=true_contig_table,
            tmp_dir=tmp_dir,
            script_path=script_path
        )
        print(f"Created BLAST database: {processed_true_base}")
    else:
        processed_true_base = true_base
    
    # Process false bases
    processed_false_bases = []
    for i, false_base in enumerate(false_bases):
        if is_fasta_directory(false_base):
            dir_name = os.path.basename(false_base.rstrip('/'))
            processed_false_base = os.path.join(blast_db_dir, f"false_{i}_{dir_name}")
            false_contig_table = os.path.join(blast_db_dir, f"contigs_false_{i}.tsv")
            
            print(f"Preparing BLAST database from FASTA directory: {false_base}")
            prepare_blast_database(
                fasta_dir=false_base,
                output_db_path=processed_false_base,
                contig_table_path=false_contig_table,
                tmp_dir=tmp_dir,
                script_path=script_path
            )
            print(f"Created BLAST database: {processed_false_base}")
            processed_false_bases.append(processed_false_base)
        else:
            processed_false_bases.append(false_base)
    
    return processed_true_base, processed_false_bases
