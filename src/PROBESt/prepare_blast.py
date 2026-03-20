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
Module for preparing BLAST databases from FASTA directories.

This module provides functions to detect whether input paths are directories
containing FASTA files (instead of pre-built BLAST databases) and to convert
them to BLAST databases using the prep_db.sh script.
"""

import os
import subprocess
import glob
from typing import Dict, List, Tuple, Optional


# Supported FASTA file extensions
FASTA_EXTENSIONS = {'.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz'}


def deduplicate_contig_table(contig_table_path: str, contig_id_col: int = 1) -> None:
    """
    Deduplicate a contig table in-place by contig ID (second column by default).

    Lines are tab-separated as produced by ``prep_db.sh``:
    column 0 = genome/FASTA stem, column 1 = contig header id (BLAST ``sseqid``).

    ``prep_db.sh`` appends via ``>>``; repeated runs can duplicate rows. For each
    contig id, the **last** row wins so re-runs stay idempotent.

    Args:
        contig_table_path: Path to the TSV file to read and rewrite.
        contig_id_col: Zero-based index of the contig id column (default ``1``).
    """
    if not os.path.isfile(contig_table_path) or os.path.getsize(contig_table_path) == 0:
        return

    with open(contig_table_path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()

    last_idx_by_key: Dict[str, int] = {}
    line_by_key: Dict[str, str] = {}

    for idx, line in enumerate(lines):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) <= contig_id_col:
            parts = line.split()
        if len(parts) <= contig_id_col:
            key = "__unparsed__%d" % idx
        else:
            key = parts[contig_id_col]

        last_idx_by_key[key] = idx
        line_by_key[key] = line

    ordered_keys = sorted(last_idx_by_key.keys(), key=lambda k: last_idx_by_key[k])
    out_lines = [line_by_key[k] for k in ordered_keys]
    text = ("\n".join(out_lines) + "\n") if out_lines else ""
    with open(contig_table_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)


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

    deduplicate_contig_table(contig_table_path)

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
