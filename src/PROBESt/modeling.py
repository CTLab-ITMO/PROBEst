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
Modeling module for PROBESt pipeline.

This module performs BLAST mapping of final probes against the input FASTA database
and extracts the best hits with extended sequences (+-10 nucleotides).
"""

import os
import subprocess
import pandas as pd
import tempfile
import warnings
import random
# Suppress Biopython deprecation warning for pairwise2 before importing
warnings.filterwarnings("ignore", category=DeprecationWarning, module="Bio.pairwise2")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import numpy as np
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from typing import Dict, Tuple, Optional
import RNA
from PROBESt.rna_structure import calculate_dimer_G


def to_RNA(dna_string: str) -> str:
    """
    Convert DNA string to RNA by replacing T with U.
    
    Args:
        dna_string: DNA sequence string
        
    Returns:
        RNA sequence string with T replaced by U
    """
    return dna_string.upper().replace('T', 'U')


def replace_degenerate_nucleotides(sequence: str, seq_type: str = "DNA") -> str:
    """
    Replace degenerate/ambiguous nucleotides with valid nucleotides randomly.
    
    IUPAC degenerate nucleotide codes:
    - N -> A, T, G, or C (randomly)
    - R (A/G) -> A or G
    - Y (C/T) -> C or T
    - S (G/C) -> G or C
    - W (A/T) -> A or T
    - K (G/T) -> G or T
    - M (A/C) -> A or C
    - B (C/G/T) -> C, G, or T
    - D (A/G/T) -> A, G, or T
    - H (A/C/T) -> A, C, or T
    - V (A/C/G) -> A, C, or G
    
    Args:
        sequence: Sequence with potentially degenerate nucleotides
        seq_type: "DNA" or "RNA" (for RNA, T is replaced with U)
        
    Returns:
        Sequence with degenerate nucleotides replaced
    """
    # Define degenerate nucleotide mappings
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
    
    result = []
    for char in sequence.upper():
        if char in degenerate_map:
            # Randomly select one of the possible bases
            replacement = random.choice(degenerate_map[char])
            # For RNA, convert T to U
            if seq_type == "RNA" and replacement == 'T':
                replacement = 'U'
            result.append(replacement)
        else:
            # Keep valid nucleotides as-is
            result.append(char)
    
    return ''.join(result)


def reverse_complement_RNA(rna_string: str) -> str:
    """
    Compute reverse complement of RNA sequence.
    
    Args:
        rna_string: RNA sequence string
        
    Returns:
        Reverse complement RNA sequence string
    """
    from Bio.Seq import Seq
    # Convert RNA to DNA temporarily for reverse complement, then back to RNA
    dna_seq = rna_string.replace('U', 'T')
    rev_comp_dna = str(Seq(dna_seq).reverse_complement())
    # Convert back to RNA
    rev_comp_rna = rev_comp_dna.replace('T', 'U')
    return rev_comp_rna


def prepare_blast_db(input_fasta: str, output_db: str) -> str:
    """
    Prepare a BLAST database from input FASTA file.
    
    Args:
        input_fasta: Path to input FASTA file
        output_db: Path for output BLAST database (without extension)
        
    Returns:
        Path to the created BLAST database
    """
    # Create directory for database if needed
    db_dir = os.path.dirname(output_db)
    if db_dir and not os.path.exists(db_dir):
        os.makedirs(db_dir, exist_ok=True)
    
    # Create BLAST database
    cmd = [
        "makeblastdb",
        "-in", input_fasta,
        "-out", output_db,
        "-dbtype", "nucl"
    ]
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )
    
    return output_db


def get_probes_from_output(output_fa: str) -> Dict[str, str]:
    """
    Extract probe sequences from the output FASTA file.
    Converts DNA sequences to RNA when parsing.
    
    Args:
        output_fa: Path to output FASTA file
        
    Returns:
        Dictionary mapping probe IDs to RNA sequences
    """
    probes = {}
    
    if not os.path.exists(output_fa):
        raise FileNotFoundError(f"Output FASTA file not found: {output_fa}")
    
    for record in SeqIO.parse(output_fa, "fasta"):
        # Remove the hit count prefix (H{number}_) if present
        probe_id = record.id
        if probe_id.startswith("H") and "_" in probe_id:
            # Extract the original probe ID after H{number}_
            parts = probe_id.split("_", 1)
            if len(parts) > 1:
                probe_id = parts[1]
        
        probe_seq = str(record.seq).upper()
        # Convert DNA to RNA when parsing
        probes[probe_id] = to_RNA(probe_seq)
    
    return probes


def run_blast_mapping(
    probes: Dict[str, str],
    blast_db: str,
    args,
    output_file: str,
    extend_nuc: int = 10
) -> str:
    """
    Run BLAST mapping of probes against the database.
    Converts RNA sequences back to DNA for BLAST (BLAST expects DNA).
    
    Args:
        probes: Dictionary of probe IDs to RNA sequences
        blast_db: Path to BLAST database
        args: Arguments object with BLAST parameters
        output_file: Path to output BLAST results file
        extend_nuc: Number of nucleotides to extend on each side (default: 10)
        
    Returns:
        Path to the output file
    """
    # Create temporary FASTA file with probes
    # Convert RNA back to DNA for BLAST (BLAST expects DNA sequences)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_fasta:
        for probe_id, seq_rna in probes.items():
            # Convert RNA back to DNA for BLAST
            seq_dna = seq_rna.replace('U', 'T')
            tmp_fasta.write(f">{probe_id}\n{seq_dna}\n")
        tmp_fasta_path = tmp_fasta.name
    
    try:
        # Build BLAST command with extended output format
        # Format 6 with: qseqid, sseqid, pident, length, mismatch, gapopen, 
        # qstart, qend, sstart, send, evalue, bitscore, qseq, sseq
        blastn_cmd = [
            args.blastn,
            "-query", tmp_fasta_path,
            "-db", blast_db,
            "-num_threads", args.threads,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq",
            "-word_size", args.word_size,
            "-reward", args.reward,
            "-penalty", args.penalty,
            "-gapopen", args.gapopen,
            "-gapextend", args.gapextend,
            "-evalue", args.evalue,
            "-out", output_file
        ]
        
        result = subprocess.run(
            blastn_cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
    finally:
        # Clean up temporary file
        if os.path.exists(tmp_fasta_path):
            os.unlink(tmp_fasta_path)
    
    return output_file


def extract_extended_sequence(
    input_fasta: str,
    sseqid: str,
    sstart: int,
    send: int,
    extend: int = 10
) -> Optional[str]:
    """
    Extract sequence from FASTA file with extended flanking regions.
    
    Args:
        input_fasta: Path to input FASTA file
        sseqid: Sequence ID from BLAST result
        sstart: Start position (1-based)
        send: End position (1-based)
        extend: Number of nucleotides to extend on each side
        
    Returns:
        Extended sequence or None if not found
    """
    try:
        # Read the FASTA file
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Check if this is the matching sequence
            # BLAST may return just the ID part before any spaces
            record_id = record.id.split()[0] if ' ' in record.id else record.id
            if record_id == sseqid or record.id == sseqid:
                seq = str(record.seq)
                seq_len = len(seq)
                
                # Determine actual start and end positions (convert to 0-based)
                # Handle reverse complement (if send < sstart)
                if send < sstart:
                    # Reverse complement case
                    actual_start = max(0, send - 1 - extend)
                    actual_end = min(seq_len, sstart + extend)
                    extracted = seq[actual_start:actual_end]
                    # Reverse complement
                    from Bio.Seq import Seq
                    extracted = str(Seq(extracted).reverse_complement())
                else:
                    # Forward strand
                    actual_start = max(0, sstart - 1 - extend)
                    actual_end = min(seq_len, send + extend)
                    extracted = seq[actual_start:actual_end]
                
                return extracted.upper()
        
        return None
        
    except Exception as e:
        print(f"Error extracting extended sequence: {str(e)}")
        return None


def extract_best_hits(
    blast_output: str,
    input_fasta: str,
    extend_nuc: int = 10
) -> pd.DataFrame:
    """
    Extract best hit for each probe from BLAST results.
    
    Args:
        blast_output: Path to BLAST output file (tab-separated)
        input_fasta: Path to input FASTA file for sequence extraction
        extend_nuc: Number of nucleotides to extend on each side
        
    Returns:
        DataFrame with columns: probe_seq, target_seq, target_seq_reverse_complement, 
        vienna_rna_mfe, dna_dna_duplex_dg, dna_rna_duplex_dg, rna_rna_duplex_dg
    """
    # Column names for BLAST output format 6
    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qseq", "sseq"
    ]
    
    # Expected output columns
    output_columns = [
        "probe_id", "probe_seq", "target_seq", "target_seq_reverse_complement",
        "vienna_rna_mfe", "dna_dna_duplex_dg", "dna_rna_duplex_dg", "rna_rna_duplex_dg"
    ]
    
    # Read BLAST results
    if not os.path.exists(blast_output) or os.path.getsize(blast_output) == 0:
        return pd.DataFrame(columns=output_columns)
    
    df = pd.read_csv(blast_output, sep="\t", names=columns, header=None)
    
    if df.empty:
        return pd.DataFrame(columns=output_columns)
    
    # Get best hit for each probe (highest bitscore, then lowest evalue)
    df_sorted = df.sort_values(["qseqid", "bitscore", "evalue"], ascending=[True, False, True])
    best_hits = df_sorted.groupby("qseqid").first().reset_index()
    
    # Extract extended sequences
    results = []
    for _, row in best_hits.iterrows():
        probe_id = row["qseqid"]
        probe_seq = row["qseq"]
        sseqid = row["sseqid"]
        sstart = int(row["sstart"])
        send = int(row["send"])
        
        # Try to get extended sequence from input FASTA
        extended_seq = extract_extended_sequence(
            input_fasta,
            sseqid,
            sstart,
            send,
            extend_nuc
        )
        
        # Use extended sequence if available, otherwise use sseq from BLAST
        if extended_seq:
            target_seq = extended_seq
        else:
            # Fallback to sseq from BLAST (may not have +-10 nucleotides)
            target_seq = row["sseq"]
        
        # Convert only probe sequence to RNA (target sequences remain as DNA)
        probe_seq_rna = to_RNA(probe_seq)
        # Target sequences stay as DNA
        target_seq_dna = target_seq.upper()
        
        # Clean sequences of degenerate nucleotides before energy calculations
        probe_seq_rna_clean = replace_degenerate_nucleotides(probe_seq_rna, seq_type="RNA")
        target_seq_dna_clean = replace_degenerate_nucleotides(target_seq_dna, seq_type="DNA")
        
        # Compute reverse complement of target sequence (DNA)
        from Bio.Seq import Seq
        target_seq_dna_seq = Seq(target_seq_dna_clean)
        target_seq_rev_comp_dna = str(target_seq_dna_seq.reverse_complement())
        
        # Calculate all energies
        # 1. Vienna RNA MFE (RNA single strand energy)
        try:
            md = RNA.md()
            md.noLP = 1
            fc = RNA.fold_compound(probe_seq_rna_clean, md)
            (structure, rna_single_dg) = fc.mfe()
        except Exception as e:
            # Try with original sequence if cleaned one fails
            try:
                fc = RNA.fold_compound(probe_seq_rna, md)
                (structure, rna_single_dg) = fc.mfe()
            except:
                rna_single_dg = 0.0
        
        # 2. DNA-DNA duplex dG (target DNA with its reverse complement)
        try:
            if len(target_seq_dna_clean) > 0 and len(target_seq_rev_comp_dna) > 0:
                dna_dna_duplex_dg = calculate_dimer_G(target_seq_dna_clean, target_seq_rev_comp_dna, type1="DNA", type2="DNA")
            else:
                dna_dna_duplex_dg = 0.0
        except Exception as e:
            # Try with original sequences if cleaned ones fail
            try:
                target_seq_dna_seq_orig = Seq(target_seq_dna)
                target_seq_rev_comp_dna_orig = str(target_seq_dna_seq_orig.reverse_complement())
                if len(target_seq_dna) > 0 and len(target_seq_rev_comp_dna_orig) > 0:
                    dna_dna_duplex_dg = calculate_dimer_G(target_seq_dna, target_seq_rev_comp_dna_orig, type1="DNA", type2="DNA")
                else:
                    dna_dna_duplex_dg = 0.0
            except:
                dna_dna_duplex_dg = 0.0
        
        # 3. DNA-RNA duplex dG (probe RNA with target DNA)
        try:
            if len(probe_seq_rna_clean) > 0 and len(target_seq_dna_clean) > 0:
                dna_rna_duplex_dg = calculate_dimer_G(probe_seq_rna_clean, target_seq_dna_clean, type1="RNA", type2="DNA")
            else:
                dna_rna_duplex_dg = 0.0
        except Exception as e:
            # Try with original sequences if cleaned ones fail
            try:
                if len(probe_seq_rna) > 0 and len(target_seq_dna) > 0:
                    dna_rna_duplex_dg = calculate_dimer_G(probe_seq_rna, target_seq_dna, type1="RNA", type2="DNA")
                else:
                    dna_rna_duplex_dg = 0.0
            except:
                dna_rna_duplex_dg = 0.0
        
        # 4. RNA-RNA duplex dG (probe RNA with its reverse complement - stability)
        try:
            if len(probe_seq_rna_clean) > 0:
                probe_rna_rc = reverse_complement_RNA(probe_seq_rna_clean)
                rna_rna_duplex_dg = calculate_dimer_G(probe_seq_rna_clean, probe_rna_rc, type1="RNA", type2="RNA")
            else:
                rna_rna_duplex_dg = 0.0
        except Exception as e:
            # Try with original sequence if cleaned one fails
            try:
                if len(probe_seq_rna) > 0:
                    probe_rna_rc = reverse_complement_RNA(probe_seq_rna)
                    rna_rna_duplex_dg = calculate_dimer_G(probe_seq_rna, probe_rna_rc, type1="RNA", type2="RNA")
                else:
                    rna_rna_duplex_dg = 0.0
            except:
                rna_rna_duplex_dg = 0.0
        
        results.append({
            "probe_id": probe_id,
            "probe_seq": probe_seq_rna,  # RNA
            "target_seq": target_seq_dna,  # DNA
            "target_seq_reverse_complement": target_seq_rev_comp_dna,  # DNA
            "vienna_rna_mfe": round(rna_single_dg, 2),
            "dna_dna_duplex_dg": round(dna_dna_duplex_dg, 2),
            "dna_rna_duplex_dg": round(dna_rna_duplex_dg, 2),
            "rna_rna_duplex_dg": round(rna_rna_duplex_dg, 2)
        })
    
    return pd.DataFrame(results)


def run_modeling(args, input_fasta: str, output_fa: str, output_table: str, probe_hits: Dict[str, Dict[str, int]] = None) -> str:
    """
    Main function to run the modeling pipeline.
    
    Args:
        args: Arguments object with BLAST parameters
        input_fasta: Path to input FASTA file
        output_fa: Path to output FASTA file with probes
        output_table: Path to output table file
        
    Returns:
        Path to the output table file
    """
    print("\n---- Modeling module ----")
    
    # Step 0: Prepare BLAST database from input.fasta
    print("Preparing BLAST database from input FASTA...")
    db_path = os.path.join(args.output, ".tmp", "modeling_db")
    prepare_blast_db(input_fasta, db_path)
    print(f"BLAST database created: {db_path}")
    
    # Step 1: Get probes from output
    print("Extracting probes from output...")
    probes = get_probes_from_output(output_fa)
    print(f"Found {len(probes)} probes")
    
    # Step 2: Run BLAST mapping
    print("Running BLAST mapping...")
    blast_output = os.path.join(args.output, ".tmp", "modeling_blast.tsv")
    os.makedirs(os.path.dirname(blast_output), exist_ok=True)
    run_blast_mapping(probes, db_path, args, blast_output, extend_nuc=10)
    print("BLAST mapping completed")
    
    # Step 3: Extract best hits
    print("Extracting best hits...")
    results_df = extract_best_hits(blast_output, input_fasta, extend_nuc=10)
    
    # Add per-probe hit counts (true/false databases) if provided
    if probe_hits is None:
        probe_hits = dict()
    results_df["true_db_hits"] = results_df["probe_id"].map(
        lambda pid: probe_hits.get(pid, {}).get("true_hits", 0)
    )
    results_df["false_db_hits"] = results_df["probe_id"].map(
        lambda pid: probe_hits.get(pid, {}).get("false_hits", 0)
    )
    
    # Step 4: Output table
    print(f"Writing results to {output_table}...")
    os.makedirs(os.path.dirname(output_table) if os.path.dirname(output_table) else ".", exist_ok=True)
    results_df.to_csv(output_table, sep="\t", index=False)
    print(f"Results written: {len(results_df)} probe-target pairs")
    
    # Step 5: Create visualizations (if enabled)
    # Check if visualization is enabled (default is True)
    visualize = getattr(args, 'visualize', True)
    if visualize:
        try:
            print("Creating visualizations...")
            vis_dir = create_visualizations_from_table(output_table, args.output)
            print(f"Visualizations created in: {vis_dir}")
        except Exception as e:
            print(f"Warning: Visualization creation failed: {str(e)}")
            print("Continuing without visualizations...")
    else:
        print("Visualizations disabled (--no-visualize flag set)")
    
    print("Modeling module completed\n")
    
    return output_table


def create_visualization(
    probe_seq: str,
    target_seq: str,
    output_path: str
) -> str:
    """
    Create visualization showing:
    1. Sequence alignment with energy calculations
    2. ViennaRNA predicted 2D structure of the probe
    
    Args:
        probe_seq: RNA sequence of the probe
        target_seq: DNA sequence of the target
        output_path: Path to save the visualization PNG
        
    Returns:
        Path to the saved visualization file
    """
    # Ensure probe is RNA (convert if needed)
    probe_rna = probe_seq if 'U' in probe_seq.upper() else to_RNA(probe_seq)
    # Target is DNA, keep as is
    target_dna = target_seq.upper()
    
    # Convert RNA to DNA for alignment (BioPython works better with DNA)
    probe_dna = probe_rna.replace('U', 'T')
    
    # Perform pairwise alignment
    alignments = pairwise2.align.globalxx(probe_dna, target_dna)
    if not alignments:
        # Fallback: use sequences as-is if alignment fails
        aligned_probe = probe_dna
        aligned_target = target_dna
        match_line = ''.join(['|' if p == t else ' ' for p, t in zip(probe_dna[:len(target_dna)], target_dna[:len(probe_dna)])])
    else:
        best_alignment = alignments[0]
        aligned_probe = best_alignment.seqA
        aligned_target = best_alignment.seqB
        match_line = ''.join(['|' if a == b else ' ' for a, b in zip(aligned_probe, aligned_target)])
    
    # Calculate energies
    # Probe is already RNA (probe_rna), target is DNA (target_dna)
    
    # Clean sequences of degenerate nucleotides before energy calculations
    probe_rna_clean = replace_degenerate_nucleotides(probe_rna, seq_type="RNA")
    target_dna_clean = replace_degenerate_nucleotides(target_dna, seq_type="DNA")
    
    # 1. Vienna RNA MFE (RNA single strand energy)
    try:
        md = RNA.md()
        md.noLP = 1  # No lonely pairs
        fc = RNA.fold_compound(probe_rna_clean, md)
        (structure, rna_single_dg) = fc.mfe()
        structure_dot_bracket = structure
    except Exception as e:
        # Try with original sequence if cleaned one fails
        try:
            fc = RNA.fold_compound(probe_rna, md)
            (structure, rna_single_dg) = fc.mfe()
            structure_dot_bracket = structure
        except:
            structure_dot_bracket = '.' * len(probe_rna)
            rna_single_dg = 0.0
    
    # 2. DNA-DNA duplex dG (target DNA with its reverse complement)
    try:
        # Get reverse complement of target DNA
        from Bio.Seq import Seq
        target_dna_seq = Seq(target_dna_clean)
        target_dna_rc = str(target_dna_seq.reverse_complement())
        # Ensure sequences are not empty
        if len(target_dna_clean) > 0 and len(target_dna_rc) > 0:
            # Calculate DNA-DNA duplex energy using calculate_dimer_G
            # This will convert DNA to RNA internally and use RNA.cofold
            dna_dna_duplex_dg = calculate_dimer_G(target_dna_clean, target_dna_rc, type1="DNA", type2="DNA")
        else:
            dna_dna_duplex_dg = 0.0
    except Exception as e:
        dna_dna_duplex_dg = 0.0
        # Try with original sequences if cleaned ones fail
        try:
            target_dna_seq = Seq(target_dna)
            target_dna_rc = str(target_dna_seq.reverse_complement())
            if len(target_dna) > 0 and len(target_dna_rc) > 0:
                dna_dna_duplex_dg = calculate_dimer_G(target_dna, target_dna_rc, type1="DNA", type2="DNA")
        except:
            pass
    
    # 3. DNA-RNA duplex dG (probe RNA with target DNA)
    try:
        # Ensure sequences are not empty
        if len(probe_rna_clean) > 0 and len(target_dna_clean) > 0:
            dna_rna_duplex_dg = calculate_dimer_G(probe_rna_clean, target_dna_clean, type1="RNA", type2="DNA")
        else:
            dna_rna_duplex_dg = 0.0
    except Exception as e:
        dna_rna_duplex_dg = 0.0
        # Try with original sequences if cleaned ones fail
        try:
            if len(probe_rna) > 0 and len(target_dna) > 0:
                dna_rna_duplex_dg = calculate_dimer_G(probe_rna, target_dna, type1="RNA", type2="DNA")
        except:
            pass
    
    # 4. RNA single dG (already calculated above as rna_single_dg, which is the MFE)
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.2, 1], hspace=0.3, wspace=0.3, width_ratios=[1.5, 1])
    
    # Plot 1: Sequence alignment
    ax1 = fig.add_subplot(gs[0, :])
    ax1.axis('off')
    
    # Display alignment
    font_size = 10
    y_pos = 0.85
    
    # Probe sequence
    ax1.text(0.05, y_pos, 'Probe:', fontsize=font_size, fontweight='bold', transform=ax1.transAxes)
    ax1.text(0.15, y_pos, aligned_probe, fontsize=font_size, family='monospace', transform=ax1.transAxes)
    
    # Match line
    y_pos -= 0.15
    ax1.text(0.15, y_pos, match_line, fontsize=font_size, family='monospace', color='green', transform=ax1.transAxes)
    
    # Target sequence
    y_pos -= 0.15
    ax1.text(0.05, y_pos, 'Target:', fontsize=font_size, fontweight='bold', transform=ax1.transAxes)
    ax1.text(0.15, y_pos, aligned_target, fontsize=font_size, family='monospace', transform=ax1.transAxes)
    
    # Energy information (moved lower to avoid intersection)
    y_pos -= 0.35
    energy_text = f'Energy Calculations:\n'
    energy_text += f'  Vienna RNA MFE: {rna_single_dg:.2f} kcal/mol\n'
    energy_text += f'  DNA-DNA Duplex dG: {dna_dna_duplex_dg:.2f} kcal/mol\n'
    energy_text += f'  DNA-RNA Duplex dG: {dna_rna_duplex_dg:.2f} kcal/mol\n'
    energy_text += f'  RNA Single dG: {rna_single_dg:.2f} kcal/mol'
    ax1.text(0.05, y_pos, energy_text, fontsize=font_size, transform=ax1.transAxes,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5), family='monospace')
    
    ax1.set_title('Sequence Alignment and Energy Calculations', fontsize=14, fontweight='bold', pad=20)
    
    # Plot 2: ViennaRNA 2D structure visualization (left, larger)
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.axis('off')
    
    # Display structure information
    structure_text = f'2D Structure (MFE: {rna_single_dg:.2f} kcal/mol)\n\n'
    structure_text += f'Sequence: {probe_rna}\n'
    structure_text += f'Structure: {structure_dot_bracket}\n\n'
    
    # Parse base pairs from dot-bracket notation
    pairs = []
    stack = []
    for i, char in enumerate(structure_dot_bracket):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs.append((j, i))
    
    # Display structure text
    ax2.text(0.05, 0.95, structure_text, fontsize=font_size, family='monospace',
            transform=ax2.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    
    # Create visualization of base pairs as arcs
    if pairs:
        # Create a circular/arc representation
        seq_len = len(probe_rna)
        center_x = 0.55
        center_y = 0.3
        radius = 0.2
        
        # Draw arcs for base pairs
        for start, end in pairs[:min(30, len(pairs))]:  # Limit to 30 pairs for clarity
            # Calculate angles
            angle_start = 2 * 3.14159 * start / seq_len
            angle_end = 2 * 3.14159 * end / seq_len
            
            # Draw arc
            arc = mpatches.Arc((center_x, center_y), radius*2, radius*2,
                              angle=0, theta1=angle_start*180/3.14159, 
                              theta2=angle_end*180/3.14159,
                              color='blue', linewidth=1.5, alpha=0.6,
                              transform=ax2.transAxes)
            ax2.add_patch(arc)
        
        # Draw sequence circle
        circle = plt.Circle((center_x, center_y), radius*0.8, 
                           fill=False, edgecolor='black', linewidth=2,
                           transform=ax2.transAxes)
        ax2.add_patch(circle)
        
        # Add sequence labels around circle (simplified - show every 5th base)
        for i in range(0, seq_len, max(1, seq_len//10)):
            angle = 2 * 3.14159 * i / seq_len
            x = center_x + radius * 0.9 * np.cos(angle)
            y = center_y + radius * 0.9 * np.sin(angle)
            ax2.text(x, y, probe_rna[i], fontsize=7, ha='center', va='center',
                    transform=ax2.transAxes, fontweight='bold')
    
    ax2.set_title('ViennaRNA Predicted 2D Structure', fontsize=12, fontweight='bold', pad=10)
    
    # Plot 3: ViennaRNA structure plot (right, smaller)
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')
    
    # Try to create a simple structure visualization using ViennaRNA's plotting
    try:
        # Create a simple text representation of the structure
        structure_viz_text = f'Structure Visualization\n\n'
        structure_viz_text += f'Sequence length: {len(probe_rna)}\n'
        structure_viz_text += f'Base pairs: {len(pairs)}\n'
        structure_viz_text += f'MFE: {rna_single_dg:.2f} kcal/mol\n\n'
        
        # Show structure in a compact format
        if len(structure_dot_bracket) <= 80:
            structure_viz_text += f'{structure_dot_bracket}'
        else:
            structure_viz_text += f'{structure_dot_bracket[:40]}...\n{structure_dot_bracket[40:80]}...'
        
        ax3.text(0.1, 0.5, structure_viz_text, fontsize=9, family='monospace',
                transform=ax3.transAxes, verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    except Exception as e:
        ax3.text(0.5, 0.5, 'Structure visualization\nnot available', 
                fontsize=10, ha='center', va='center', transform=ax3.transAxes,
                bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.5))
    
    # Save figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_path


def create_visualizations_from_table(
    modeling_table: str,
    output_dir: str
) -> str:
    """
    Create visualizations for all probe-target pairs in the modeling results table.
    
    Args:
        modeling_table: Path to modeling_results.tsv file
        output_dir: Directory to save visualizations
        
    Returns:
        Path to the directory containing visualizations
    """
    # Read the table
    df = pd.read_csv(modeling_table, sep="\t")
    
    # Create output directory
    vis_dir = os.path.join(output_dir, "visualizations")
    os.makedirs(vis_dir, exist_ok=True)
    
    # Create visualization for each probe-target pair
    for idx, row in df.iterrows():
        probe_seq = row["probe_seq"]  # RNA
        target_seq = row["target_seq"]  # DNA
        
        # Create output filename
        output_filename = f"probe_{idx+1}_visualization.png"
        output_path = os.path.join(vis_dir, output_filename)
        
        try:
            create_visualization(probe_seq, target_seq, output_path)
            print(f"Created visualization: {output_filename}")
        except Exception as e:
            print(f"Error creating visualization for probe {idx+1}: {str(e)}")
            continue
    
    return vis_dir

