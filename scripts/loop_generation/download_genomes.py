#!/usr/bin/env python3
"""
Download genomes for bacterial species from NCBI.
Downloads up to 100 genomes per species.
Uses NCBI assembly summary to get genome IDs, then downloads via ncbi-genome-download.
"""

import os
import sys
import subprocess
import argparse
import shutil
import tempfile
import urllib.request
import gzip
from pathlib import Path

# List of bacterial species
SPECIES_LIST = [
    "Escherichia coli", "Bacillus subtilis", "Mycoplasma genitalium", 
    "Vibrio fischeri", "Caulobacter crescentus", "Pseudomonas aeruginosa",
    "Salmonella enterica", "Mycobacterium tuberculosis", "Staphylococcus aureus",
    "Streptococcus pneumoniae", "Lactococcus lactis", "Lactobacillus plantarum",
    "Sinorhizobium meliloti", "Agrobacterium tumefaciens", "Helicobacter pylori",
    "Campylobacter jejuni", "Neisseria meningitidis", "Neisseria gonorrhoeae",
    "Clostridium difficile", "Clostridium acetobutylicum", "Corynebacterium glutamicum",
    "Streptomyces coelicolor", "Thermus aquaticus", "Deinococcus radiodurans",
    "Synechocystis sp. PCC 6803", "Anabaena sp. PCC 7120", "Rhodobacter sphaeroides",
    "Shewanella oneidensis", "Geobacter sulfurreducens", "Mycoplasma pneumoniae",
    "Mycobacterium smegmatis", "Pseudomonas putida", "Acinetobacter baylyi",
    "Zymomonas mobilis", "Bdellovibrio bacteriovorus", "Myxococcus xanthus",
    "Bacillus anthracis", "Yersinia pestis", "Listeria monocytogenes",
    "Enterococcus faecalis", "Bacillus cereus", "Vibrio cholerae",
    "Shigella flexneri", "Klebsiella pneumoniae", "Francisella tularensis",
    "Brucella abortus", "Legionella pneumophila", "Rickettsia prowazekii",
    "Chlamydia trachomatis", "Treponema pallidum"
]


def sanitize_species_name(species_name):
    """Convert species name to a safe directory name."""
    return species_name.replace(" ", "_").replace(".", "").replace("/", "_")


def get_assembly_accessions_from_ncbi(species_name, max_genomes=100, assembly_level="all", use_refseq=True):
    """
    Query NCBI assembly summary to get assembly accessions for a species.
    
    Args:
        species_name: Name of the species
        max_genomes: Maximum number of genomes to retrieve
        assembly_level: Assembly level filter (all, complete, chromosome, scaffold, contig)
        use_refseq: If True, use refseq summary (GCF_), else use genbank (GCA_)
    
    Returns:
        List of assembly accessions (e.g., ['GCF_000005825.2', ...] for refseq or ['GCA_000005825.2', ...] for genbank)
    """
    print(f"Querying NCBI for assembly accessions for: {species_name}")
    
    # Download assembly summary file for bacteria
    if use_refseq:
        assembly_summary_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
        accession_prefix = "GCF_"
    else:
        assembly_summary_url = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
        accession_prefix = "GCA_"
    
    try:
        # Download assembly summary
        print(f"Downloading NCBI assembly summary from {'refseq' if use_refseq else 'genbank'}...")
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt') as tmp_file:
            tmp_path = tmp_file.name
        
        urllib.request.urlretrieve(assembly_summary_url, tmp_path)
        
        # Parse assembly summary
        # Format: tab-separated, skip header lines (start with #)
        # Column order: assembly_accession, bioproject, biosample, wgs_master, refseq_category,
        #               taxid, species_taxid, organism_name, infraspecific_name, isolate, version_status,
        #               assembly_level, release_type, genome_rep, seq_rel_date, asm_name, submitter,
        #               gbrs_paired_asm, paired_asm_comp, ftp_path, excluded_from_refseq, relation_to_type_material
        accessions = []
        species_parts = species_name.lower().split()
        
        with open(tmp_path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                # Column 0: assembly_accession
                # Column 7: organism_name
                # Column 11: assembly_level
                accession = parts[0].strip() if len(parts) > 0 else ""
                organism_name = parts[7].strip() if len(parts) > 7 else ""
                asm_level = parts[11].strip().lower() if len(parts) > 11 else ""
                
                # Check if organism name matches species
                if accession and accession.startswith(accession_prefix):
                    organism_lower = organism_name.lower()
                    # Check if all parts of species name are in organism name
                    if all(sp in organism_lower for sp in species_parts):
                        # Check assembly level if specified
                        if assembly_level == "all" or asm_level == assembly_level.lower():
                            accessions.append(accession)
                
                # Limit results early for performance
                if len(accessions) >= max_genomes * 2:  # Get more to filter
                    break
        
        # Clean up
        os.unlink(tmp_path)
        
        # Remove duplicates and limit
        unique_accessions = list(dict.fromkeys(accessions))[:max_genomes]
        
        if unique_accessions:
            print(f"Found {len(unique_accessions)} assembly accessions")
            print(f"Sample accessions: {', '.join(unique_accessions[:3])}")
        else:
            print(f"WARNING: No accessions found for {species_name}")
        
        return unique_accessions
        
    except Exception as e:
        print(f"ERROR: Failed to query NCBI assembly summary: {e}")
        return []


def download_genomes_for_species(species_name, output_dir, max_genomes=100, test_mode=False):
    """
    Download genomes for a single species using ncbi-genome-download.
    
    Args:
        species_name: Name of the species
        output_dir: Base output directory
        max_genomes: Maximum number of genomes to download
        test_mode: If True, only download a few genomes for testing
    """
    sanitized_name = sanitize_species_name(species_name)
    species_dir = os.path.join(output_dir, sanitized_name)
    os.makedirs(species_dir, exist_ok=True)
    
    # Check if we already have enough genomes (count only uncompressed files)
    existing_files = []
    for pattern in ["*.fna", "*.fa", "*.fasta"]:
        existing_files.extend([f for f in Path(species_dir).glob(pattern) if not str(f).endswith('.gz')])
    limit = 3 if test_mode else max_genomes
    if len(existing_files) >= limit:
        print(f"Already have {len(existing_files)} genomes for {species_name}, skipping download")
        return True
    
    print(f"\n{'='*60}")
    print(f"Downloading genomes for: {species_name}")
    print(f"Output directory: {species_dir}")
    print(f"{'='*60}\n")
    
    # Check if ncbi-genome-download is available
    try:
        result = subprocess.run(
            ["ncbi-genome-download", "--version"],
            capture_output=True,
            text=True,
            check=True
        )
        print(f"Using ncbi-genome-download: {result.stdout.strip()}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ERROR: ncbi-genome-download not found!")
        print("Please install it with: pip install ncbi-genome-download")
        return False
    
    # Get assembly accessions from NCBI
    limit = 3 if test_mode else max_genomes
    accessions = get_assembly_accessions_from_ncbi(species_name, max_genomes=limit)
    
    if not accessions:
        print(f"WARNING: No assembly accessions found for {species_name}")
        # Fallback: try using --genera with fuzzy search
        print("Trying fallback method with --genera and --fuzzy-genus...")
        cmd = [
            "ncbi-genome-download",
            "bacteria",
            "--output-folder", output_dir,
            "--genera", species_name,
            "--fuzzy-genus",
            "--format", "fasta",
            "--parallel", "8",
            "--retries", "3",
            "--section", "refseq"
        ]
        acc_file_path = None
    else:
        # Use assembly accessions
        print(f"Downloading {len(accessions)} genomes using assembly accessions...")
        print(f"Accessions: {', '.join(accessions[:5])}{'...' if len(accessions) > 5 else ''}")
        # Create temporary file with accessions (one per line, as required by ncbi-genome-download)
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as acc_file:
            for acc in accessions:
                acc_file.write(acc + '\n')
            acc_file_path = acc_file.name
        
        cmd = [
            "ncbi-genome-download",
            "bacteria",
            "--output-folder", output_dir,
            "--assembly-accessions", acc_file_path,
            "--format", "fasta",
            "--parallel", "8",
            "--retries", "3",
            "--section", "refseq"
        ]
    
    print(f"Running command: {' '.join(cmd[:10])}...")  # Don't print full accessions list
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if result.stdout:
            print("Download output:", result.stdout)
        if result.stderr:
            print("Download warnings/errors:", result.stderr)
        
        # Debug: List what was actually created in output_dir
        print(f"\nChecking output directory structure: {output_dir}")
        if os.path.exists(output_dir):
            for item in os.listdir(output_dir):
                item_path = os.path.join(output_dir, item)
                if os.path.isdir(item_path):
                    print(f"  Found directory: {item}")
                    # List subdirectories
                    try:
                        for subitem in os.listdir(item_path):
                            subitem_path = os.path.join(item_path, subitem)
                            if os.path.isdir(subitem_path):
                                print(f"    - {subitem}/")
                            else:
                                print(f"    - {subitem}")
                    except:
                        pass
        
        # Clean up accession file if created
        if acc_file_path and os.path.exists(acc_file_path):
            try:
                os.unlink(acc_file_path)
            except:
                pass
        
        # Move downloaded genomes to species-specific directory
        # ncbi-genome-download creates structure: output_dir/{refseq,genbank}/bacteria/...
        # Check both refseq and genbank directories
        moved_count = 0
        for section_name in ["refseq", "genbank"]:
            section_dir = os.path.join(output_dir, section_name, "bacteria")
            if not os.path.exists(section_dir):
                continue
            
            print(f"Searching for downloaded files in {section_dir}...")
            
            # Find all FASTA files (when using --assembly-accessions, files are organized by accession)
            fasta_files = []
            for root, dirs, files in os.walk(section_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    # Look for FASTA files (both compressed and uncompressed)
                    if any(file.endswith(ext) for ext in [".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz"]):
                        # Also check if it's actually a FASTA file (not a metadata file)
                        if not any(file.endswith(ext) for ext in [".txt", ".tsv", ".md5", ".report"]):
                            fasta_files.append(file_path)
            
            print(f"Found {len(fasta_files)} FASTA files in {section_dir}")
            
            # Sort and limit
            fasta_files = sorted(fasta_files)[:limit]
            
            # Move files to species directory and decompress if needed
            for fasta_file in fasta_files:
                filename = os.path.basename(fasta_file)
                dest = os.path.join(species_dir, filename)
                if not os.path.exists(dest):
                    try:
                        shutil.move(fasta_file, dest)
                        moved_count += 1
                        print(f"  Moved: {filename}")
                        
                        # Decompress if it's a gzipped file
                        if dest.endswith('.gz'):
                            decompressed_dest = dest[:-3]  # Remove .gz extension
                            if not os.path.exists(decompressed_dest):
                                print(f"  Decompressing: {filename}")
                                with gzip.open(dest, 'rb') as f_in:
                                    with open(decompressed_dest, 'wb') as f_out:
                                        shutil.copyfileobj(f_in, f_out)
                                # Remove the compressed file
                                os.remove(dest)
                                print(f"  Decompressed to: {os.path.basename(decompressed_dest)}")
                    except Exception as e:
                        print(f"  WARNING: Failed to move/decompress {filename}: {e}")
            
            # Clean up the section directory structure after moving files
            if moved_count > 0 and os.path.exists(section_dir):
                try:
                    shutil.rmtree(section_dir)
                    print(f"Cleaned up {section_dir}")
                except Exception as e:
                    print(f"WARNING: Failed to clean up {section_dir}: {e}")
        
        print(f"Moved {moved_count} genome files to {species_dir}")
        
        if moved_count == 0:
            print(f"WARNING: No genomes were moved. Check if species name '{species_name}' matches NCBI taxonomy.")
            return False
        
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Failed to download genomes for {species_name}")
        print(f"Error output: {e.stderr}")
        # Clean up accession file if created
        if acc_file_path and os.path.exists(acc_file_path):
            try:
                os.unlink(acc_file_path)
            except:
                pass
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Download genomes for bacterial species from NCBI"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="test_out/modeling/genomes",
        help="Output directory for downloaded genomes"
    )
    parser.add_argument(
        "--max-genomes",
        type=int,
        default=100,
        help="Maximum number of genomes per species (default: 100)"
    )
    parser.add_argument(
        "--species",
        type=str,
        nargs="+",
        help="Specific species to download (default: all)"
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Test mode: download only 3 genomes per species"
    )
    
    args = parser.parse_args()
    
    # Get project root
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "../.."))
    output_dir = os.path.join(project_root, args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    # Determine which species to download
    species_to_download = args.species if args.species else SPECIES_LIST
    
    if args.test:
        print("TEST MODE: Downloading only 3 genomes per species")
        species_to_download = species_to_download[:3]  # Limit to first 3 for testing
    
    print(f"Will download genomes for {len(species_to_download)} species")
    print(f"Output directory: {output_dir}\n")
    
    # Download genomes for each species
    success_count = 0
    for species in species_to_download:
        if download_genomes_for_species(
            species, 
            output_dir, 
            max_genomes=args.max_genomes,
            test_mode=args.test
        ):
            success_count += 1
    
    print(f"\n{'='*60}")
    print(f"Download complete: {success_count}/{len(species_to_download)} species succeeded")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
