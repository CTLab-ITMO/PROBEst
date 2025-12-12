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


import pandas as pd
from tqdm import tqdm
import re
from typing import List, Dict, Set
import os
import sys
import subprocess
import tempfile
from pathlib import Path
import shutil

# Add the project root directory to Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.PROBESt.genome_operations import genome_parse

def clean_taxon(taxon: str) -> str:
    """Clean a taxonomy entry."""
    # Remove content in parentheses and brackets
    cleaned = re.sub(r'\([^)]*\)', '', taxon)
    cleaned = re.sub(r'\[[^\]]*\]', '', cleaned)
    # Remove special characters and normalize whitespace
    cleaned = re.sub(r'[^a-zA-Z\s]', ' ', cleaned)
    cleaned = re.sub(r'\s+', ' ', cleaned).strip()
    return cleaned

def extract_species_from_taxonomy(taxonomy: str) -> List[str]:
    """
    Extract species names from taxonomy string.
    
    Args:
        taxonomy (str): Semicolon-separated taxonomy string
        
    Returns:
        List[str]: List of valid species names
    """
    if pd.isna(taxonomy):
        return []
    
    # Split taxonomy into levels and clean each level
    taxa = [clean_taxon(t) for t in taxonomy.split(';') if t.strip()]
    
    species_candidates = []
    genera = set()
    
    for taxon in taxa:
        words = taxon.split()
        
        # Case 1: Look for species-level entries (Genus species)
        if len(words) == 2 and words[0][0].isupper() and words[1][0].islower():
            species_candidates.append(f"{words[0]} {words[1]}")
        
        # Case 2: Store potential genus names for fallback
        if len(words) == 1 and words[0][0].isupper():
            genera.add(words[0])
    
    return species_candidates

def check_existing_genomes(genus: str, genome_dir: str) -> str:
    """
    Check if a genome file already exists for a given genus.
    
    Args:
        genus (str): Genus name
        genome_dir (str): Directory containing genome files
        
    Returns:
        str: Path to existing genome file or None if not found
    """
    genome_path = Path(genome_dir) / f"{genus}_reference.fna.gz"
    if genome_path.exists():
        return str(genome_path)
    return None

def download_reference_genomes(genera: Set[str], output_dir: str) -> Dict[str, str]:
    """
    Download reference genomes for given genera using ncbi-genome-download.
    First tries to download reference genomes, then falls back to any available genome.
    
    Args:
        genera (Set[str]): Set of genus names
        output_dir (str): Directory to save the genomes
        
    Returns:
        Dict[str, str]: Mapping of genus to genome file path
    """
    genome_files = {}
    genera_to_download = set()
    
    # Check which genomes already exist
    for genus in genera:
        existing_genome = check_existing_genomes(genus, output_dir)
        if existing_genome:
            genome_files[genus] = existing_genome
            print(f"Using existing genome for {genus}")
        else:
            genera_to_download.add(genus)
    
    if not genera_to_download:
        return genome_files
    
    try:
        # Create a temporary directory for downloads
        with tempfile.TemporaryDirectory() as temp_dir:
            # Try each domain in sequence (using supported group names)
            domains = ['bacteria', 'archaea', 'fungi']
            
            # First try with reference genomes only
            for domain in domains:
                if not genera_to_download:
                    break
                    
                print(f"\nTrying to download reference {domain} genomes for remaining genera: {', '.join(sorted(genera_to_download))}")
                genera_str = ','.join(genera_to_download)
                
                cmd = [
                    "ncbi-genome-download",
                    domain,
                    "-F", "fasta",
                    "--genera", genera_str,
                    "--assembly-levels", "complete",
                    "--refseq-categories", "reference",
                    "--output-folder", temp_dir,
                    "--parallel", "3",
                    "--retries", "3",
                    "--verbose"
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True)
                print(result.stdout)
                
                if result.returncode != 0:
                    print(f"Error downloading {domain} genomes: {result.stderr}")
                    continue
                
                # Find downloaded genome files
                genome_dir = Path(temp_dir) / "refseq" / domain
                if genome_dir.exists():
                    found_genera = set()
                    for genus in list(genera_to_download):  # Create a copy to iterate over
                        # Look for files matching this genus (case-insensitive)
                        genome_files_for_genus = []
                        for path in genome_dir.glob("*/*_genomic.fna.gz"):
                            try:
                                # Read the first few lines of the file to check the organism name
                                with subprocess.Popen(['zcat', str(path)], stdout=subprocess.PIPE) as proc:
                                    header = proc.stdout.readline().decode('utf-8')
                                    if genus.lower() in header.lower():
                                        genome_files_for_genus.append(path)
                            except Exception as e:
                                print(f"Error reading genome file {path}: {str(e)}")
                                continue
                        
                        if genome_files_for_genus:
                            # Take the first genome file for this genus
                            genome_file = genome_files_for_genus[0]
                            output_path = Path(output_dir) / f"{genus}_reference.fna.gz"
                            subprocess.run(["cp", str(genome_file), str(output_path)])
                            genome_files[genus] = str(output_path)
                            found_genera.add(genus)
                            print(f"Found reference {domain} genome for {genus}")
                    
                    # Remove found genera from the download list
                    genera_to_download -= found_genera
            
            # If there are still genera without genomes, try downloading any available genome
            if genera_to_download:
                print("\nNo reference genomes found for some genera. Trying to download any available genome...")
                
                # Clean up the temporary directory
                for item in Path(temp_dir).glob("*"):
                    if item.is_file():
                        item.unlink()
                    elif item.is_dir():
                        shutil.rmtree(item)
                
                for domain in domains:
                    if not genera_to_download:
                        break
                        
                    print(f"\nTrying to download any {domain} genome for remaining genera: {', '.join(sorted(genera_to_download))}")
                    genera_str = ','.join(genera_to_download)
                    
                    cmd = [
                        "ncbi-genome-download",
                        domain,
                        "-F", "fasta",
                        "--genera", genera_str,
                        "--assembly-levels", "complete",
                        "--output-folder", temp_dir,
                        "--parallel", "3",
                        "--retries", "3",
                        "--verbose"
                    ]
                    
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    print(result.stdout)
                    
                    if result.returncode != 0:
                        print(f"Error downloading {domain} genomes: {result.stderr}")
                        continue
                    
                    # Find downloaded genome files
                    genome_dir = Path(temp_dir) / "refseq" / domain
                    if genome_dir.exists():
                        found_genera = set()
                        for genus in list(genera_to_download):
                            genome_files_for_genus = []
                            for path in genome_dir.glob("*/*_genomic.fna.gz"):
                                try:
                                    with subprocess.Popen(['zcat', str(path)], stdout=subprocess.PIPE) as proc:
                                        header = proc.stdout.readline().decode('utf-8')
                                        if genus.lower() in header.lower():
                                            genome_files_for_genus.append(path)
                                except Exception as e:
                                    print(f"Error reading genome file {path}: {str(e)}")
                                    continue
                            
                            if genome_files_for_genus:
                                genome_file = genome_files_for_genus[0]
                                output_path = Path(output_dir) / f"{genus}_reference.fna.gz"
                                subprocess.run(["cp", str(genome_file), str(output_path)])
                                genome_files[genus] = str(output_path)
                                found_genera.add(genus)
                                print(f"Found {domain} genome for {genus}")
                        
                        genera_to_download -= found_genera
    
    except Exception as e:
        print(f"Error downloading genomes: {str(e)}")
    
    if genera_to_download:
        print(f"\nCould not find genomes for the following genera: {', '.join(sorted(genera_to_download))}")
    
    return genome_files

def create_blast_db(genome_file: str) -> str:
    """
    Create a BLAST database for a genome file if it doesn't exist.
    
    Args:
        genome_file (str): Path to the genome file
        
    Returns:
        str: Path to the BLAST database (without extension)
    """
    # Remove .gz extension if present
    base_name = genome_file[:-3] if genome_file.endswith('.gz') else genome_file
    db_path = base_name + '.blastdb'
    
    # Check if database already exists
    if os.path.exists(db_path + '.nhr'):
        return db_path
    
    # Create temporary uncompressed file if input is gzipped
    if genome_file.endswith('.gz'):
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
            subprocess.run(['zcat', genome_file], stdout=temp_file)
            genome_file = temp_file.name
    
    # Create BLAST database
    cmd = [
        "makeblastdb",
        "-in", genome_file,
        "-dbtype", "nucl",
        "-out", db_path
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Clean up temporary file if created
    if genome_file.endswith('.gz'):
        os.unlink(genome_file)
    
    if result.returncode != 0:
        raise Exception(f"Error creating BLAST database: {result.stderr}")
    
    return db_path

def extract_flanking_sequence(genome_file: str, start: int, end: int, strand: str, probe_length: int, flank_size: int = 15) -> Dict[str, str]:
    """
    Extract sequence with flanking regions from genome file.
    
    Args:
        genome_file (str): Path to genome file
        start (int): Start position (1-based)
        end (int): End position (1-based)
        strand (str): Strand direction ('plus' or 'minus')
        probe_length (int): Length of the probe sequence
        flank_size (int): Flanking size (default is 15)
        
    Returns:
        Dict[str, str]: Dictionary containing left_flank, target, and right_flank sequences
    """
    try:
        # If start and end are the same, use probe length to determine end
        if start == end:
            end = start + probe_length - 1
        
        # Adjust positions to include flanking regions
        left_start = max(1, start - flank_size)
        left_end = start - 1
        target_start = start
        target_end = end
        right_start = end + 1
        right_end = end + flank_size
        
        # Extract sequence using zcat and grep
        cmd = f"zcat {genome_file} | grep -v '^>' | tr -d '\\n'"
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise Exception(f"Error extracting sequence: {result.stderr}")
        
        # Get the full sequence
        full_sequence = result.stdout.strip()
        
        # Extract the regions of interest (adjusting for 0-based indexing)
        left_flank = full_sequence[left_start-1:left_end] if left_start <= left_end else ""
        target = full_sequence[target_start-1:target_end]
        right_flank = full_sequence[right_start-1:right_end] if right_start <= right_end else ""
        
        # Reverse complement if on minus strand
        if strand == 'minus':
            left_flank = left_flank.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            target = target.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            right_flank = right_flank.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            # Swap left and right flanks for minus strand
            left_flank, right_flank = right_flank, left_flank
        
        return {
            'left_flank': left_flank,
            'target': target,
            'right_flank': right_flank
        }
    
    except Exception as e:
        print(f"Error extracting flanking sequence: {str(e)}")
        return None

def search_probe_in_genome(probe: str, genome_file: str, genus: str) -> Dict:
    """
    Search for a probe sequence in a genome file.
    
    Args:
        probe (str): Probe sequence to search for
        genome_file (str): Path to genome file
        genus (str): Genus name
        
    Returns:
        Dict: Dictionary with search results
    """
    try:
        # Clean the probe sequence
        probe = probe.strip().upper()
        probe = re.sub(r'[^ACGTUI]', '', probe)
        
        if not all(base in 'ACGTUI' for base in probe):
            raise Warning(f"Invalid sequence for probe in {genus}: {probe}")
            return None
        
        # Create a temporary file for the probe
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as probe_file:
            # Write the probe sequence in FASTA format
            probe_file.write(f">probe\n{probe}\n")
            probe_path = probe_file.name
        
        # Create a temporary file for BLAST results
        with tempfile.NamedTemporaryFile(delete=False) as blast_file:
            blast_path = blast_file.name
        
        # Create BLAST database if it doesn't exist
        try:
            db_path = create_blast_db(genome_file)
        except Exception as e:
            print(f"Error creating BLAST database for {genus}: {str(e)}")
            return None
        
        # Run BLAST
        cmd = [
            "blastn",
            "-query", probe_path,
            "-db", db_path,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand",
            "-out", blast_path, 
            "-word_size", "7",  
            "-evalue", "1e-1"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error running BLAST for {genus}: {result.stderr}")
            return None
        
        # Read BLAST results
        with open(blast_path, 'r') as f:
            hits = f.readlines()
        
        # Clean up temporary files
        os.unlink(probe_path)
        os.unlink(blast_path)
        
        if hits:
            # Parse the best hit
            best_hit = hits[0].strip().split('\t')
            start = int(best_hit[8])  # sstart
            end = int(best_hit[9])    # send
            strand = 'minus' if best_hit[12] == 'minus' else 'plus'
            
            # Extract sequence with flanking regions
            flanking_sequences = extract_flanking_sequence(genome_file, start, end, strand, len(probe))
            
            if flanking_sequences:
                return {
                    'genus': genus,
                    'identity': float(best_hit[2]),
                    'length': int(best_hit[3]),
                    'mismatches': int(best_hit[4]),
                    'evalue': float(best_hit[10]),
                    'bitscore': float(best_hit[11]),
                    'left_flank': flanking_sequences['left_flank'],
                    'target': flanking_sequences['target'],
                    'right_flank': flanking_sequences['right_flank']
                }
    
    except Exception as e:
        print(f"Error searching probe in genome for {genus}: {str(e)}")
    
    return None

def process_probeBase(input_file: str, output_file: str, download_genomes: bool = True) -> None:
    """
    Process probeBase data and perform genome operations.
    
    Args:
        input_file (str): Path to input probeBase CSV file
        output_file (str): Path to save results
        download_genomes (bool): Whether to attempt downloading genomes (default: True)
    """
    print("Reading and processing probeBase data...")
    
    # Create output directory for genomes
    genome_dir = "data/genomes"
    os.makedirs(genome_dir, exist_ok=True)
    
    # Read and pivot the data
    df = pd.read_csv(input_file, header=None)
    df.columns = ['Attribute', 'Value', 'id']
    df_wide = df.pivot(index='id', columns='Attribute', values='Value')
    df_wide.reset_index(inplace=True)
    
    # Extract required columns
    df_filtered = df_wide[['id', 'Sequence', 'Taxonomy']].copy()
    
    # Parse species names from taxonomy
    print("Extracting species names from taxonomy...")
    df_filtered['species_list'] = df_filtered['Taxonomy'].apply(extract_species_from_taxonomy)
    
    # Explode the species list to create one row per species
    df_exploded = df_filtered.explode('species_list')
    df_exploded = df_exploded.rename(columns={'species_list': 'species'})
    
    # Filter out rows with no species
    df_exploded = df_exploded.dropna(subset=['species'])
    
    # Remove duplicates
    df_exploded = df_exploded.drop_duplicates(['id', 'species'])
    
    # Extract all unique genera
    df_exploded['genus'] = df_exploded['species'].apply(lambda x: x.split()[0])
    all_genera = set(df_exploded['genus'].unique())
    print(f"\nFound {len(all_genera)} unique genera")
    
    # Download reference genomes for all genera if requested
    genome_files = {}
    if download_genomes:
        genome_files = download_reference_genomes(all_genera, genome_dir)
        if not genome_files:
            print("No genomes were downloaded or found.")
    else:
        print("Skipping genome download as requested.")
        # Check for existing genomes only
        for genus in all_genera:
            existing_genome = check_existing_genomes(genus, genome_dir)
            if existing_genome:
                genome_files[genus] = existing_genome
                print(f"Using existing genome for {genus}")
    
    # Initialize results list
    results = []
    
    # Process each unique probe
    for probe_id in tqdm(df_exploded['id'].unique(), desc="Processing probes"):
        probe_data = df_exploded[df_exploded['id'] == probe_id].iloc[0]
        species = probe_data['species']
        sequence = probe_data['Sequence']
        genus = species.split()[0]
        
        if genus in genome_files:
            # Search for the probe in the genus-level genome
            blast_result = search_probe_in_genome(sequence, genome_files[genus], genus)
            
            if blast_result:
                results.append({
                    'id': probe_id,
                    'species': species,
                    'genus': genus,
                    'probe_sequence': sequence,
                    'identity': blast_result['identity'],
                    'length': blast_result['length'],
                    'mismatches': blast_result['mismatches'],
                    'evalue': blast_result['evalue'],
                    'bitscore': blast_result['bitscore'],
                    'left_flank': blast_result['left_flank'],
                    'target': blast_result['target'],
                    'right_flank': blast_result['right_flank']
                })
    
    # Create DataFrame from results
    if results:
        df_results = pd.DataFrame(results)
        df_results.to_csv(output_file, index=False)
        print(f"\nResults saved to {output_file}")
        print(f"\nFound matches for {len(df_results)} probes")
    else:
        print("\nNo results were generated")

if __name__ == "__main__":
    input_file = "data/databases/open/probeBase.csv"
    output_file = "data/databases/open/probeBase_genome_results.csv"
    process_probeBase(input_file, output_file, download_genomes=False)