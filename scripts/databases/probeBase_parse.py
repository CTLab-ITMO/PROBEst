import pandas as pd
from tqdm import tqdm
import re
from typing import List, Dict, Set
import os
import sys
import subprocess
import tempfile
from pathlib import Path

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
            # Run ncbi-genome-download for remaining genera
            genera_str = ','.join(genera_to_download)
            cmd = [
                "ncbi-genome-download",
                "bacteria",
                "-F", "fasta",
                "--genera", genera_str,
                "--assembly-levels", "complete",
                "--refseq-categories", "reference",
                "--output-folder", temp_dir,
                "--parallel", "3",
                "--retries", "3",
                "--verbose"
            ]
            
            print(f"Downloading reference genomes for: {genera_str}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            print(result.stdout)
            
            if result.returncode != 0:
                print(f"Error downloading genomes: {result.stderr}")
                return genome_files
            
            # Find downloaded genome files
            genome_dir = Path(temp_dir) / "refseq" / "bacteria"
            if genome_dir.exists():
                for genus in genera_to_download:
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
                        print(f"Found genome for {genus}")
    
    except Exception as e:
        print(f"Error downloading genomes: {str(e)}")
    
    return genome_files

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
        if not all(base in 'ACGT' for base in probe):
            print(f"Invalid sequence for probe in {genus}: {probe}")
            return None
        
        # Create a temporary file for the probe
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as probe_file:
            # Write the probe sequence in FASTA format
            probe_file.write(f">probe\n{probe}\n")
            probe_path = probe_file.name
        
        # Create a temporary file for BLAST results
        with tempfile.NamedTemporaryFile(delete=False) as blast_file:
            blast_path = blast_file.name
        
        # Run BLAST
        cmd = [
            "blastn",
            "-query", probe_path,
            "-subject", genome_file,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-out", blast_path,
            "-dust", "no",  # Disable dust filtering
            "-word_size", "7",  # Use smaller word size for short sequences
            "-evalue", "1e-5"  # Set a reasonable e-value threshold
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
            return {
                'genus': genus,
                'identity': float(best_hit[2]),
                'length': int(best_hit[3]),
                'mismatches': int(best_hit[4]),
                'evalue': float(best_hit[10]),
                'bitscore': float(best_hit[11])
            }
    
    except Exception as e:
        print(f"Error searching probe in genome for {genus}: {str(e)}")
    
    return None

def process_probeBase(input_file: str, output_file: str) -> None:
    """
    Process probeBase data and perform genome operations.
    
    Args:
        input_file (str): Path to input probeBase CSV file
        output_file (str): Path to save results
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
    
    # Download reference genomes for all genera
    genome_files = download_reference_genomes(all_genera, genome_dir)
    
    if not genome_files:
        print("No genomes were downloaded or found. Exiting.")
        return
    
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
                    'bitscore': blast_result['bitscore']
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
    process_probeBase(input_file, output_file)
