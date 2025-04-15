"""Module for genome operations including fetching, BLAST search and parsing."""

from Bio import Entrez, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import pandas as pd
import os
import tempfile
from typing import List, Tuple

# Set your email for Entrez
Entrez.email = "your.email@example.com"  # TODO: Make this configurable

def genome_fetch(species: str) -> str:
    """
    Fetch genome for a given species from NCBI using Entrez.
    
    Args:
        species (str): Species name to fetch genome for
        
    Returns:
        str: Path to the downloaded genome file
    """
    # Search for the species genome
    handle = Entrez.esearch(db="nucleotide", term=f"{species}[Organism] AND complete genome[Title]")
    record = Entrez.read(handle)
    handle.close()
    
    if not record["IdList"]:
        raise ValueError(f"No genome found for species: {species}")
    
    # Get the first genome record
    genome_id = record["IdList"][0]
    handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
    
    # Save to file
    output_file = f"{species}.fasta"
    with open(output_file, "w") as out_handle:
        out_handle.write(handle.read())
    handle.close()
    
    return output_file

def genome_blastn(genome_file: str, probe: str, extend: int = 0) -> str:
    """
    Perform BLAST search of a probe against a genome.
    
    Args:
        genome_file (str): Path to the genome FASTA file
        probe (str): DNA probe sequence
        extend (int, optional): Number of bases to extend on each side. Defaults to 0.
        
    Returns:
        str: Best hit sequence extended by specified number of bases
    """
    # Create temporary files for BLAST
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as probe_file:
        probe_file.write(f">probe\n{probe}\n")
        probe_path = probe_file.name
    
    # Make BLAST database
    os.system(f"makeblastdb -in {genome_file} -dbtype nucl")
    
    # Run BLAST
    output_file = "blast_results.xml"
    blastn_cline = NcbiblastnCommandline(query=probe_path, 
                                       db=genome_file,
                                       outfmt=5,
                                       out=output_file)
    blastn_cline()
    
    # Parse results
    best_hit = ""
    with open(output_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        first_record = next(blast_records)
        if first_record.alignments:
            alignment = first_record.alignments[0]
            hsp = alignment.hsps[0]
            
            # Get the sequence with extensions
            genome_record = SeqIO.read(genome_file, "fasta")
            start = max(0, hsp.sbjct_start - extend - 1)
            end = min(len(genome_record.seq), hsp.sbjct_end + extend)
            best_hit = str(genome_record.seq[start:end])
    
    # Cleanup
    os.unlink(probe_path)
    os.unlink(output_file)
    for ext in [".nhr", ".nin", ".nsq"]:
        if os.path.exists(genome_file + ext):
            os.unlink(genome_file + ext)
    
    return best_hit

def genome_parse(species: List[str], probe: str, extend: int = 0) -> pd.DataFrame:
    """
    Process multiple species genomes with a probe.
    
    Args:
        species (List[str]): List of species names
        probe (str): DNA probe sequence
        extend (int, optional): Number of bases to extend on each side. Defaults to 0.
        
    Returns:
        pd.DataFrame: DataFrame with results for all species
    """
    results = []
    
    for sp in species:
        try:
            # Fetch genome
            genome_file = genome_fetch(sp)
            
            # Run BLAST
            hit_sequence = genome_blastn(genome_file, probe, extend)
            
            # Store results
            results.append({
                'species': sp,
                'probe_id': f"probe_{len(probe)}bp",
                'species_dna_string': hit_sequence,
                'extend': extend,
                'probe_dna_string': probe
            })
            
            # Cleanup
            os.unlink(genome_file)
            
        except Exception as e:
            print(f"Error processing {sp}: {str(e)}")
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Save to CSV
    output_file = "genome_parse_results.csv"
    df.to_csv(output_file, index=False)
    
    return df 