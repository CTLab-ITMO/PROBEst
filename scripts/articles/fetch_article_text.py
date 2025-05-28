#!/usr/bin/env python3

import os
import sys
import requests
import time
from pathlib import Path
from typing import List, Optional
import logging
from tqdm import tqdm
import argparse

def setup_logging(silent: bool) -> logging.Logger:
    """Set up logging based on mode."""
    if silent:
        # Disable all logging output
        logging.basicConfig(level=logging.CRITICAL)
        return logging.getLogger(__name__)
    else:
        # Set up normal logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)

def read_pubmed_ids(input_file: str) -> List[str]:
    """Read PubMed IDs from input file."""
    with open(input_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def get_article_info(pubmed_id: str) -> Optional[dict]:
    """Get article information from PubMed API."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    params = {
        'db': 'pubmed',
        'id': pubmed_id,
        'retmode': 'json',
        'rettype': 'abstract',
        'tool': 'my_tool',
        'email': 'your_email@example.com'  # Replace with your email
    }
    
    try:
        response = requests.get(f"{base_url}esummary.fcgi", params=params)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        logger.error(f"Error fetching article info for {pubmed_id}: {e}")
        return None

def download_pdf(url: str, output_path: Path) -> bool:
    """Download PDF from URL."""
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        return True
    except requests.exceptions.RequestException as e:
        logger.error(f"Error downloading PDF from {url}: {e}")
        return False

def fetch_article(pubmed_id: str, outdir: Path, silent: bool) -> None:
    """Fetch article and supplementary materials for a given PubMed ID."""
    # Create output directory if it doesn't exist
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Get article information
    article_info = get_article_info(pubmed_id)
    if not article_info:
        return
    
    try:
        result = article_info['result'][pubmed_id]
        
        # Extract DOI from articleids
        doi = None
        for article_id in result.get('articleids', []):
            if article_id.get('idtype') == 'doi':
                doi = article_id.get('value')
                break
        
        if doi:
            # Try to get full text from various sources
            pdf_urls = [
                f"https://doi.org/{doi}",  # Direct DOI link
                f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{result.get('pmcid', '')}/pdf/",  # PMC link
            ]
            if not silent:
                print(f"Attempting to download from URLs: {pdf_urls}")
            
            # Try to download main article
            main_pdf_path = outdir / f"{pubmed_id}.pdf"
            for url in pdf_urls:
                if download_pdf(url, main_pdf_path):
                    if not silent:
                        logger.info(f"Successfully downloaded main article for {pubmed_id}")
                    break
            
            # Try to download supplementary materials
            suppl_pdf_path = outdir / f"{pubmed_id}_suppl.pdf"
            suppl_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{result.get('pmcid', '')}/bin/supp_{pubmed_id}.pdf"
            if download_pdf(suppl_url, suppl_pdf_path):
                if not silent:
                    logger.info(f"Successfully downloaded supplementary materials for {pubmed_id}")
        else:
            if not silent:
                logger.warning(f"No DOI found for PubMed ID {pubmed_id}")
        
    except (KeyError, TypeError) as e:
        if not silent:
            logger.error(f"Error processing article info for {pubmed_id}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Fetch article text and supplementary materials from PubMed IDs.')
    parser.add_argument('input_file', help='File containing PubMed IDs (one per line)')
    parser.add_argument('output_directory', help='Directory to save downloaded files')
    parser.add_argument('--silent', action='store_true', help='Run in silent mode with progress bar')
    args = parser.parse_args()

    # Set up logging based on mode
    global logger
    logger = setup_logging(args.silent)
    
    input_file = args.input_file
    outdir = Path(args.output_directory)
    
    # Read PubMed IDs
    pubmed_ids = read_pubmed_ids(input_file)
    if not args.silent:
        logger.info(f"Found {len(pubmed_ids)} PubMed IDs to process")
    
    # Process each PubMed ID
    if args.silent:
        for pubmed_id in tqdm(pubmed_ids, desc="Downloading articles"):
            fetch_article(pubmed_id, outdir, args.silent)
            time.sleep(.1)
    else:
        for pubmed_id in pubmed_ids:
            logger.info(f"Processing PubMed ID: {pubmed_id}")
            fetch_article(pubmed_id, outdir, args.silent)
            time.sleep(.1)

if __name__ == "__main__":
    main()

# python scripts/articles/fetch_article_text.py data/articles/pmid-nucleotide-set.txt data/articles/ --silent