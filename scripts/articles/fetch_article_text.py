import os
import sys
import requests
import time
from pathlib import Path
from typing import List, Optional, Tuple
import logging
from tqdm import tqdm
import argparse
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET
import ftplib
import urllib.parse
import magic  # for file type detection

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

def read_ids(input_file: str) -> Tuple[List[str], bool]:
    """Read IDs from input file and determine if they are PMC IDs."""
    with open(input_file, 'r') as f:
        ids = [line.strip() for line in f if line.strip()]
    
    # Check if the first ID is a PMC ID (starts with PMC)
    is_pmc = ids and ids[0].startswith('PMC')
    return ids, is_pmc

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

def is_valid_pdf(file_path: Path) -> bool:
    """Check if a file is a valid PDF."""
    try:
        # Check if file exists and has content
        if not file_path.exists() or file_path.stat().st_size == 0:
            return False
            
        # Check file type using python-magic
        file_type = magic.from_file(str(file_path), mime=True)
        if file_type != 'application/pdf':
            return False
            
        # Check PDF header
        with open(file_path, 'rb') as f:
            header = f.read(5)
            if header != b'%PDF-':
                return False
                
        return True
    except Exception as e:
        logger.error(f"Error validating PDF {file_path}: {e}")
        return False

def download_from_ftp(ftp_url: str, output_path: Path) -> bool:
    """Download file from FTP URL."""
    try:
        # Parse FTP URL
        parsed = urllib.parse.urlparse(ftp_url)
        host = parsed.netloc
        path = parsed.path
        
        # Connect to FTP server
        with ftplib.FTP(host) as ftp:
            ftp.login()  # Anonymous login
            
            # Get the directory and filename
            directory = os.path.dirname(path)
            filename = os.path.basename(path)
            
            # Change to the directory
            ftp.cwd(directory)
            
            # Get file size
            try:
                file_size = ftp.size(filename)
                if file_size == 0:
                    logger.warning(f"FTP file {filename} has zero size")
                    return False
            except:
                pass  # Some FTP servers don't support SIZE command
            
            # Download the file
            with open(output_path, 'wb') as f:
                ftp.retrbinary(f'RETR {filename}', f.write)
            
            # Verify the downloaded file
            if not is_valid_pdf(output_path):
                logger.warning(f"Downloaded file from FTP is not a valid PDF: {output_path}")
                output_path.unlink()  # Delete invalid file
                return False
                
            return True
    except Exception as e:
        logger.error(f"Error downloading from FTP {ftp_url}: {e}")
        if output_path.exists():
            output_path.unlink()  # Clean up failed download
        return False

def download_pdf(url: str, output_path: Path, headers: Optional[dict] = None, max_retries: int = 3) -> bool:
    """Download PDF from URL with retries and validation."""
    if headers is None:
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }
    
    for attempt in range(max_retries):
        try:
            # Handle FTP URLs
            if url.startswith('ftp://'):
                if download_from_ftp(url, output_path):
                    return True
                continue
            
            # Handle HTTP URLs
            response = requests.get(url, stream=True, headers=headers)
            response.raise_for_status()
            
            # Check content type
            content_type = response.headers.get('content-type', '').lower()
            if 'pdf' not in content_type and 'octet-stream' not in content_type:
                logger.warning(f"URL {url} does not point to a PDF (content-type: {content_type})")
                continue
            
            # Check if we got HTML instead of PDF
            content = response.content[:1024]  # Check first 1KB
            if b'<!DOCTYPE html>' in content or b'<html' in content:
                logger.warning(f"Received HTML instead of PDF from {url}")
                continue
            
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            
            # Verify the downloaded file
            if not is_valid_pdf(output_path):
                logger.warning(f"Downloaded file is not a valid PDF: {output_path}")
                output_path.unlink()  # Delete invalid file
                continue
                
            return True
            
        except requests.exceptions.RequestException as e:
            logger.error(f"Error downloading PDF from {url} (attempt {attempt + 1}/{max_retries}): {e}")
            if output_path.exists():
                output_path.unlink()  # Clean up failed download
            if attempt < max_retries - 1:
                time.sleep(1)  # Wait before retry
            continue
    
    return False

def get_pmc_pdf_url(pmc_id: str) -> Optional[str]:
    """Get PDF URL for a PMC article using the OA Web Service API."""
    base_url = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi"
    
    # Remove 'PMC' prefix if present
    pmc_id = pmc_id.replace('PMC', '')
    
    params = {
        'id': f'PMC{pmc_id}',
        'format': 'pdf'
    }
    
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        
        # Parse XML response
        root = ET.fromstring(response.content)
        
        # Look for PDF link in the response
        for record in root.findall('.//record'):
            for link in record.findall('.//link'):
                if link.get('format') == 'pdf':
                    return link.get('href')
        
        # If no PDF found in OA API, try direct PMC URL
        return f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/pdf/"
        
    except (requests.exceptions.RequestException, ET.ParseError) as e:
        logger.error(f"Error fetching PDF URL for PMC{pmc_id}: {e}")
        return None

def fetch_pmc_article(pmc_id: str, outdir: Path, silent: bool) -> None:
    """Fetch article using PMC OA Web Service API."""
    # Create output directory if it doesn't exist
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Remove 'PMC' prefix if present
    pmc_id = pmc_id.replace('PMC', '')
    
    # Check if main article already exists and is valid
    main_pdf_path = outdir / f"PMC{pmc_id}.pdf"
    if main_pdf_path.exists():
        if is_valid_pdf(main_pdf_path):
            if not silent:
                logger.info(f"Skipping PMC{pmc_id} - valid file already exists")
            return
        else:
            if not silent:
                logger.warning(f"Existing file for PMC{pmc_id} is invalid, will redownload")
            main_pdf_path.unlink()  # Delete invalid file
    
    # Get PDF URL using OA Web Service API
    pdf_url = get_pmc_pdf_url(pmc_id)
    
    if pdf_url:
        if download_pdf(pdf_url, main_pdf_path):
            if not silent:
                logger.info(f"Successfully downloaded main article for PMC{pmc_id}")
        else:
            if not silent:
                logger.warning(f"Failed to download PDF for PMC{pmc_id} from {pdf_url}")
    else:
        if not silent:
            logger.warning(f"Could not find PDF URL for PMC{pmc_id}")

def fetch_pubmed_article(pubmed_id: str, outdir: Path, silent: bool) -> None:
    """Fetch article and supplementary materials for a given PubMed ID."""
    # Create output directory if it doesn't exist
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Check if main article already exists
    main_pdf_path = outdir / f"{pubmed_id}.pdf"
    if main_pdf_path.exists():
        if not silent:
            logger.info(f"Skipping {pubmed_id} - file already exists")
        return
    
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

def get_pmc_id_from_pubmed(pubmed_id: str) -> Optional[str]:
    """Convert PubMed ID to PMC ID using Entrez API elink."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    params = {
        'dbfrom': 'pubmed',
        'db': 'pmc',
        'id': pubmed_id,
        'retmode': 'json',
        'tool': 'my_tool',
        'email': 'dvsmutin@gmail.com'
    }
    
    try:
        response = requests.get(f"{base_url}elink.fcgi", params=params)
        response.raise_for_status()
        data = response.json()
        
        # Extract PMC ID from the response
        linksets = data.get('linksets', [])
        if linksets and 'linksetdbs' in linksets[0]:
            for linksetdb in linksets[0]['linksetdbs']:
                if linksetdb.get('dbto') == 'pmc' and 'links' in linksetdb:
                    pmc_ids = linksetdb['links']
                    if pmc_ids:
                        return f"PMC{pmc_ids[0]}"  # Return first PMC ID found
        
        logger.warning(f"No PMC ID found for PubMed ID {pubmed_id}")
        return None
            
    except (requests.exceptions.RequestException, KeyError, ValueError) as e:
        logger.error(f"Error converting PubMed ID {pubmed_id} to PMC ID: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description='Fetch article text and supplementary materials from PubMed or PMC IDs.')
    parser.add_argument('input_file', help='File containing PubMed IDs or PMC IDs (one per line)')
    parser.add_argument('output_directory', help='Directory to save downloaded files')
    parser.add_argument('--silent', action='store_true', help='Run in silent mode with progress bar')
    args = parser.parse_args()

    # Set up logging based on mode
    global logger
    logger = setup_logging(args.silent)
    
    input_file = args.input_file
    outdir = Path(args.output_directory)
    
    # Read IDs and determine if they are PMC IDs
    ids, is_pmc = read_ids(input_file)
    if not args.silent:
        logger.info(f"Found {len(ids)} {'PMC' if is_pmc else 'PubMed'} IDs to process")
    
    # Process each ID
    if args.silent:
        for id in tqdm(ids, desc="Downloading articles"):
            if is_pmc:
                fetch_pmc_article(id, outdir, args.silent)
            else:
                # Convert PubMed ID to PMC ID
                pmc_id = get_pmc_id_from_pubmed(id)
                if pmc_id:
                    fetch_pmc_article(pmc_id, outdir, args.silent)
                else:
                    logger.warning(f"Could not convert PubMed ID {id} to PMC ID")
            time.sleep(.02)  # Be nice to the API
    else:
        for id in ids:
            if is_pmc:
                logger.info(f"Processing PMC ID: {id}")
                fetch_pmc_article(id, outdir, args.silent)
            else:
                logger.info(f"Converting PubMed ID {id} to PMC ID")
                pmc_id = get_pmc_id_from_pubmed(id)
                if pmc_id:
                    logger.info(f"Found PMC ID: {pmc_id}")
                    fetch_pmc_article(pmc_id, outdir, args.silent)
                else:
                    logger.warning(f"Could not convert PubMed ID {id} to PMC ID")
            time.sleep(.02)  # Be nice to the API

if __name__ == "__main__":
    main()

# python scripts/articles/fetch_article_text.py data/articles/pmid-nucleotide-set.txt data/articles2/ --silent
# python scripts/articles/fetch_article_text.py data/articles/pmc_result.txt data/articles2/ --silent