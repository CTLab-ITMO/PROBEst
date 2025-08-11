#!/usr/bin/env python3
"""
PDF Nucleotide Sequence Checker

This module provides functionality to check PDF articles for the presence of
nucleotide sequences. It can be used both as an imported module and as a
standalone script.

Usage as module:
    from check_probe_pdf import check_pdf_for_sequences
    
    has_sequences = check_pdf_for_sequences("path/to/article.pdf")

Usage as script:
    python check_probe_pdf.py /path/to/article/folder
"""

from cmd import Cmd
import os
import subprocess
import sys
import re
import argparse
from pathlib import Path
from typing import List, Tuple, Dict, Optional
import logging

# Try to import PDF processing libraries
try:
    import PyPDF2
    PYPDF2_AVAILABLE = True
except ImportError:
    PYPDF2_AVAILABLE = False

try:
    import pdfplumber
    PDFPLUMBER_AVAILABLE = True
except ImportError:
    PDFPLUMBER_AVAILABLE = False

try:
    import fitz  # PyMuPDF
    PYMUPDF_AVAILABLE = True
except ImportError:
    PYMUPDF_AVAILABLE = False

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def extract_text_from_pdf(pdf_path: str) -> str:
    """
    Extract raw text from a PDF file using available libraries.
    
    Args:
        pdf_path: Path to the PDF file
        
    Returns:
        Extracted text as string
        
    Raises:
        RuntimeError: If no PDF processing library is available
        FileNotFoundError: If PDF file doesn't exist
        Exception: For other PDF processing errors
    """
    if not os.path.exists(pdf_path):
        raise FileNotFoundError(f"PDF file not found: {pdf_path}")
    
    # Try PyMuPDF first (most reliable)
    if PYMUPDF_AVAILABLE:
        try:
            doc = fitz.open(pdf_path)
            text = ""
            for page in doc:
                text += page.get_text()
            doc.close()
            return text
        except Exception as e:
            logger.warning(f"PyMuPDF failed for {pdf_path}: {e}")
    
    # Try pdfplumber
    if PDFPLUMBER_AVAILABLE:
        try:
            with pdfplumber.open(pdf_path) as pdf:
                text = ""
                for page in pdf.pages:
                    if page.extract_text():
                        text += page.extract_text()
                return text
        except Exception as e:
            logger.warning(f"pdfplumber failed for {pdf_path}: {e}")
    
    # Try PyPDF2
    if PYPDF2_AVAILABLE:
        try:
            with open(pdf_path, 'rb') as file:
                reader = PyPDF2.PdfReader(file)
                text = ""
                for page in reader.pages:
                    text += page.extract_text() or ""
                return text
        except Exception as e:
            logger.warning(f"PyPDF2 failed for {pdf_path}: {e}")
    
    raise RuntimeError("No PDF processing library available. Please install one of: PyMuPDF, pdfplumber, or PyPDF2")


def clean_text(text: str) -> str:
    """
    Clean extracted text by removing newlines and extra whitespace.
    
    Args:
        text: Raw extracted text
        
    Returns:
        Cleaned text
    """
    # Remove newlines and replace with spaces
    text = re.sub(r'\n+', ' ', text)
    # Remove extra whitespace
    text = re.sub(r'\s+', ' ', text)
    # Strip leading/trailing whitespace
    return text.strip()


def find_nucleotide_sequences(text: str, min_length: int = 7) -> List[str]:
    """
    Find nucleotide sequences in text using regex patterns.
    
    Args:
        text: Text to search in
        min_length: Minimum sequence length (default: 7)
        
    Returns:
        List of found nucleotide sequences
    """
    # Extended nucleotide pattern including degenerate bases
    # Standard nucleotides: A, T, G, C, U
    # Degenerate bases: R (A/G), Y (C/T), S (G/C), W (A/T), K (G/T), M (A/C), B (C/G/T), D (A/G/T), H (A/C/T), V (A/C/G), N (any)
    # RNA nucleotides: I (inosine)
    #nucleotide_pattern = r'[ATGCIURYSWKMBDHVNatgciuryswkmbdhvn]{' + str(min_length) + r',}'
    nucleotide_pattern = '[ATGCUatgcu -]{'+ str(min_length) + '}+'
    
    sequences = re.findall(nucleotide_pattern, text)
    
    return sequences


def check_pdf_for_sequences(pdf_path: str, min_length: int = 10) -> Tuple[bool, List[str]]:
    """
    Check if a PDF contains nucleotide sequences.
    
    Args:
        pdf_path: Path to the PDF file
        min_length: Minimum sequence length to consider
        
    Returns:
        Tuple of (has_sequences, sequences_found)
    """
    try:
        # Extract text from PDF
        raw_text = extract_text_from_pdf(pdf_path)
        
        # Clean the text
        clean_text_content = clean_text(raw_text)
        
        # Find nucleotide sequences
        sequences = find_nucleotide_sequences(clean_text_content, min_length)
        
        has_sequences = len(sequences) > 0
        
        return has_sequences, sequences
        
    except Exception as e:
        logger.error(f"Error processing {pdf_path}: {e}")
        return False, []


def check_folder_for_sequences(folder_path: str, min_length: int = 7, output: str = "pdf_checked") -> Dict[str, bool]:
    """
    Check all PDF files in a folder for nucleotide sequences.
    
    Args:
        folder_path: Path to folder containing PDF files
        min_length: Minimum sequence length to consider
        output: Output file with the space-separated filename and potential number of sequences
        
    Returns:
        Dictionary mapping PDF filenames to boolean indicating presence of sequences
    """
    folder = Path(folder_path)
    
    if not folder.exists():
        raise FileNotFoundError(f"Folder not found: {folder_path}")
    
    if not folder.is_dir():
        raise ValueError(f"Path is not a directory: {folder_path}")
    
    # Find all PDF files
    pdf_files = list(folder.glob("*.pdf"))
    
    if not pdf_files:
        logger.warning(f"No PDF files found in {folder_path}")
        return {}
    
    results = {}
    
    for pdf_file in pdf_files:
        logger.info(f"Processing: {pdf_file.name}")
        has_sequences, sequences = check_pdf_for_sequences(str(pdf_file), min_length)
        results[pdf_file.name] = has_sequences
        
        if has_sequences:
            logger.info(f"  Found {len(sequences)} sequences")
            CMD = "echo '" + pdf_file.name + " " + str(len(sequences)) + "' >> " + output
            subprocess.run(CMD, shell = True)
    
    return results


def main():
    """Main function for command-line execution."""
    parser = argparse.ArgumentParser(
        description="Check PDF articles for nucleotide sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python check_probe_pdf.py /path/to/articles/
  python check_probe_pdf.py /path/to/articles/ --min-length 10 --output pdf_checked.txt
        """
    )
    
    parser.add_argument(
        "folder_path",
        help="Path to folder containing PDF articles"
    )
    
    parser.add_argument(
        "--min-length",
        type=int,
        default=10,
        help="Minimum nucleotide sequence length (default: 7)"
    )
    
    parser.add_argument(
        "--output", "-o",
        default="pdf_checked",
        help="Output file with the pdf files, potentially containing nucleotide probes"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Check if required libraries are available
        if not any([PYPDF2_AVAILABLE, PDFPLUMBER_AVAILABLE, PYMUPDF_AVAILABLE]):
            logger.error("No PDF processing library available!")
            logger.error("Please install one of: PyMuPDF, pdfplumber, or PyPDF2")
            logger.error("Run: pip install PyMuPDF pdfplumber PyPDF2")
            sys.exit(1)
        
        # Process the folder
        results = check_folder_for_sequences(args.folder_path, args.min_length, args.output)
        
        if not results:
            logger.warning("No PDF files found or processed")
            sys.exit(0)
        
        # Output results
        print(f"\nResults for folder: {args.folder_path}")
        print(f"Minimum sequence length: {args.min_length}")
        print("-" * 50)
        
        true_count = 0
        for filename, has_sequences in results.items():
            status = "TRUE" if has_sequences else "FALSE"
            print(f"{filename}: {status}")
            if has_sequences:
                true_count += 1
        
        print("-" * 50)
        print(f"Total PDFs: {len(results)}")
        print(f"PDFs with sequences: {true_count}")
        print(f"PDFs without sequences: {len(results) - true_count}")
        
        # Exit with appropriate code
        if true_count > 0:
            sys.exit(0)  # Success - found sequences
        else:
            sys.exit(1)  # No sequences found
            
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()