#!/bin/bash
# Script: test_run_generation.sh
# Description: This script prepares BLAST databases and a contigs table for a given set of FASTA files,
#              then executes a Python pipeline for probe generation and analysis. It is designed to handle
#              multiple input datasets (true and false bases) and generate output for further analysis.
#
# Usage: ./test_run_generation.sh
#
# Workflow:
#   1. Cleans up previous output directories.
#   2. Prepares BLAST databases and contigs tables for:
#      - A "true base" dataset (positive control).
#      - Two "false base" datasets (negative controls).
#   3. Executes a Python pipeline (`pipeline.py`) to analyze the input sequences and generate probes.

# Function to prepare BLAST databases and contigs table
prepare () {
    echo "---- Preparation stage ----\n"

    # Clean up previous output directories
    rm -rf data/test/general/contigs
    rm -rf data/test/general/blastn_base
    rm -rf data/test/general/output
    
    # Prepare BLAST database and contigs table for the "true base" dataset
    bash scripts/generator/prep_db.sh \
        -n data/test/general/blastn_base/true_base \
        -c data/test/general/contigs \
        -t data/test/general/.tmp \
        data/test/general/fasta_base/true_base/*
    
    # Prepare BLAST database and contigs table for the first "false base" dataset
    bash scripts/generator/prep_db.sh \
        -n data/test/general/blastn_base/false_base_1 \
        -c data/test/general/contigs \
        -t data/test/general/tmp \
        data/test/general/fasta_base/false_base_1/*
    
    # Prepare BLAST database and contigs table for the second "false base" dataset
    bash scripts/generator/prep_db.sh \
        -n data/test/general/blastn_base/false_base_2 \
        -c data/test/general/contigs \
        -t data/test/general/tmp \
        data/test/general/fasta_base/false_base_2/*
}

# Execute the preparation function
prepare

echo -e "\nAll BLASTn database created" 

# Execute the Python pipeline for probe generation and analysis
python pipeline.py \
    -i data/test/general/test.fna \
    -o data/test/general/output \
    -tb data/test/general/blastn_base/true_base \
    -fb data/test/general/blastn_base/false_base_1 \
    data/test/general/blastn_base/false_base_2 \
    -c data/test/general/contigs \
    -a FISH \
    --PRIMER_PICK_PRIMER 1 \
    --PRIMER_NUM_RETURN 1