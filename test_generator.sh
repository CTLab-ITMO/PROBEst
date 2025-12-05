#!/bin/bash
# Description: This script prepares BLAST databases and a contigs table for a given set of FASTA files,
#              then executes a Python pipeline for probe generation and analysis. It is designed to handle
#              multiple input datasets (true and false bases) and generate output for further analysis.
#
# Usage: bash test_generator.sh
#
# Workflow:
#   1. Cleans up previous output directories.
#   2. Prepares BLAST databases and contigs tables for:
#      - A "true base" dataset (positive control).
#      - Two "false base" datasets (negative controls).
#   3. Executes a Python pipeline (`pipeline.py`) to analyze the input sequences and generate probes.

# Function to prepare BLAST databases and contigs table
prepare () {
    echo -e "---- Preparation stage ----\n"

    # Clean up previous output directories
    rm -rf data/test/general/contigs
    rm -rf data/test/general/blastn_base
    rm -rf data/test/general/output
    rm -rf data/test/general/output_oligominer
    
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

echo -e "All BLASTn database created" 

# Execute the Python pipeline for probe generation and analysis (using Primer3)
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

# Optional: Test with OligoMiner (uncomment if OligoMiner is installed)
#conda activate probeMining
if [ -n "$OLIGOMINER_PATH" ] && [ -d "$OLIGOMINER_PATH" ]; then
     echo -e "\n---- Testing with OligoMiner ----\n"
    
    # Clean up previous OligoMiner output
    rm -rf data/test/general/output_oligominer
    
    # Execute the Python pipeline with OligoMiner
    python pipeline.py \
        -i data/test/general/test.fna \
        -o data/test/general/output_oligominer \
        -tb data/test/general/blastn_base/true_base \
        -fb data/test/general/blastn_base/false_base_1 \
        data/test/general/blastn_base/false_base_2 \
        -c data/test/general/contigs \
        -a FISH \
        --initial_generator oligominer \
        --oligominer_path "$OLIGOMINER_PATH" \
        --oligominer_probe_length 25 \
        --oligominer_temperature 58
    
    echo -e "\nOligoMiner test completed"
else
    echo -e "\nSkipping OligoMiner test (set OLIGOMINER_PATH environment variable to enable)"
fi