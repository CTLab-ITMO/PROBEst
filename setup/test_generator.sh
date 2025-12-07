#!/bin/bash
# Description: This script prepares BLAST databases and a contigs table for a given set of FASTA files,
#              then executes a Python pipeline for probe generation and analysis. It is designed to handle
#              multiple input datasets (true and false bases) and generate output for further analysis.
#
# Usage: bash setup/test_generator.sh
#
# Workflow:
#   1. Cleans up previous output directories.
#   2. Prepares BLAST databases and contigs tables for:
#      - A "true base" dataset (positive control).
#      - Two "false base" datasets (negative controls).
#   3. Executes a Python pipeline (`pipeline.py`) to analyze the input sequences and generate probes.
#
# This script automatically detects and uses the PROBESt conda environment if available.
# Set PROBEST_ENV_NAME environment variable to use a custom conda environment name (default: probest).

# Get project root (parent directory of setup folder)
SETUP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SETUP_DIR/.." && pwd)"

# Determine which Python command to use
PROBEST_ENV_NAME="${PROBEST_ENV_NAME:-probest}"
PYTHON_CMD="python"

# Check if conda is available and environment exists
if command -v conda &> /dev/null; then
    if conda env list | grep -q "^${PROBEST_ENV_NAME} "; then
        echo "Using conda environment: $PROBEST_ENV_NAME"
        PYTHON_CMD="conda run -n $PROBEST_ENV_NAME python"
    else
        echo "Warning: Conda environment '$PROBEST_ENV_NAME' not found."
        echo "  Using system python. For conda installation, run: bash setup/install.sh"
        echo ""
    fi
else
    echo "Conda not found. Using system python."
    echo "  For conda installation, run: bash setup/install.sh"
    echo ""
fi

# Function to prepare BLAST databases and contigs table
prepare () {
    echo -e "---- Preparation stage ----\n"

    # Clean up previous output directories
    rm -rf "$PROJECT_ROOT/data/test/general/contigs"
    rm -rf "$PROJECT_ROOT/data/test/general/blastn_base"
    rm -rf "$PROJECT_ROOT/data/test/general/output"
    rm -rf "$PROJECT_ROOT/data/test/general/output_oligominer"
    
    # Prepare BLAST database and contigs table for the "true base" dataset
    bash "$PROJECT_ROOT/scripts/generator/prep_db.sh" \
        -n "$PROJECT_ROOT/data/test/general/blastn_base/true_base" \
        -c "$PROJECT_ROOT/data/test/general/contigs" \
        -t "$PROJECT_ROOT/data/test/general/.tmp" \
        "$PROJECT_ROOT/data/test/general/fasta_base/true_base"/*
    
    # Prepare BLAST database and contigs table for the first "false base" dataset
    bash "$PROJECT_ROOT/scripts/generator/prep_db.sh" \
        -n "$PROJECT_ROOT/data/test/general/blastn_base/false_base_1" \
        -c "$PROJECT_ROOT/data/test/general/contigs" \
        -t "$PROJECT_ROOT/data/test/general/tmp" \
        "$PROJECT_ROOT/data/test/general/fasta_base/false_base_1"/*
    
    # Prepare BLAST database and contigs table for the second "false base" dataset
    bash "$PROJECT_ROOT/scripts/generator/prep_db.sh" \
        -n "$PROJECT_ROOT/data/test/general/blastn_base/false_base_2" \
        -c "$PROJECT_ROOT/data/test/general/contigs" \
        -t "$PROJECT_ROOT/data/test/general/tmp" \
        "$PROJECT_ROOT/data/test/general/fasta_base/false_base_2"/*
}

# Change to project root for running pipeline
cd "$PROJECT_ROOT"

# Execute the Python pipeline for probe generation and analysis (using Primer3)

test1(){
    echo "Running pipeline with Primer3..."
    $PYTHON_CMD pipeline.py \
        -i data/test/general/test.fna \
        -o data/test/general/output \
        -tb data/test/general/blastn_base/true_base \
        -fb data/test/general/blastn_base/false_base_1 \
        data/test/general/blastn_base/false_base_2 \
        -c data/test/general/contigs \
        -a FISH \
        --PRIMER_PICK_PRIMER 1 \
        --PRIMER_NUM_RETURN 1
}

# Test with OligoMiner (now obligatory, should be in PROJECT_ROOT/OligoMiner)
test2(){
    OLIGOMINER_DIR=""
    # Check if OLIGOMINER_PATH environment variable is set (from conda environment)
    if [ -n "$OLIGOMINER_PATH" ] && [ -d "$OLIGOMINER_PATH" ]; then
        OLIGOMINER_DIR="$OLIGOMINER_PATH"
    # Check for OligoMiner in project root (standard location)
    elif [ -d "$PROJECT_ROOT/OligoMiner" ]; then
        OLIGOMINER_DIR="$PROJECT_ROOT/OligoMiner"
    fi

    echo "Running pipeline with OligoMiner..."
    $PYTHON_CMD pipeline.py \
        -i data/test/general/test.fna \
        -o data/test/general/output_oligominer \
        -tb data/test/general/blastn_base/true_base \
        -fb data/test/general/blastn_base/false_base_1 \
        data/test/general/blastn_base/false_base_2 \
        -c data/test/general/contigs \
        -a FISH \
        --initial_generator oligominer \
        --oligominer_path "$OLIGOMINER_DIR" \
        --oligominer_probe_length 25 \
        --oligominer_temperature 58
}



# Execute the preparation function
prepare

echo -e "All BLASTn database created" 
#test1
test2