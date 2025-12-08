#!/bin/bash
# Description: Run probe generation loop for all species
#              For each species, uses its own genomes as true_base and all other species as false_base
#
# Usage: bash scripts/loop_generation/run_generation_loop.sh [species_list]

source ~/.bashrc
conda activate probest
# Get project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Default directories
GENOMES_DIR="$PROJECT_ROOT/test_out/all_genomes"
BLASTDB_DIR="$PROJECT_ROOT/test_out/all_blastdb"
OUTPUT_DIR="$PROJECT_ROOT/test_out/all_output"
CONTIGS_DIR="$PROJECT_ROOT/test_out/modeling/contigs"
cat $PROJECT_ROOT/test_out/modeling/contigs/* > $PROJECT_ROOT/test_out/fullcontig
CONTIG_DIR="$PROJECT_ROOT/test_out/fullcontig"
# Temporary input sampling (user requested path with double 'l')
TMP_INPUT_DIR="$PROJECT_ROOT/test_out/.tmpinput"

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
        echo "  Using system python."
    fi
else
    echo "Conda not found. Using system python."
fi

# Change to project root
cd "$PROJECT_ROOT"

# Get list of species to process
if [ $# -gt 0 ]; then
    # Use provided species list
    SPECIES_LIST=("$@")
else
    # Get all species from genomes directory
    SPECIES_LIST=()
    for species_dir in "$GENOMES_DIR"/*; do
        echo $species_dir
        if [ -d "$species_dir" ]; then
            SPECIES_LIST+=("$(basename "$species_dir")")
        fi
    done
fi

if [ ${#SPECIES_LIST[@]} -eq 0 ]; then
    echo "ERROR: No species found in $GENOMES_DIR"
    #exit 1
fi

echo "=========================================="
echo "Running generation loop for species"
echo "=========================================="
echo "Number of species: ${#SPECIES_LIST[@]}"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Process each species
for species_name in "${SPECIES_LIST[@]}"; do
    echo "----------------------------------------"
    echo "Processing species: $species_name"
    echo "----------------------------------------"
    
    # Check if species directories exist
    species_genomes_dir="$GENOMES_DIR/$species_name"
    species_blastdb="$BLASTDB_DIR/$species_name"
    
    if [ ! -d "$species_genomes_dir" ]; then
        echo "WARNING: Genomes directory not found: $species_genomes_dir, skipping..."
        continue
    fi
    
    if [ ! -f "${species_blastdb}.nhr" ]; then
        echo "WARNING: BLAST database not found: $species_blastdb, skipping..."
        continue
    fi
    
    # Create output directory for this species
    species_output_dir="$OUTPUT_DIR/$species_name"
    mkdir -p "$species_output_dir"
    
    # Build false_base list (all other species)
    false_bases=()
    for other_species in "${SPECIES_LIST[@]}"; do
        if [ "$other_species" != "$species_name" ]; then
            other_blastdb="$BLASTDB_DIR/$other_species"
            if [ -f "${other_blastdb}.nhr" ]; then
                false_bases+=("$other_blastdb")
            fi
        fi
    done
    
    if [ ${#false_bases[@]} -eq 0 ]; then
        echo "WARNING: No false_base databases found, skipping..."
        continue
    fi
    
    echo "True base: $species_blastdb"
    echo "False bases: ${#false_bases[@]} databases"
    echo "Input FASTA: $species_genomes_dir"
    echo "Output: $species_output_dir"
    echo ""
    
    # Prepare temporary sampled input (10 random genomes)
    rm -rf "$TMP_INPUT_DIR"
    mkdir -p "$TMP_INPUT_DIR"
    find "$species_genomes_dir" -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" -o -name "*.fna.gz" -o -name "*.fa.gz" -o -name "*.fasta.gz" \) | shuf | head -n 5 | while read -r f; do
        cp "$f" "$TMP_INPUT_DIR"/
    done
    echo "Sampled inputs placed in: $TMP_INPUT_DIR"
    
    # Run pipeline
    echo "Running pipeline..."
    $PYTHON_CMD pipeline.py \
        -i "$TMP_INPUT_DIR" \
        -o "$species_output_dir" \
        -tb "$species_blastdb" \
        -fb "${false_bases[@]}" \
        -c "$CONTIG_DIR" \
        -a FISH \
        --PRIMER_PICK_PRIMER 10 \
        --PRIMER_NUM_RETURN 10 \
        --top 20 \
        --set_size 100 \
        --iterations 15 \
        --dedegeneration_set_size 25 \
        -t 8

    
    if [ $? -eq 0 ]; then
        echo "✓ Successfully processed $species_name"
    else
        echo "✗ ERROR: Failed to process $species_name"
    fi
    
    echo ""
done

echo "=========================================="
echo "Generation loop complete"
echo "=========================================="
