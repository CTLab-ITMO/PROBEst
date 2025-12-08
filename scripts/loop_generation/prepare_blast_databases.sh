#!/bin/bash
# Description: Prepare BLAST databases for all downloaded species genomes
#
# Usage: bash scripts/loop_generation/prepare_blast_databases.sh [genomes_dir] [blastdb_dir]

# Get project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Default directories
GENOMES_DIR="${1:-$PROJECT_ROOT/test_out/modeling/genomes}"
BLASTDB_DIR="${2:-$PROJECT_ROOT/test_out/modeling/blastdb}"
CONTIGS_DIR="$PROJECT_ROOT/test_out/modeling/contigs"
TMP_DIR="$PROJECT_ROOT/test_out/modeling/.tmp"

# Create output directories
mkdir -p "$BLASTDB_DIR"
mkdir -p "$CONTIGS_DIR"
mkdir -p "$TMP_DIR"

echo "=========================================="
echo "Preparing BLAST databases for all species"
echo "=========================================="
echo "Genomes directory: $GENOMES_DIR"
echo "BLAST database directory: $BLASTDB_DIR"
echo ""

# Check if prep_db.sh exists
PREP_DB_SCRIPT="$PROJECT_ROOT/scripts/generator/prep_db.sh"
if [ ! -f "$PREP_DB_SCRIPT" ]; then
    echo "ERROR: prep_db.sh not found at $PREP_DB_SCRIPT"
    exit 1
fi

# Process each species directory
species_count=0
for species_dir in "$GENOMES_DIR"/*; do
    if [ ! -d "$species_dir" ]; then
        continue
    fi
    
    species_name=$(basename "$species_dir")
    echo "----------------------------------------"
    echo "Processing species: $species_name"
    echo "----------------------------------------"
    
    # Count FASTA files
    fasta_count=$(find "$species_dir" -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" -o -name "*.fna.gz" -o -name "*.fa.gz" -o -name "*.fasta.gz" \) | wc -l)
    
    if [ "$fasta_count" -eq 0 ]; then
        echo "WARNING: No FASTA files found in $species_dir, skipping..."
        continue
    fi
    
    echo "Found $fasta_count FASTA file(s)"
    
    # Create BLAST database for this species
    species_blastdb="$BLASTDB_DIR/${species_name}"
    species_contigs="$CONTIGS_DIR/${species_name}_contigs.tsv"
    species_tmp="$TMP_DIR/${species_name}_tmp"
    
    # Clean up previous database if exists
    rm -rf "$species_blastdb"*
    rm -f "$species_contigs"
    rm -rf "$species_tmp"
    
    # Prepare BLAST database
    echo "Creating BLAST database..."
    # Find all FASTA files
    fasta_files=$(find "$species_dir" -type f \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" -o -name "*.fna.gz" -o -name "*.fa.gz" -o -name "*.fasta.gz" \))
    
    if [ -z "$fasta_files" ]; then
        echo "WARNING: No FASTA files found, skipping..."
        continue
    fi
    
    bash "$PREP_DB_SCRIPT" \
        -n "$species_blastdb" \
        -c "$species_contigs" \
        -t "$species_tmp" \
        $fasta_files 2>/dev/null
    
    if [ $? -eq 0 ] && [ -f "${species_blastdb}.nhr" ]; then
        echo "✓ BLAST database created successfully: $species_blastdb"
        ((species_count++))
    else
        echo "✗ ERROR: Failed to create BLAST database for $species_name"
    fi
    
    echo ""
done

echo "=========================================="
echo "BLAST database preparation complete"
echo "Successfully processed: $species_count species"
echo "=========================================="
