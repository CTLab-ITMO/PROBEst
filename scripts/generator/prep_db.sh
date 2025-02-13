#!/bin/bash
# Script: merge_fasta_and_create_blastdb.sh
# Description: This script merges multiple FASTA files into a single file, processes sequence headers,
#              and creates a BLAST database from the merged sequences. It supports gzipped FASTA files
#              and allows customization of the output database name and temporary directory.
#
# Usage: ./merge_fasta_and_create_blastdb.sh -n <database_name> -c <contig_name> -t <tmp_dir> <fasta_files>
#
# Options:
#   -n <database_name>: Name of the output BLAST database (required).
#   -c <contig_name>:   Output file to store contig names and their corresponding sequence headers (required).
#   -t <tmp_dir>:      Temporary directory for intermediate files (optional, defaults to ./.tmp).
#   <fasta_files>:     List of input FASTA files (gzipped or uncompressed).
#
# Example:
#   ./merge_fasta_and_create_blastdb.sh -n my_database -c contig_names.txt -t ./temp file1.fa.gz file2.fasta

# Parse command-line options
while getopts n:c:t: flag
do
    case "${flag}" in
        n) database_name=${OPTARG};;  # Name of the BLAST database
        c) contig_name=${OPTARG};;   # File to store contig names and sequence headers
        t) tmp=${OPTARG};;           # Temporary directory (optional)
    esac
done
shift $(( OPTIND - 1 ))  # Shift arguments to process input files

# Set default temporary directory if not provided
if [ -z $tmp ]; then
    tmp=./.tmp
fi
mkdir -p $tmp  # Create temporary directory if it doesn't exist

# Create directory for the BLAST database if it doesn't exist
dbdir=$(echo $database_name | sed "s/\/.*/\//g")
mkdir -p $dbdir

# Process each input FASTA file
for fname in "$@"; do
    # Extract the base name of the file (without path and extension)
    clear_name=$(echo $fname | sed "s/.*\///g" | sed "s/\..*//g")
    echo "Processing $clear_name"

    # Decompress gzipped files (if applicable)
    if [[ $fname == *.gz ]]; then
        zcat $fname |\
            sed 's/ .*//g' >> $tmp/$clear_name.fa
    else
        sed 's/ .*//g' $fname >> $tmp/$clear_name.fa
    fi

    # Extract sequence headers and prepend the clear_name for identification
    grep "^>" $tmp/$clear_name.fa |\
    sed "s/^>/${clear_name}\t/" >> $contig_name
done

# Merge all processed FASTA files into a single file
cat $tmp/*.fa 1>> $tmp/merged.fa 2>/dev/null
cat $tmp/*.fasta 1>> $tmp/merged.fa 2>/dev/null

# Create a BLAST database from the merged FASTA file
makeblastdb -in $tmp/merged.fa -out $database_name -dbtype nucl > /dev/null

# Clean up temporary files and directory
rm -r $tmp

echo "BLAST database created successfully: $database_name/n"