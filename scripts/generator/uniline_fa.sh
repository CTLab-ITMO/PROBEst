#!/bin/bash
# Script: clean_fasta.sh
# Description: This script processes a FASTA file to clean up sequence headers and format the sequences.
#              It removes additional information from headers (e.g., after spaces or pipes) and ensures
#              that sequences are properly formatted (no line breaks within sequences).
#
# Usage: ./clean_fasta.sh -i <input_fasta> -o <output_fasta>
#
# Options:
#   -i <input_fasta>: Path to the input FASTA file (required).
#   -o <output_fasta>: Path to the output FASTA file (required).
#
# Example:
#   ./clean_fasta.sh -i input.fasta -o cleaned.fasta

# Parse command-line options
while getopts i:o: flag
do
    case "${flag}" in
        i) input=${OPTARG};;  # Input FASTA file
        o) output=${OPTARG};; # Output FASTA file
    esac
done

# Process the input FASTA file
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' \
    < $input |\
    sed "s/ .*//g" |\
    sed "s/>.*[|]/>/g" |\
    sed '/^$/d' > $output 