#!/bin/bash
# Script: fasta2table.sh
# Description: This script converts a FASTA file into a tabular format.
#              Each sequence header and its corresponding sequence are separated by a tab,
#              making it easier to process and analyze sequence data in tools like Excel or R.
#
# Usage: ./fasta2table.sh -i <input_fasta> -o <output_table>
#
# Options:
#   -i <input_fasta>: Path to the input FASTA file (required).
#   -o <output_table>: Path to the output tabular file (required).
#
# Example:
#   ./fasta2table.sh -i input.fasta -o output_table.txt

# Parse command-line options
while getopts i:o: flag
do
    case "${flag}" in
        i) input=${OPTARG};;  # Input FASTA file
        o) output=${OPTARG};; # Output tabular file
    esac
done

# Convert FASTA to tabular format

awk 'BEGIN{RS=">"}{print ">"$1"\t"$2;}' $input |\
    sort |\
    sed 's/>\t.*//g' |\
    sed 's/\t/\n/g' |\
    awk '!/^[[:space:]]*$/' > $output